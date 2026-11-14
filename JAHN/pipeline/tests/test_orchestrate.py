"""The orchestrator: policy allowlist + rules, agent fallback, probe parsing, and the round driver.

The end-to-end group uses the toy plugin (no abguard, no data); the `needs_data` group exercises the
real de novo VHH generator (ANDERS' tool) end to end.
"""
import json
import sys
import textwrap
from datetime import datetime, timedelta
from pathlib import Path

import polars as pl
import pytest

from abpipe import orchestrate, policy as P, probe as PR
from abpipe.policy import AgentPolicy, PolicyError, apply_loop_patch, apply_module_patch, rules_decision, validate_patch

# --------------------------------------------------------------------- toy plugin
PLUGIN = textwrap.dedent("""
    import polars as pl
    from abpipe.base import Filter, Generator, Scorer
    from abpipe.registry import register

    @register("generator", "t2_gen")
    class TGen(Generator):
        \"\"\"toy sequences of increasing length\"\"\"
        def generate(self, n, seed, ctx):
            return pl.DataFrame([dict(sequence_id=f"t2_{seed}_{i}", kind="vhh", heavy_sequence="E" * (10 + i),
                                      light_sequence="", design_description="toy") for i in range(n)])

    @register("filter", "t2_min_len")
    class TMinLen(Filter):
        def __init__(self, min_len=12):
            self.min_len = min_len
        def check(self, df, ctx):
            return df.select("sequence_id", reason=pl.when(pl.col("heavy_sequence").str.len_chars() < self.min_len)
                             .then(pl.lit("too short")).otherwise(None))

    @register("scorer", "t2_len")
    class TLen(Scorer):
        outputs = {"score__len": "higher"}
        def score(self, df, ctx):
            return df.select("sequence_id", score__len=pl.col("heavy_sequence").str.len_chars().cast(pl.Float64))
""")


@pytest.fixture(scope="module")
def plugin(tmp_path_factory):
    d = tmp_path_factory.mktemp("plugin2")
    (d / "abpipe_orch_test_plugin.py").write_text(PLUGIN)
    sys.path.insert(0, str(d))
    yield "abpipe_orch_test_plugin"
    sys.path.remove(str(d))


def toy_cfg(tmp_path, plugin, *, n=8, rounds=None):
    raw = {"run": {"name": "t", "out_dir": str(tmp_path / "runs"), "seed": 1},
           "data_dir": str(tmp_path), "repo_root": str(tmp_path), "plugins": [plugin],
           "generators": [{"name": "t2_gen", "n": n}],
           "filters": [{"name": "t2_min_len", "params": {"min_len": 12}}],
           "scorers": [{"name": "t2_len"}],
           "aggregator": {"name": "borda", "params": {"columns": ["score__len"]}}}
    if rounds:
        raw["orchestrate"] = rounds
    return raw


BASE_STATE = {"loop": {"select": 60, "min_pair_hamming": 3, "max_sub_share": 0.25, "pool_frac": 0.10},
              "effective": {"scored_evolution": {"n": 20000, "min_subs": 4, "max_subs": 8, "eta": 0.5,
                                                "temperature": 1.0, "frac_sample": 0.4},
                            "prefilter_gate": {"min_p_binder": 0.8}},
              "rounds": [], "patience": 2, "primary": "score__pred_log10_kd", "primary_direction": "lower",
              "deadline_remaining_s": None, "est_round_seconds": None,
              "probe": {"gpu": None, "containers": [], "heavy_images": [], "heavy_job_running": False,
                        "disk_free_gb": 100.0, "load1": 1.0, "japanfold_balance_usd": None}}


def state(**kw):
    s = json.loads(json.dumps(BASE_STATE))          # deep copy
    s.update(kw)
    return s


# --------------------------------------------------------------------- patch allowlist
def test_validate_patch_accepts_and_rejects():
    validate_patch({"loop": {"select": 80}, "modules": {"scored_evolution": {"n": 30000, "eta": 0.3}}})
    for bad in ({"loop": {"select": 5}}, {"loop": {"select": "x"}}, {"loop": {"bogus": 1}},
                {"modules": {"scored_evolution": {"min_subs": 3}}},            # the 70%-identity floor
                {"modules": {"scored_evolution": {"eta": 1.5}}},
                {"modules": {"window_novelty": {"min_dist": 2}}},               # the >2-mutations rule
                {"modules": {"nope": {"x": 1}}},
                {"modules": {"prefilter_gate": {"min_p_binder": 0.2}}},
                {"what": 1}):
        with pytest.raises(PolicyError):
            validate_patch(bad)


def test_apply_module_patch_moves_n_and_params_and_keeps_consistency():
    raw = {"generators": [{"name": "scored_evolution", "n": 100, "params": {"eta": 0.5}}],
           "filters": [{"name": "prefilter_gate", "params": {"min_p_binder": 0.8}}]}
    out = apply_module_patch(raw, {"modules": {"scored_evolution": {"n": 30000, "eta": 0.3},
                                               "prefilter_gate": {"min_p_binder": 0.7}}})
    assert raw["generators"][0]["n"] == 100 and raw["generators"][0]["params"]["eta"] == 0.5   # untouched
    assert out["generators"][0]["n"] == 30000 and out["generators"][0]["params"]["eta"] == 0.3
    assert out["filters"][0]["params"]["min_p_binder"] == 0.7
    with pytest.raises(PolicyError, match="does not have"):
        apply_module_patch(raw, {"modules": {"vhh_variants": {"n": 40}}})
    with pytest.raises(PolicyError, match="min_subs"):
        apply_module_patch({"generators": [{"name": "scored_evolution", "n": 10,
                                            "params": {"min_subs": 4, "max_subs": 8}}]},
                           {"modules": {"scored_evolution": {"max_subs": 3}}})


def test_apply_loop_patch():
    lp = apply_loop_patch({"select": 60, "min_pair_hamming": 3, "max_sub_share": 0.25}, {"loop": {"select": 80}})
    assert lp["select"] == 80 and lp["min_pair_hamming"] == 3
    with pytest.raises(PolicyError):
        apply_loop_patch({}, {"loop": {"min_pair_hamming": 2}})                  # the diversity floor


# --------------------------------------------------------------------- rules policy
def test_cold_start_and_deadline_stops():
    d = rules_decision(state())
    assert not d.stop and d.patch == {"loop": {}, "modules": {}}
    assert rules_decision(state(deadline_remaining_s=30)).stop
    assert rules_decision(state(deadline_remaining_s=120, est_round_seconds=300)).stop
    assert rules_decision(state(probe={**BASE_STATE["probe"], "disk_free_gb": 1.0})).stop


def test_starved_funnel_raises_the_batch_and_relaxes_the_gate():
    hist = [{"round": 0, "survivors": 0, "new_in_selected": 0, "selected_medians": {"score__pred_log10_kd": 1.0}}]
    d = rules_decision(state(rounds=hist))
    assert d.patch["modules"]["scored_evolution"]["n"] == 30000
    assert d.patch["modules"]["prefilter_gate"]["min_p_binder"] == 0.75
    assert any("funnel is starved" in r for r in d.reasoning)


def test_plateau_explores_then_stops_when_capped():
    flat = {"round": 0, "survivors": 500, "new_in_selected": 0,
            "selected_medians": {"score__pred_log10_kd": 0.5}, "mean_pairwise_hamming": 5.0}
    hist = [dict(flat, round=r) for r in range(3)]
    d = rules_decision(state(rounds=hist))
    assert d.patch["modules"]["scored_evolution"] == pytest.approx({"temperature": 1.25, "eta": 0.35, "max_subs": 9})
    capped = state(rounds=hist, effective={"scored_evolution": {"n": 20000, "min_subs": 4, "max_subs": 10,
                                                               "eta": 0.1, "temperature": 3.0, "frac_sample": 0.4},
                                           "prefilter_gate": {"min_p_binder": 0.8}})
    d2 = rules_decision(capped)
    assert d2.stop and d2.stop_reason == "converged"


def test_improvement_or_new_blood_keeps_the_loop_going():
    hist = [{"round": r, "survivors": 500, "new_in_selected": 0 if r < 2 else 3,
             "selected_medians": {"score__pred_log10_kd": 1.0 - 0.1 * r}} for r in range(3)]
    d = rules_decision(state(rounds=hist))
    assert not d.stop and "scored_evolution" not in d.patch["modules"]


def test_diversity_collapse_pulls_towards_fresh_samples():
    hist = [{"round": 0, "survivors": 500, "new_in_selected": 5,
             "selected_medians": {"score__pred_log10_kd": 0.5}, "mean_pairwise_hamming": 3.2}]
    d = rules_decision(state(rounds=hist))
    assert d.patch["modules"]["scored_evolution"]["frac_sample"] == 0.55
    assert d.patch["modules"]["scored_evolution"]["eta"] == 0.4


def test_short_selection_widens_the_pool_walk():
    hist = [{"round": 0, "survivors": 600, "selected": 41, "new_in_selected": 5,
             "selected_medians": {"score__pred_log10_kd": 0.5}, "mean_pairwise_hamming": 6.0}]
    d = rules_decision(state(rounds=hist))
    assert d.patch["loop"]["pool_frac"] == 0.15
    assert any("came back short" in r for r in d.reasoning)


# --------------------------------------------------------------------- agent policy
def agent_script(tmp_path, body):
    p = tmp_path / "agent.sh"
    p.write_text(f"#!/bin/bash\ncat > /dev/null\n{body}\n")
    p.chmod(0o755)
    return f"bash {p}"


def test_agent_patch_is_applied_when_valid(tmp_path):
    cmd = agent_script(tmp_path, """echo '{"reasoning":["batch up"],"patch":{"modules":{"scored_evolution":{"n":123}}},"stop":false}'""")
    d = AgentPolicy(agent_cmd=cmd, timeout_s=60)(state())
    assert d.source == "agent" and d.patch["modules"]["scored_evolution"]["n"] == 123
    assert d.reasoning == ["batch up"] and not d.stop


def test_agent_falls_back_to_rules_on_garbage_bounds_or_timeouts(tmp_path):
    for body in ("echo 'I think we should stop now, no JSON here'",
                 """echo '{"reasoning":["x"],"patch":{"modules":{"scored_evolution":{"min_subs":2}}}}'""",
                 "echo '{}'"):
        d = AgentPolicy(agent_cmd=agent_script(tmp_path, body), timeout_s=60)(state(rounds=[
            {"round": 0, "survivors": 0, "new_in_selected": 0, "selected_medians": {"score__pred_log10_kd": 1.0}}]))
        assert d.source == "rules(agent-failed)" and not d.stop
        assert d.patch["modules"]["scored_evolution"]["n"] == 30000        # the rules fallback did the adapting
    d = AgentPolicy(agent_cmd="sleep 5", timeout_s=1)(state())
    assert d.source == "rules(agent-failed)"


def test_extract_json_finds_the_object_in_prose():
    assert P._extract_json("thinking... {\"a\": 1} done") == {"a": 1}


# --------------------------------------------------------------------- probe
def test_parse_gpu_and_docker():
    assert PR.parse_gpu("37, 4300, 79600\n") == {"util_pct": 37, "mem_used_mb": 4300, "mem_total_mb": 79600}
    assert PR.parse_gpu("") is None and PR.parse_gpu(None) is None
    d = PR.parse_docker_ps("quay.io/apheris/oslo-hackathon:apheris-boltz-v1.3.0-boltz2\nnginx\n")
    assert d["heavy_job_running"] and d["heavy_images"] == ["quay.io/apheris/oslo-hackathon:apheris-boltz-v1.3.0-boltz2"]
    assert not PR.parse_docker_ps("nginx")["heavy_job_running"]


def test_parse_deadline():
    assert orchestrate.parse_deadline(None) is None
    assert orchestrate.parse_deadline("17:30").strftime("%H:%M") == "17:30"
    assert orchestrate.parse_deadline("2026-11-14T17:30").hour == 17
    with pytest.raises(ValueError):
        orchestrate.parse_deadline("sometime")


# --------------------------------------------------------------------- generic selection
def test_select_generic_diverse_keeps_a_levenshtein_floor():
    ranked = pl.DataFrame({"sequence_id": [f"c{i}" for i in range(6)],
                           "heavy_sequence": ["AAAA", "AAAB", "AABB", "ABBB", "BBBB", "CCCA"],
                           "score__x": [6, 5, 4, 3, 2, 1]}).with_row_index("rank", offset=1)
    sel = orchestrate.select_generic_diverse(ranked, 6, min_dist=3)
    seqs = sel["heavy_sequence"].to_list()
    from rapidfuzz.distance import Levenshtein
    assert all(Levenshtein.distance(a, b) >= 3 for i, a in enumerate(seqs) for b in seqs[i + 1:])
    assert sel["heavy_sequence"].to_list()[0] == "AAAA"          # best-first


# --------------------------------------------------------------------- round driver, end to end
def test_orchestrate_runs_rounds_and_journals_every_decision(tmp_path, plugin):
    out = tmp_path / "orch"
    doc = orchestrate.run_orchestrated(toy_cfg(tmp_path, plugin), rounds=3, select=10, out_dir=out,
                                       patience=2, log=lambda m: None)
    assert [r["round"] for r in doc["rounds"]] == [0, 1, 2] and "stopped" not in doc
    for r in range(3):
        rd = out / f"round_{r:02d}"
        assert (rd / "selected.parquet").exists() and (rd / "selected.csv").exists()
        assert (rd / "loop_round.json").exists() and (rd / "manifest.json").exists()
    lines = (out / "decisions.jsonl").read_text().strip().splitlines()
    assert len(lines) == 3 and all(json.loads(l)["reasoning"] for l in lines)
    md = (out / "DECISIONS.md").read_text()
    assert md.count("## Round") == 3 and md.count("- why:") >= 3
    assert all("score__len" in json.loads(l).get("selected_medians", {}) or
               json.loads(l).get("stop") for l in lines)
    s0 = pl.read_parquet(out / "round_00" / "selected.parquet")
    s2 = pl.read_parquet(out / "round_02" / "selected.parquet")
    assert s2["score__len"].max() >= s0["score__len"].max()        # elitism: the pool always contains round 0's best
    assert doc["primary"] == "score__len" and doc["primary_direction"] == "higher"


def test_orchestrate_resumes_finished_rounds(tmp_path, plugin):
    out = tmp_path / "orch"
    orchestrate.run_orchestrated(toy_cfg(tmp_path, plugin), rounds=2, select=10, out_dir=out, log=lambda m: None)
    doc = orchestrate.run_orchestrated(toy_cfg(tmp_path, plugin), rounds=4, select=10, out_dir=out, log=lambda m: None)
    assert [r["round"] for r in doc["rounds"]] == [0, 1, 2, 3]
    assert len((out / "decisions.jsonl").read_text().strip().splitlines()) == 4


def test_orchestrate_stops_on_a_tight_deadline(tmp_path, plugin):
    out = tmp_path / "orch"
    doc = orchestrate.run_orchestrated(toy_cfg(tmp_path, plugin), rounds=3, select=10, out_dir=out,
                                       deadline=datetime.now() + timedelta(seconds=30), log=lambda m: None)
    assert doc.get("stopped") and "deadline" in doc["stopped"]
    assert not (out / "round_00").exists() or not (out / "round_00" / "selected.parquet").exists()


def test_orchestrate_dry_run_executes_nothing(tmp_path, plugin):
    out = tmp_path / "orch"
    doc = orchestrate.run_orchestrated(toy_cfg(tmp_path, plugin), rounds=3, select=10, out_dir=out,
                                       dry_run=True, log=lambda m: None)
    assert "stopped" not in doc and not (out / "round_00").exists()
    assert (out / "DECISIONS.md").exists()                        # the header, but no round sections
    assert "## Round" not in (out / "DECISIONS.md").read_text()


# --------------------------------------------------------------------- the de novo generator
def test_vhh_generator_needs_a_parent_list():
    from abpipe.generators.vhh_variants import VhhVariants
    with pytest.raises(ValueError, match="ranked"):
        VhhVariants()
    VhhVariants(ranked="some.csv")


@pytest.mark.needs_data
def test_vhh_generator_end_to_end_from_the_demo_parents(data_dir, repo_root, tmp_path):
    from abpipe.base import Context
    from abpipe.generators.vhh_variants import VhhVariants
    ctx = Context(data_dir=data_dir, repo_root=repo_root, run_dir=tmp_path, cache_dir=tmp_path, seed=0)
    g = VhhVariants(ranked=str(repo_root / "ANDERS/demo/ranked.csv"))
    df = g.generate(8, seed=3, ctx=ctx)
    assert 0 < df.height <= 8
    assert set(df["kind"]) == {"vhh"} and df["sequence_id"].n_unique() == df.height
    assert df["heavy_sequence"].str.len_chars().max() <= 250
    assert {"info__strategy", "info__parent_score", "info__identity"} <= set(df.columns)
    assert all(r <= 0.62 for r in df["info__identity"].to_list())
