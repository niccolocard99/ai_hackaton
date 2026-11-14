"""Pre-filter scorer, novelty index, scored_evolution generator and the round driver.

The first group needs neither abguard nor the organiser data (stubbed model and tables). The `needs_data` group runs
against the real pre-filter and checks the wrapper against Niccolo's own published numbers.
"""
import json
import re

import numpy as np
import polars as pl
import pytest

from abpipe import _niccolo, data, novelty
from abpipe.base import Context, Parent
from abpipe.generators.scored_evolution import ScoredEvolution, file_sha256, history_digest
from abpipe.iterate import select_diverse

WT = "WGGDGFYAMD"
HEAVY = "Q" * 96 + "SR" + WT + "YWG" + "Q" * 9            # window at 98..107, flanks SR / YWG, 120 residues
PARENT = Parent(id="tz", heavy=HEAVY, light="D" * 107, cdrh3="SRWGGDGFYAMDY", window=WT, window_start=98)
SEQUON = re.compile(r"(?=N[^P][ST])")


def _nm(w):
    return sum(a != b for a, b in zip(w, WT))


def _wins(df):
    return [h[98:108] for h in df["heavy_sequence"].to_list()]


@pytest.fixture
def ctx(tmp_path):
    return Context(data_dir=tmp_path, repo_root=tmp_path, run_dir=tmp_path, cache_dir=tmp_path, seed=0, parent=PARENT)


@pytest.fixture
def stub(monkeypatch):
    """A fitted-model stand-in with random weights, a toy FACS table, and a small 'provided windows' index."""
    rng = np.random.default_rng(1)
    fit = _niccolo.Fitted(pf=None, mod=None, tag="stub", n_facs_dropped=0, clf_w=rng.normal(size=(10, 20)),
                          rank_w=rng.normal(size=(10, 20)), wt_window=WT)
    near = []
    for _ in range(3000):                                  # provided windows: 1-3 substitutions from wild type
        w = list(WT)
        for p in rng.choice(10, size=rng.integers(1, 4), replace=False):
            w[p] = _niccolo.AA[rng.integers(0, 20)]
        near.append("".join(w))
    facs = pl.DataFrame({"cdrh3": near, "binding_class": rng.choice(["high", "medium", "low"], size=len(near)).tolist(),
                         "edit_distance": [_nm(w) for w in near]})
    monkeypatch.setattr(_niccolo, "get_prefilter", lambda **kw: fit)
    monkeypatch.setattr(data, "facs_table", lambda c: facs)
    monkeypatch.setattr(novelty, "provided_index", lambda c: novelty.WindowIndex(near))
    return SimpleStub(fit, near)


class SimpleStub:
    def __init__(self, fit, provided):
        self.fit, self.provided = fit, provided


# ----------------------------------------------------------------------------- novelty index
def test_encode_decode_roundtrip_and_validation():
    rng = np.random.default_rng(0)
    ws = ["".join(r) for r in rng.choice(list(_niccolo.AA), size=(50, 10))]
    assert novelty.decode(novelty.encode(ws)) == ws
    with pytest.raises(ValueError, match="10 residues"):
        novelty.encode(["AAA"])
    with pytest.raises(ValueError, match="non-standard"):
        novelty.encode(["AAAAAAAAAX"])


def test_window_index_matches_brute_force_hamming():
    rng = np.random.default_rng(0)
    refs = ["".join(r) for r in rng.choice(list("ACDE"), size=(400, 10))]    # small alphabet -> many near neighbours
    R = novelty.encode(refs)
    Q = R[rng.integers(0, len(R), 600)].copy()
    for i in range(len(Q)):
        for p in rng.choice(10, size=rng.integers(0, 5), replace=False):
            Q[i, p] = rng.integers(0, 4)
    brute = np.minimum((Q[:, None, :] != R[None, :, :]).sum(2).min(1), 3)
    idx = novelty.WindowIndex(refs)
    assert (idx.distance_class(Q) == brute).all()
    assert (idx.within(Q, 2) == (brute <= 2)).all() and (idx.within(Q, 1) == (brute <= 1)).all()
    assert set(brute.tolist()) == {0, 1, 2, 3}, "the test must exercise every distance class"


# ----------------------------------------------------------------------------- selection
def _ranked(windows):
    rows = [{"sequence_id": f"c{i}", "heavy_sequence": HEAVY[:98] + w + HEAVY[108:]} for i, w in enumerate(windows)]
    return pl.DataFrame(rows).with_row_index("rank", offset=1)


def test_select_diverse_enforces_pairwise_distance_and_substitution_cap():
    rng = np.random.default_rng(2)
    windows = []
    for _ in range(400):
        w = list(WT)
        for p in rng.choice(10, size=5, replace=False):
            w[p] = rng.choice(list("ADEKRST"))
        windows.append("".join(w))
    sel = select_diverse(_ranked(windows), PARENT, 25, min_pair_hamming=3, max_sub_share=0.25, pool_frac=1.0)
    wins = _wins(sel)
    assert 5 <= len(wins) <= 25 and sel["sequence_id"][0] == "c0"
    assert sel["rank"].to_list() == sorted(sel["rank"].to_list())
    assert min(sum(a != b for a, b in zip(x, y)) for i, x in enumerate(wins) for y in wins[i + 1:]) >= 3
    from collections import Counter
    cnt = Counter((p, w[p]) for w in wins for p in range(10) if w[p] != WT[p])
    assert max(cnt.values()) <= max(2, int(np.ceil(0.25 * len(wins)))), "the cap is a share of the RETURNED list"


def test_select_diverse_never_pads_the_list_from_the_middle_of_the_pool():
    """Regression: with a cap that grew with the list, a pool whose strong candidates share a few favoured residues
    was walked down to its median (53 of 60 picks came from pool ranks ~3250 of 6556). The top slice here is one tight
    cluster, so at most one of it can be kept; the 450 diverse mid-pool candidates must NOT be used to fill the list."""
    rng = np.random.default_rng(3)
    base = list(WT)
    for p, a in ((1, "K"), (2, "S"), (3, "E"), (7, "E"), (9, "I")):
        base[p] = a
    cluster = []
    for _ in range(50):                                            # each differs from the base by <= 1 edit -> pairwise <= 2
        w = list(base)
        w[rng.integers(0, 10)] = rng.choice(list("ADEKRST"))
        cluster.append("".join(w))
    rest = []
    for _ in range(450):
        w = list(WT)
        for p in rng.choice(10, size=5, replace=False):
            w[p] = rng.choice(list("ACDEFGHIKLMNPQRSTVWY"))
        rest.append("".join(w))
    sel = select_diverse(_ranked(cluster + rest), PARENT, 25, min_pair_hamming=3, pool_frac=0.10)
    assert sel.height == 1 and sel["rank"][0] == 1, "only the best of the tight top cluster may be kept"
    padded = select_diverse(_ranked(cluster + rest), PARENT, 25, min_pair_hamming=3, pool_frac=1.0)
    assert padded.height > 10 and padded["rank"].max() > 50, "sanity: with the pool limit off, the walk does reach the mid-pool"


def test_select_diverse_cap_is_a_share_of_the_list_that_comes_back():
    """Regression: with the cap tied to the requested size (25% of 60 = 15), 15 of a 17-long list shared A106E. The
    pool here has 40 strong candidates that all share one substitution and are otherwise mutually far apart, and
    nothing else good: the list must come back short, with that substitution in at most a quarter of it."""
    rng = np.random.default_rng(5)
    shared = []
    for i in range(40):                                  # same G103K, other positions differ in every member
        w = list(WT)
        w[3] = "K"
        pos = [p for p in range(10) if p != 3]
        for p in rng.choice(pos, size=3, replace=False):
            w[p] = rng.choice(list("ADEFHILMNQRSTVWY"))
        shared.append("".join(w))
    sel = select_diverse(_ranked(shared), PARENT, 60, min_pair_hamming=3, max_sub_share=0.25, pool_frac=1.0)
    wins = _wins(sel)
    n_shared = sum(w[3] == "K" for w in wins)
    assert 1 <= len(wins) < 20, "a homogeneous pool must give a short list"
    assert n_shared <= max(2, int(np.ceil(0.25 * len(wins)))), (n_shared, len(wins))


# ----------------------------------------------------------------------------- generator
def test_cold_start_obeys_every_rule(stub, ctx):
    df = ScoredEvolution(min_subs=4, max_subs=7).generate(500, seed=3, ctx=ctx)
    wins = _wins(df)
    assert len(df) == len(set(wins)) == 500 and df["sequence_id"].n_unique() == 500
    assert 4 <= min(map(_nm, wins)) and max(map(_nm, wins)) <= 7
    assert not any("C" in w for w in wins)
    assert (novelty.WindowIndex(stub.provided).distance_class(novelty.encode(wins)) >= 3).all()
    assert not any(SEQUON.search("SR" + w + "YWG") for w in wins)
    assert all(h[:98] == HEAVY[:98] and h[108:] == HEAVY[108:] for h in df["heavy_sequence"])
    assert set(df["light_sequence"]) == {PARENT.light} and set(df["origin"]) == {"scored_evolution:sample"}


def test_same_seed_same_output_different_seed_differs(stub, ctx):
    g = ScoredEvolution()
    a, b, c = g.generate(200, 5, ctx), g.generate(200, 5, ctx), g.generate(200, 6, ctx)
    assert _wins(a) == _wins(b) and _wins(a) != _wins(c)


def test_prior_follows_the_classifier_weights(stub, ctx):
    w = np.zeros((10, 20)); w[3, _niccolo.AA.index("W")] = 10.0              # WT window has D at index 3
    stub.fit.clf_w[:] = w
    df = ScoredEvolution(w_clf=1.0, w_rank=0.0, w_facs=0.0, temperature=0.5).generate(400, 1, ctx)
    at3 = [x[3] for x in _wins(df) if x[3] != WT[3]]
    assert len(at3) > 50 and np.mean([a == "W" for a in at3]) > 0.9


def test_no_new_sequons_even_when_the_prior_asks_for_them(stub, ctx):
    stub.fit.clf_w[:] = 0.0
    for a in "NST":
        stub.fit.clf_w[:, _niccolo.AA.index(a)] = 5.0
    df = ScoredEvolution(w_clf=1.0, w_rank=0.0, w_facs=0.0, temperature=0.5).generate(300, 2, ctx)
    wins = _wins(df)
    assert len(wins) > 100 and np.mean([sum(c in "NST" for c in w) for w in wins]) > 3     # the pull towards N/S/T is real
    assert not any(SEQUON.search("SR" + w + "YWG") for w in wins)


def test_feedback_history_and_integrity(stub, ctx, tmp_path):
    r0 = tmp_path / "r0.parquet"
    df0 = ScoredEvolution().generate(300, 1, ctx)
    df0.write_parquet(r0)
    sel = tmp_path / "sel.parquet"
    df0.head(20).with_row_index("rank", offset=1).write_parquet(sel)
    elites = _wins(df0.head(20))
    g1 = ScoredEvolution(parents=str(sel), parents_sha256=file_sha256(sel), history=[str(r0)],
                         history_sha256=history_digest([r0]), n_elite=20)
    df1 = g1.generate(400, 2, ctx)
    assert not set(_wins(df1)) & set(_wins(df0)), "history windows must never be proposed again"
    ops = set(df1["origin"].to_list())
    assert {"scored_evolution:mutate", "scored_evolution:cross", "scored_evolution:sample"} <= ops
    for w in _wins(df1.filter(pl.col("origin") == "scored_evolution:mutate")):
        assert min(sum(a != b for a, b in zip(w, e)) for e in elites) <= 2
    sel.write_bytes(sel.read_bytes() + b"x")                                  # the parents file changes under the config
    with pytest.raises(ValueError, match="sha256"):
        g1.generate(10, 3, ctx)


def test_bad_parameters_are_rejected():
    for kw in ({"min_subs": 5, "max_subs": 4}, {"eta": 1.5}, {"temperature": 0}, {"min_dist_provided": 5},
               {"frac_sample": 0, "frac_mutate": 0, "frac_cross": 0}):
        with pytest.raises(ValueError):
            ScoredEvolution(**kw)


# ----------------------------------------------------------------------------- real pre-filter
def _fit(data_dir, repo_root, tmp_path):
    try:
        return _niccolo.get_prefilter(data_dir=data_dir, cache_dir=tmp_path, repo_root=repo_root)
    except _niccolo.PrefilterUnavailable as e:
        pytest.skip(str(e))


@pytest.mark.needs_data
def test_wrapper_reproduces_the_published_validation(data_dir, repo_root, tmp_path):
    """NICCOLO/eda/AGENT_GUIDE.md section 1: classifier AUROC 0.927 (0.913 at >=4 substitutions, n=556),
    ranker Spearman +0.687 (+0.628 at >=4, n=131). If these move, the model changed: update the guide first."""
    fit = _fit(data_dir, repo_root, tmp_path)
    spr, facs = fit.mod.load()
    v = fit.mod.validate(spr, facs)
    assert v["clf"]["n"] == 620 and v["clf"]["pos"] == 179 and v["clf"]["n_ge4"] == 556
    assert v["clf"]["auroc"] == pytest.approx(0.927, abs=0.003)
    assert v["clf"]["auroc_ge4"] == pytest.approx(0.913, abs=0.003)
    assert v["rank"]["n"] == 179 and v["rank"]["n_ge4"] == 131
    assert v["rank"]["spearman"] == pytest.approx(0.687, abs=0.005)
    assert v["rank"]["spearman_ge4"] == pytest.approx(0.628, abs=0.005)


@pytest.mark.needs_data
def test_scorer_equals_a_direct_prefilter_call(data_dir, repo_root, tmp_path):
    from abpipe.scorers.niccolo_prefilter import NiccoloPrefilter
    fit = _fit(data_dir, repo_root, tmp_path)
    parent = data.load_parent({"table": _niccolo.SPR_REL, "window_column": "cdrh3", "id": "trastuzumab"}, data_dir)
    c = Context(data_dir=data_dir, repo_root=repo_root, run_dir=tmp_path, cache_dir=tmp_path, seed=0, parent=parent)
    rng = np.random.default_rng(4)
    wins = []
    for _ in range(40):
        w = list(parent.window)
        for p in rng.choice(10, size=rng.integers(4, 8), replace=False):
            w[p] = rng.choice(list("ADEFGHIKLMNPQRSTVWY"))
        wins.append("".join(w))
    s = parent.window_start
    heavy = [parent.heavy[:s] + w + parent.heavy[s + 10:] for w in wins]
    heavy += [parent.heavy[:s] + "WGGDGFYAMX" + parent.heavy[s + 10:], parent.heavy[:-3]]      # unscorable rows
    df = pl.DataFrame({"sequence_id": [f"s{i}" for i in range(len(heavy))], "heavy_sequence": heavy})
    out = NiccoloPrefilter().score(df, c)
    direct = fit.pf.score(wins)
    np.testing.assert_allclose(out["score__p_binder"][:40].to_numpy(), direct["p_binder"].to_numpy(), rtol=0, atol=1e-12)
    np.testing.assert_allclose(out["score__pred_log10_kd"][:40].to_numpy(), direct["pred_log10_kd"].to_numpy(), rtol=0, atol=1e-12)
    assert out["score__p_binder"][40:].null_count() == 2 and out["score__pred_log10_kd"][40:].null_count() == 2
    assert set(out.columns) == {"sequence_id", "score__p_binder", "score__pred_log10_kd"}      # FACS stays out unless asked


@pytest.mark.needs_data
def test_two_rounds_end_to_end_learn_from_each_other(data_dir, repo_root, tmp_path):
    """Real filters (incl. the official-checker wrapper) and the real pre-filter, two rounds, small n."""
    from abpipe import iterate
    try:
        _niccolo.find_splits()
    except _niccolo.PrefilterUnavailable as e:
        pytest.skip(str(e))
    cfg = {"run": {"name": "t", "out_dir": str(tmp_path / "runs"), "seed": 5}, "data_dir": str(data_dir),
           "repo_root": str(repo_root),
           "parent": {"table": _niccolo.SPR_REL, "window_column": "cdrh3", "id": "trastuzumab"},
           "generators": [{"name": "scored_evolution", "n": 1500, "params": {}}],
           "filters": [{"name": "framework_integrity"}, {"name": "window_novelty", "params": {"min_dist": 3}},
                       {"name": "prefilter_gate", "params": {"min_p_binder": 0.8}}, {"name": "abguard_hard"},
                       {"name": "dedupe"}, {"name": "new_liabilities", "params": {"severities": ["high"], "max_new": 0}},
                       {"name": "track_check"}],
           "scorers": [{"name": "niccolo_prefilter"}],
           "aggregator": {"name": "borda", "params": {"columns": ["score__pred_log10_kd"]}}}
    out = tmp_path / "loop"
    summ = iterate.run_loop(cfg, rounds=2, n_select=15, out_dir=out, pool_frac=0.3, log=lambda m: None)
    assert [r["round"] for r in summ["rounds"]] == [0, 1] and "stopped" not in summ
    g1 = pl.read_parquet(out / "round_01" / "01_generated.parquet")
    assert {"scored_evolution:mutate", "scored_evolution:cross"} <= set(g1["origin"].to_list()), "round 1 must use the elites"
    s0, s1 = (pl.read_parquet(out / f"round_0{r}" / "selected.parquet") for r in (0, 1))
    assert s1.height >= 3 and s1["sequence_id"].n_unique() == s1.height
    assert s1["score__pred_log10_kd"].min() <= s0["score__pred_log10_kd"].min() + 1e-12      # elitism: never worse
    assert (s1["score__p_binder"] >= 0.8).all() and set(s1["info__track"]) == {"novel"}
    parent = data.load_parent(cfg["parent"], data_dir)
    assert (novelty.provided_index(Context(data_dir=data_dir, repo_root=repo_root, run_dir=tmp_path,
            cache_dir=out / "_cache", seed=0, parent=parent)).distance_class(
            novelty.encode([w for w in _niccolo.windows(s1["heavy_sequence"].to_list(), parent) if w])) >= 3).all()
