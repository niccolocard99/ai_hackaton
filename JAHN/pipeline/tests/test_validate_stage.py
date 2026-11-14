"""The `validate` stage: expensive checks after ranking, on the top of the list only.

Cheap-first ordering breaks down at the end of the funnel. abguard's checks cost
~0.7 ms/design for numbering (and the novelty gate 26 ms/design unpatched), so running
them on every survivor of a 1M-candidate screen dominates the run -- while only the top
of the ranked list can ever be submitted. `validators` run after `rank`, capped by
`validate_top`, so the gate sees the designs that matter and nothing else.

These tests need neither abguard nor the organiser data.
"""
from __future__ import annotations

import json
import sys
import textwrap

import polars as pl
import pytest

from abpipe import config, runner

PLUGIN = textwrap.dedent("""
    import polars as pl
    from abpipe.base import Filter, Generator, Scorer
    from abpipe.registry import register

    @register("generator", "v_gen")
    class VGen(Generator):
        \"\"\"toy sequences of increasing length\"\"\"
        def generate(self, n, seed, ctx):
            return pl.DataFrame([dict(sequence_id=f"v{i}", kind="paired",
                                      heavy_sequence="E" * (10 + i), light_sequence="D" * 5,
                                      design_description="toy") for i in range(n)])

    @register("scorer", "v_len")
    class VLen(Scorer):
        outputs = {"score__len": "higher"}
        def score(self, df, ctx):
            return df.select("sequence_id",
                             score__len=pl.col("heavy_sequence").str.len_chars().cast(pl.Float64))

    @register("filter", "v_min_len")
    class VMinLen(Filter):
        \"\"\"stands in for abguard_hard: the expensive check we want to run late\"\"\"
        def __init__(self, min_len: int = 17):
            self.min_len = min_len
        def check(self, df, ctx):
            return df.select("sequence_id",
                             reason=pl.when(pl.col("heavy_sequence").str.len_chars() < self.min_len)
                             .then(pl.lit("too short")).otherwise(None))

    WINDOWS = ["WWWWWWWWWW", "WWWWWWWWWP", "WWWWWWWWPP", "WWWWWWWPPP",
               "WWWWWWPPPP", "WWWWWGGGGG", "WWWWGGGGGG", "WWWGGGGGGG"]

    @register("generator", "v_window_gen")
    class VWindowGen(Generator):
        \"\"\"the real geometry: one fixed scaffold, a varying window at a fixed offset\"\"\"
        def generate(self, n, seed, ctx):
            return pl.DataFrame([dict(sequence_id=f"v{i}", kind="paired",
                                      heavy_sequence="Q" * 20 + WINDOWS[i] + "Q" * 10,
                                      light_sequence="D" * 10, design_description="toy")
                                 for i in range(min(n, len(WINDOWS)))])

    @register("scorer", "v_bcount")
    class VBCount(Scorer):
        outputs = {"score__b": "higher"}
        def score(self, df, ctx):
            return df.select("sequence_id",
                             score__b=pl.col("heavy_sequence").str.slice(20, 10)
                             .str.count_matches("W").cast(pl.Float64))

    @register("filter", "v_no_x")
    class VNoX(Filter):
        \"\"\"stands in for abguard: the authoritative gate, run after ranking\"\"\"
        def check(self, df, ctx):
            return df.select("sequence_id",
                             reason=pl.when(pl.col("heavy_sequence").str.slice(20, 10).str.contains("P"))
                             .then(pl.lit("contains P")).otherwise(None))

    @register("filter", "v_counting")
    class VCounting(Filter):
        \"\"\"records how many rows it was asked to check\"\"\"
        seen = []
        def check(self, df, ctx):
            VCounting.seen.append(df.height)
            return df.select("sequence_id", reason=pl.lit(None, dtype=pl.String))
""")


@pytest.fixture(scope="module")
def plugin(tmp_path_factory):
    directory = tmp_path_factory.mktemp("vplugin")
    (directory / "abpipe_validate_plugin.py").write_text(PLUGIN)
    sys.path.insert(0, str(directory))
    yield "abpipe_validate_plugin"
    sys.path.remove(str(directory))


def _cfg(tmp_path, plugin, *, validators=None, validate_top=None, n=8):
    raw = {
        "run": {"name": "v", "out_dir": str(tmp_path / "runs"), "seed": 1},
        "data_dir": str(tmp_path), "repo_root": str(tmp_path), "plugins": [plugin],
        "generators": [{"name": "v_gen", "n": n}],
        "filters": [],
        "scorers": [{"name": "v_len"}],
        "aggregator": {"name": "borda"},
    }
    if validators is not None:
        raw["validators"] = validators
    if validate_top is not None:
        raw["validate_top"] = validate_top
    return config.load(raw)


def _validated(cfg) -> pl.DataFrame:
    return pl.read_parquet(cfg.out_dir / cfg.name / "05_validated.parquet")


def test_validators_see_only_the_top_ranked_rows(tmp_path, plugin):
    import abpipe_validate_plugin as plug

    plug.VCounting.seen.clear()
    cfg = _cfg(tmp_path, plugin, validators=[{"name": "v_counting"}], validate_top=3)
    runner.run(cfg)

    assert plug.VCounting.seen == [3], "the validator must be handed the top 3 only, not all 8"


def test_validate_keeps_the_survivors_in_rank_order(tmp_path, plugin):
    # score__len is higher-is-better, so v7 (len 17) ranks first, then v6, v5...
    cfg = _cfg(tmp_path, plugin, validators=[{"name": "v_min_len", "params": {"min_len": 16}}],
               validate_top=4)
    runner.run(cfg)

    out = _validated(cfg)
    assert out.get_column("sequence_id").to_list() == ["v7", "v6"]
    assert out.get_column("rank").to_list() == sorted(out.get_column("rank").to_list())


def test_validate_is_skipped_when_no_validators_are_configured(tmp_path, plugin):
    cfg = _cfg(tmp_path, plugin)
    manifest = runner.run(cfg)

    assert "validate" not in manifest["stages"]
    assert not (cfg.out_dir / cfg.name / "05_validated.parquet").exists()


def test_validate_without_a_cap_checks_the_whole_ranked_list(tmp_path, plugin):
    import abpipe_validate_plugin as plug

    plug.VCounting.seen.clear()
    cfg = _cfg(tmp_path, plugin, validators=[{"name": "v_counting"}])
    runner.run(cfg)

    assert plug.VCounting.seen == [8]


def test_validate_records_its_rejections_in_the_manifest(tmp_path, plugin):
    cfg = _cfg(tmp_path, plugin, validators=[{"name": "v_min_len", "params": {"min_len": 16}}],
               validate_top=4)
    manifest = runner.run(cfg)

    stage = manifest["stages"]["validate"]
    assert stage["rows_in"] == 4
    assert stage["rows_out"] == 2
    assert stage["rejections"]

    rejections = pl.read_parquet(cfg.out_dir / cfg.name / "05_validate_rejections.parquet")
    assert set(rejections.get_column("sequence_id").to_list()) == {"v5", "v4"}


def test_validate_stage_can_be_rerun_alone(tmp_path, plugin):
    cfg = _cfg(tmp_path, plugin, validators=[{"name": "v_min_len", "params": {"min_len": 16}}],
               validate_top=4)
    runner.run(cfg)
    runner.run(cfg, stages=("validate",))

    assert _validated(cfg).get_column("sequence_id").to_list() == ["v7", "v6"]


def test_validate_top_must_be_positive(tmp_path, plugin):
    with pytest.raises(config.ConfigError):
        _cfg(tmp_path, plugin, validators=[{"name": "v_counting"}], validate_top=0)


def test_the_loop_never_selects_a_design_that_failed_validation(tmp_path, plugin):
    """The reason moving abguard after `rank` is not just a reordering.

    `run_loop` picks its elites from `03_scored.parquet`, which is upstream of the
    validate stage. If selection ignored the validators, designs that abguard rejects
    would seed the next round and reach `selected.csv` -- an eligibility leak with no
    symptom until submission. Selection must be restricted to validated survivors.
    """
    import abpipe.data
    from abpipe.base import Parent
    from abpipe.iterate import run_loop

    parent = Parent(id="toy", heavy="Q" * 20 + "WWWWWWWWWW" + "Q" * 10, light="D" * 10,
                    cdrh3="WWWWWWWWWW", window="WWWWWWWWWW", window_start=20)
    monkey = pytest.MonkeyPatch()
    monkey.setattr(abpipe.data, "load_parent", lambda spec, data_dir: parent)
    try:
        raw = {
            "run": {"name": "vloop", "out_dir": str(tmp_path / "runs"), "seed": 1},
            "data_dir": str(tmp_path), "repo_root": str(tmp_path), "plugins": [plugin],
            "parent": {"table": "unused.csv"},
            "generators": [{"name": "v_window_gen", "n": 8}],
            "filters": [],
            "scorers": [{"name": "v_bcount"}],
            "aggregator": {"name": "borda"},
            # The expensive gate, after ranking.
            "validators": [{"name": "v_no_x"}],
            "validate_top": 8,
        }
        summary = run_loop(raw, rounds=1, n_select=8, out_dir=str(tmp_path / "loop"),
                           min_pair_hamming=1, max_sub_share=1.0, pool_frac=1.0,
                           log=lambda m: None)
        assert "stopped" not in summary
    finally:
        monkey.undo()

    selected = pl.read_parquet(tmp_path / "loop" / "round_00" / "selected.parquet")
    # Only the four P-free windows clear the validator; the rest must never be selected.
    assert set(selected.get_column("sequence_id").to_list()) == {"v0", "v5", "v6", "v7"}
