"""Contract test for the frozen stage 1: the columns abpipe depends on must keep existing.

Why this file is separate from test_niccolo_loop.py's `needs_data` group: that group is gated on `abguard`,
which is not installed on most machines, so it skips silently. A change to NICCOLO/eda/prefilter.py once
removed the `p_facs_high` column and broke the scorer, the gate filter and `emit_facs` with a KeyError, and
the suite still reported green. This test needs only the SPR table -- no abguard, no FACS splits -- so it
runs almost everywhere and fails loudly instead.
"""
import os
import pathlib
import shutil

import polars as pl
import pytest

from abpipe import _niccolo
from abpipe.base import Context, Parent

SPR = "trastuzumab_her2_spr_all_768.csv"
REQUIRED_SCORE_COLUMNS = {"cdrh3", "n_mut", "p_binder", "pred_log10_kd"}


def _find_spr():
    """The SPR table, wherever this machine keeps it. It is git-ignored, so the path varies."""
    repo = _niccolo.repo_guess()
    for p in (os.environ.get("ABPIPE_DATA_DIR") and pathlib.Path(os.environ["ABPIPE_DATA_DIR"]) / "her2-trastuzumab" / SPR,
              pathlib.Path(os.path.expanduser("~/hackathon-team-3/apheris/hackathon-data/her2-trastuzumab")) / SPR,
              repo / "NICCOLO/eda/data" / SPR,
              repo / "NICCOLO/eda/data/her2-trastuzumab" / SPR):
        if p and pathlib.Path(p).exists():
            return pathlib.Path(p)
    return None


@pytest.fixture(scope="module")
def staged(tmp_path_factory):
    src = _find_spr()
    if src is None:
        pytest.skip(f"{SPR} not found; set ABPIPE_DATA_DIR (see NICCOLO/eda/AGENT_GUIDE.md section 2)")
    d = tmp_path_factory.mktemp("data") / "d"
    (d / "her2-trastuzumab").mkdir(parents=True)
    shutil.copy(src, d / "her2-trastuzumab" / SPR)
    return d


@pytest.fixture(scope="module")
def fitted(staged, tmp_path_factory):
    return _niccolo.get_prefilter(data_dir=staged, cache_dir=tmp_path_factory.mktemp("cache"),
                                  repo_root=_niccolo.repo_guess())


def test_runs_without_the_facs_splits(fitted):
    """The default model does not use the FACS feature, so the 11 MB clonotype splits must not be required."""
    assert _niccolo.has_facs(fitted) is False
    assert fitted.clf_w.shape == (_niccolo.L, 20) and fitted.rank_w.shape == (_niccolo.L, 20)


def test_score_emits_the_columns_the_pipeline_reads(fitted):
    out = fitted.pf.score(["WGGDGFYAMD", "WAADGFYSIE"])
    missing = REQUIRED_SCORE_COLUMNS - set(out.columns)
    assert not missing, (f"prefilter.PreFilter.score stopped emitting {sorted(missing)}. abpipe reads these; "
                         "update abpipe/_niccolo.py and NICCOLO/eda/AGENT_GUIDE.md together, or restore them.")
    assert "p_facs_high" not in out.columns, "use_facs=False must not emit p_facs_high"


def test_score_windows_and_the_scorer_agree(fitted, staged, tmp_path):
    import pandas as pd
    row = pd.read_csv(staged / "her2-trastuzumab" / SPR).query("cdrh3 == 'WGGDGFYAMD'").iloc[0]
    heavy, light = row.vh_sequence, row.vl_sequence
    start = heavy.index("WGGDGFYAMD")
    parent = Parent(id="trastuzumab", heavy=heavy, light=light, cdrh3="SRWGGDGFYAMDY",
                    window="WGGDGFYAMD", window_start=start)
    ctx = Context(data_dir=staged, repo_root=_niccolo.repo_guess(), run_dir=tmp_path, cache_dir=tmp_path,
                  seed=0, parent=parent)
    wins = ["WAADGFYSIE", "WGGSGPYLFD", "AAAAAAAAAA"]
    df = pl.DataFrame({"sequence_id": ["a", "b", "c"],
                       "heavy_sequence": [heavy[:start] + w + heavy[start + 10:] for w in wins],
                       "light_sequence": [light] * 3})
    res = _niccolo.score_windows(fitted, _niccolo.windows(df["heavy_sequence"].to_list(), parent))
    assert {"p_binder", "pred_log10_kd"} <= set(res)

    from abpipe.scorers.niccolo_prefilter import NiccoloPrefilter
    scored = NiccoloPrefilter().score(df, ctx)
    assert set(scored.columns) == {"sequence_id", "score__p_binder", "score__pred_log10_kd"}
    direct = fitted.pf.score(wins)
    assert scored["score__p_binder"].to_list() == pytest.approx(direct["p_binder"].to_list())
    assert scored["score__pred_log10_kd"].to_list() == pytest.approx(direct["pred_log10_kd"].to_list())

    from abpipe.filters.prefilter_gate import PrefilterGate
    reasons = PrefilterGate(min_p_binder=0.5).check(df, ctx)["reason"].to_list()
    assert reasons[2] is not None, "poly-alanine must be rejected by the gate"

    with pytest.raises(_niccolo.PrefilterUnavailable, match="emit_facs"):
        NiccoloPrefilter(emit_facs=True).score(df, ctx)


def test_the_scope_guard_still_refuses_other_parents(fitted, staged, tmp_path):
    other = Parent(id="x", heavy="Q" * 120, light="D" * 107, cdrh3="SRAAAAAAAAAAY",
                   window="AAAAAAAAAA", window_start=98)
    with pytest.raises(_niccolo.PrefilterUnavailable, match="specific to trastuzumab"):
        _niccolo.check_parent(other, fitted)
