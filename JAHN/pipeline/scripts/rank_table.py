"""Score and rank any designs table with the frozen stage-1 scorers + NanoMelt tie-break.

    python scripts/rank_table.py DESIGNS.csv -o scored.csv

The table needs sequence_id, kind, heavy_sequence (light_sequence optional). What it does:
  - NICCOLO's pre-filter (classifier + ranker) on the CDRH3 window — the frozen stage-1 scores;
  - NanoMelt VH Tm and the delta against the parent (tie-break, per the PR #20 handoff);
  - borda over [pred_log10_kd, p_binder] with Tm at weight 0.25 — the saturated campaign's ranking.

Caveat carried from the training data: the pre-filter was fitted on the parent framework, so a
design that changes the framework (VL D28K, R87S) gets an approximation of its window's scores,
not a framework-aware prediction.
"""
from __future__ import annotations

import argparse

import polars as pl

from abpipe import _niccolo, registry, report
from abpipe.base import Context
from abpipe.data import load_parent
from abpipe.iterate import select_diverse


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("designs")
    ap.add_argument("-o", "--out", default="scored.csv")
    ap.add_argument("--min-pair-hamming", type=int, default=3)
    ap.add_argument("--no-nanomelt", action="store_true")
    a = ap.parse_args()

    import tempfile
    from pathlib import Path
    repo = Path(__file__).resolve().parents[3]
    data_dir = repo / "apheris/hackathon-data"
    parent = load_parent({"table": _niccolo.SPR_REL, "window_column": "cdrh3", "id": "trastuzumab"}, data_dir)

    df = pl.read_csv(a.designs)
    missing = [c for c in ("sequence_id", "kind", "heavy_sequence") if c not in df.columns]
    if missing:
        raise SystemExit(f"missing columns {missing}")
    df = df.with_columns(pl.col("light_sequence").fill_null("")) if "light_sequence" in df.columns else df
    skip = df.filter(pl.col("heavy_sequence").is_null() | (pl.col("heavy_sequence") == ""))
    if skip.height:
        print(f"note: {skip.height} row(s) without a heavy sequence are skipped: "
              f"{skip['sequence_id'].to_list()[:5]}")
    df = df.filter(pl.col("heavy_sequence").is_not_null() & (pl.col("heavy_sequence") != ""))

    with tempfile.TemporaryDirectory() as td:
        ctx = Context(data_dir=data_dir, repo_root=repo, run_dir=Path(td), cache_dir=Path(td), parent=parent)
        fit = _niccolo.get_prefilter(data_dir=data_dir, cache_dir=Path(td), repo_root=repo)
        wins = _niccolo.windows(df["heavy_sequence"].to_list(), parent)
        unscorable = [df["sequence_id"].to_list()[i] for i, w in enumerate(wins) if w is None]
        if unscorable:
            print(f"note: {len(unscorable)} design(s) not scorable by the pre-filter (framework length or "
                  f"residues): {unscorable[:5]}")
        s = _niccolo.score_windows(fit, wins)
        df = df.with_columns(pl.Series("score__p_binder", s["p_binder"]),
                             pl.Series("score__pred_log10_kd", s["pred_log10_kd"]))
        if not a.no_nanomelt:
            try:
                nm = registry.get("scorer", "nanomelt_tm")(ncpu=8)
                df = df.drop(["score__tm_pred", "score__dtm_vs_wt"]).join(
                    nm.score(df, ctx), on="sequence_id", how="left")
            except Exception as e:                                   # the binary or the ESM stage
                print(f"note: NanoMelt unavailable ({str(e)[-200:]}); ranking on the frozen columns")
                df = df.with_columns(pl.lit(None, dtype=pl.Float64).alias("score__tm_pred"),
                                     pl.lit(None, dtype=pl.Float64).alias("score__dtm_vs_wt"))
        else:
            df = df.with_columns(pl.lit(None, dtype=pl.Float64).alias("score__tm_pred"),
                                 pl.lit(None, dtype=pl.Float64).alias("score__dtm_vs_wt"))

    have_tm = ("score__dtm_vs_wt" in df.columns and df["score__dtm_vs_wt"].drop_nulls().len() > 0
               and not a.no_nanomelt)
    columns = (["score__pred_log10_kd", "score__p_binder", "score__dtm_vs_wt"] if have_tm
               else ["score__pred_log10_kd", "score__p_binder"])       # the frozen stage-1 aggregation
    agg = registry.get("aggregator", "borda")(columns=columns,
                                              weights=({"score__dtm_vs_wt": 0.25} if have_tm else None))
    directions = {"score__pred_log10_kd": "lower", "score__p_binder": "higher", "score__dtm_vs_wt": "higher"}
    ranked = agg.combine(df, directions, None)
    ranked = ranked.drop([c for c in ranked.columns if c.startswith("rank__")])
    ranked.write_csv(a.out)
    cols = [c for c in ("sequence_id", "rank", "score__pred_log10_kd", "score__p_binder",
                        "score__tm_pred", "score__dtm_vs_wt", "design_description") if c in ranked.columns]
    print(ranked.select(cols).head(25))
    print(f"scored {ranked.height} designs -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
