"""Union of the parallel orchestrated replicas: one ranked, deduplicated table from all runs.

    python scripts/union_replicas.py runs/orchestrated-optim runs/orchestrated-optim-s21 ... -o runs/union

Every replica ranks under the same frozen stage-1 scorers, so their selected sets are directly
comparable: pool them, re-rank with the same borda (pred_log10_kd, p_binder), and keep the top N
with the same pairwise-window diversity floor select_diverse enforces inside a single run. Every
design carries its run of origin; every run is a row in the generation ledger.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import polars as pl

from abpipe import _niccolo, registry
from abpipe.config import load as load_cfg
from abpipe.data import load_parent
from abpipe.iterate import select_diverse


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("runs", nargs="+", help="orchestrated run directories (each containing round_*/selected.parquet)")
    ap.add_argument("-c", "--config", default="configs/orchestrated_loop.yaml")
    ap.add_argument("-o", "--out", default="runs/union")
    ap.add_argument("-n", "--select", type=int, default=60)
    a = ap.parse_args()

    cfg = load_cfg(a.config)
    registry.load_plugins(cfg.plugins)
    ctx_parent = load_parent(cfg.parent, cfg.data_dir)
    agg = registry.get("aggregator", cfg.aggregator.name)(**cfg.aggregator.params)
    directions = {}
    for spec in cfg.scorers:
        directions.update(registry.get("scorer", spec.name)(**spec.params).outputs)

    frames = []
    for run in a.runs:
        for sel in sorted(Path(run).glob("round_*/selected.parquet")):
            df = pl.read_parquet(sel).with_columns(
                pl.lit(f"{Path(run).name}/{sel.parent.name}").alias("run_origin"))
            frames.append(df.drop([c for c in df.columns if c.startswith("rank")]))
    pool = pl.concat(frames, how="diagonal_relaxed")
    ranked = agg.combine(pool, directions, None)
    sel = select_diverse(ranked, ctx_parent, a.select)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    sel.write_parquet(out / "union_selected.parquet")
    scalar = [c for c, t in sel.schema.items() if not isinstance(t, (pl.List, pl.Struct, pl.Array))]
    sel.select(scalar).write_csv(out / "union_selected.csv")
    wins = [w for w in _niccolo.windows(sel["heavy_sequence"].to_list(), ctx_parent) if w]
    (out / "union_summary.json").write_text(json.dumps({
        "runs": [str(Path(r).name) for r in a.runs], "pooled": pool.height, "selected": sel.height,
        "distinct_windows": len(set(wins)),
        "origin_counts": sel["run_origin"].value_counts().sort("run_origin").to_dicts(),
        "score_medians": {c: sel[c].median() for c in sel.columns if c.startswith("score__")},
    }, indent=2, default=str))
    print(f"pooled {pool.height} selected rows -> {sel.height} diverse (top {a.select} requested)")
    print(sel["run_origin"].value_counts().sort("run_origin"))


if __name__ == "__main__":
    main()
