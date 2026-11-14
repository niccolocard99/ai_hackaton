"""ANDERS/vhh_variants.py as an abpipe generator: de novo VHH variants from a scored/ranked list.

The tool (ANDERS/vhh_variants.py, merged in PR #11) is run unchanged as a subprocess; this module
only maps abpipe's contract onto its CLI and keeps its metadata as info__* columns:

  n            -> --per-strategy (n / strategies, floor 2) and --max-output n
  parents      -> the previous round's selected.csv (overrides the cold-start `ranked` list), with a
                  sha256 pin like scored_evolution, so a loop round never generates from a stale table
  score_col    -> --score-col (the loop's selected.csv carries `rank_score`, which the tool detects)

Everything the tool pre-screens (identity <= identity_max, CDRH3 >= cdrh3_dist_min edits from every
reference loop and from trastuzumab) is a MARGIN inside the official de novo gate; track_check remains
the authoritative verdict, as with every other module.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import polars as pl

from abpipe.base import Context, Generator
from abpipe.generators.scored_evolution import file_sha256
from abpipe.registry import register

_TOOL = "ANDERS/vhh_variants.py"


@register("generator", "vhh_variants")
class VhhVariants(Generator):
    """De-novo-track VHH variants by recombine / extend / mutate / graft from a ranked parent list."""

    version = "1"

    def __init__(self, ranked: str | None = None, parents: str | None = None, parents_sha256: str | None = None,
                 score_col: str | None = None, lower_is_better: bool | None = None, parents_top: int = 6,
                 strategies: str = "recombine,extend,mutate,graft", identity_max: float = 0.62,
                 cdrh3_dist_min: int = 4, run_official: bool = False, python_bin: str | None = None,
                 tool: str | None = None, timeout_s: int = 3600):
        if ranked is None and parents is None:
            raise ValueError("vhh_variants needs `ranked` (a scored/ranked CSV) for the cold-start round, "
                             "or `parents` (the previous round's selected.csv) once the loop runs")
        self.ranked, self.parents, self.parents_sha256 = ranked, parents, parents_sha256
        self.score_col, self.lower_is_better = score_col, lower_is_better
        self.parents_top, self.strategies = parents_top, strategies
        self.identity_max, self.cdrh3_dist_min = identity_max, cdrh3_dist_min
        self.run_official, self.python_bin, self.tool, self.timeout_s = run_official, python_bin, tool, timeout_s

    def generate(self, n: int, seed: int, ctx: Context) -> pl.DataFrame:
        tool = Path(self.tool) if self.tool else ctx.repo_root / _TOOL
        if not tool.exists():
            raise FileNotFoundError(f"{tool} not found (ANDERS generator, PR #11)")
        src = self.parents or self.ranked
        if self.parents and self.parents_sha256 and file_sha256(self.parents) != self.parents_sha256:
            raise ValueError(f"{self.parents} changed since the config was written (sha256 mismatch)")
        n_strat = max(1, len([s for s in self.strategies.split(",") if s.strip()]))
        out = ctx.run_dir / "_vhh_variants"
        out.mkdir(parents=True, exist_ok=True)
        cmd = [self.python_bin or sys.executable, str(tool), "generate", str(Path(src).expanduser()),
               "-o", str(out), "--seed", str(seed), "--parents-top", str(self.parents_top),
               "--strategies", self.strategies, "--per-strategy", str(max(2, n // n_strat)),
               "--max-output", str(n), "--identity-max", str(self.identity_max),
               "--cdrh3-dist-min", str(self.cdrh3_dist_min)]
        if self.score_col:
            cmd += ["--score-col", self.score_col]
        if self.lower_is_better:
            cmd += ["--lower-is-better"]
        if not self.run_official:
            cmd += ["--no-official"]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=self.timeout_s, cwd=str(ctx.repo_root))
        designs, meta = out / "variants.csv", out / "variants_meta.csv"
        if proc.returncode != 0 or not designs.exists():
            raise RuntimeError(f"vhh_variants.py failed (exit {proc.returncode}):\n"
                               f"{(proc.stdout or '')[-800:]}\n{(proc.stderr or '')[-800:]}")
        df = pl.read_csv(designs, schema_overrides={"sequence_id": pl.String})
        m = pl.read_csv(meta, schema_overrides={"sequence_id": pl.String})
        df = df.with_columns(pl.col("light_sequence").fill_null(""))
        df = df.join(m.select("sequence_id", "strategy", "parent_score", "identity_to_trastuzumab", "cdrh3",
                              "cdrh3_min_distance", "cdrh3_dist_trastuzumab"), on="sequence_id", how="left")
        rows = []
        for i, r in enumerate(df.iter_rows(named=True)):
            rows.append({
                "sequence_id": f"{self.name}_s{seed}_{i:04d}", "kind": r["kind"],
                "heavy_sequence": r["heavy_sequence"], "light_sequence": r["light_sequence"],
                "design_description": r["design_description"], "parent_id": self.parents or self.ranked,
                "origin": f"{self.name}:{r['strategy']}", "seed": seed,
                "meta": json.dumps({k: r[k] for k in ("strategy", "parent_score", "identity_to_trastuzumab",
                                                      "cdrh3", "cdrh3_min_distance", "cdrh3_dist_trastuzumab")}),
                "info__strategy": r["strategy"], "info__parent_score": r["parent_score"],
                "info__identity": r["identity_to_trastuzumab"], "info__cdrh3": r["cdrh3"],
                "info__cdrh3_min_dist": r["cdrh3_min_distance"]})
        return pl.DataFrame(rows, schema={"sequence_id": pl.String, "kind": pl.String, "heavy_sequence": pl.String,
                                           "light_sequence": pl.String, "design_description": pl.String,
                                           "parent_id": pl.String, "origin": pl.String, "seed": pl.Int64,
                                           "meta": pl.String, "info__strategy": pl.String,
                                           "info__parent_score": pl.Float64, "info__identity": pl.Float64,
                                           "info__cdrh3": pl.String, "info__cdrh3_min_dist": pl.Int64})
