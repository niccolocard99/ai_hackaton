"""Shared plumbing for check_denovo.py and check_optimisation.py.

Runs abguard's own annotation (immunum numbering + the shared hard checks) and
hands back one row per design with the novelty columns, so each track script only
adds its track-specific rules. Nothing here re-implements abguard's checks.
"""
from __future__ import annotations

import argparse
import glob
import os
import shutil
import sys
from pathlib import Path


def _ensure_abguard() -> None:
    """Re-exec under the Python that has abguard installed, if this one lacks it."""
    try:
        import abguard  # noqa: F401
        return
    except ImportError:
        pass
    candidates = [os.environ.get("ABGUARD_PY", "")]
    exe = shutil.which("abguard")
    if exe:
        candidates.append(str(Path(exe).resolve().parent / "python"))
    # newest uv-cache env that has the novelty gate (the older build lacks --reference-db)
    envs = glob.glob(os.path.expanduser("~/.cache/uv/archive-v0/*/lib/python3*/site-packages/abguard/novelty.py"))
    for novelty in sorted(envs, key=os.path.getmtime, reverse=True):
        candidates.append(str(Path(novelty).parents[4] / "bin" / "python"))
    for py in candidates:
        if py and os.path.exists(py) and os.environ.get("_TRACK_REEXEC") != "1":
            os.environ["_TRACK_REEXEC"] = "1"
            os.execv(py, [py, *sys.argv])
    sys.exit("abguard is not importable. Install it, or set ABGUARD_PY=/path/to/python that has it.")


_ensure_abguard()

import polars as pl  # noqa: E402
from abguard.annotate import annotate  # noqa: E402
from abguard.io import read_table  # noqa: E402
from abguard.novelty import trastuzumab_cdrh3  # noqa: E402
from abguard.reference import bundled_reference, cdrh3_reference  # noqa: E402
from abguard.config import PipelineConfig  # noqa: E402

IDENTITY_THRESHOLD = 0.70


def levenshtein(a: str, b: str) -> int:
    if a == b:
        return 0
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def cdrh3_identity(a: str, b: str) -> float:
    return 1 - levenshtein(a, b) / max(len(a), len(b))


def build_parser(description: str) -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=description)
    p.add_argument("submission", type=Path,
                   help="table with sequence_id, kind, heavy_sequence, light_sequence, design_description")
    p.add_argument("--reference-db", type=Path, default=None,
                   help="CDRH3 reference (Parquet bundle or DuckDB). Default: the one bundled with abguard")
    p.add_argument("-o", "--output", type=Path, default=None, help="write the per-design report (csv/parquet/json)")
    p.add_argument("--strict-exit", action="store_true", help="exit 1 if any design fails")
    return p


def annotated(path: Path, reference_db: Path | None) -> pl.DataFrame:
    config = PipelineConfig()
    ref = cdrh3_reference(reference_db or bundled_reference(), config=config)
    # gate_novelty=False: abguard's own novelty gate is track-agnostic; the track
    # scripts apply the track-specific thresholds themselves.
    return annotate(read_table(path), config, include_liabilities=False, reference=ref, gate_novelty=False)


def finish(df: pl.DataFrame, track: str, args: argparse.Namespace) -> int:
    cols = ["sequence_id", "kind", "track_assigned", "identity_to_trastuzumab", "cdrh3",
            "cdrh3_min_distance", "cdrh3_identity_to_trastuzumab", "passes", "failures"]
    out = df.select([c for c in cols if c in df.columns])
    if args.output:
        suffix = args.output.suffix
        if suffix == ".parquet":
            out.write_parquet(args.output)
        elif suffix == ".json":
            out.write_json(args.output)
        else:
            out.with_columns(pl.col("failures").list.join("; ")).write_csv(args.output)
    n_pass = int(df["passes"].sum())
    for row in df.iter_rows(named=True):
        status = "PASS" if row["passes"] else "FAIL"
        why = "" if row["passes"] else "  <- " + "; ".join(row["failures"])
        print(f"{status}  {row['sequence_id']:<24} assigned={row['track_assigned'] or '?':<8}{why}")
    print(f"\n{track} track: {n_pass}/{df.height} design(s) eligible")
    return 1 if (args.strict_exit and n_pass < df.height) else 0
