#!/usr/bin/env python3
"""BASE VALIDATION: run the team's real track checkers on designs, in parallel shards.

This does not re-implement abguard or the track rules. Every verdict comes from
NICCOLO/scripts/check_<track>.py, which wraps abguard (numbering, format, length,
unpaired cysteines, novelty) and adds that track's own rules. Only speed changes.

The optimisation checker's `--provided` test runs a pure-Python Levenshtein of the
whole VH and VL against all 768 provided antibodies, about 20 million DP cells per
design. rapidfuzz computes the identical unit-cost edit distance (insert, delete,
substitute each cost 1) in C. We swap it in before the checker module is imported,
the same way JAHN/pipeline does, and touch nothing else.

A design is eligible only if it shows PASS for the track abguard ASSIGNED it. FAIL
in the other track's script is expected, because a design is never reclassified.
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from her2lib import TZ_VL, subs, vh_with_loop  # noqa: E402

SCRIPTS = Path("/home/hackathon/wt/leo-design-loop/NICCOLO/scripts")
PROVIDED = Path("/home/hackathon/hackathon-team-3/apheris/hackathon-data/her2-trastuzumab/"
                "trastuzumab_her2_spr_all_768.csv")

WRAPPER = """
import importlib, sys
scripts_dir, module, *argv = sys.argv[1:]
sys.path.insert(0, scripts_dir)
import _track_common
from rapidfuzz.distance import Levenshtein
_track_common.levenshtein = Levenshtein.distance
mod = importlib.import_module(module)
sys.argv = [module + ".py", *argv]
raise SystemExit(mod.main())
"""


def design_id(loop: str) -> str:
    return f"LEO-{loop}"


def _table(loops: list[str], path: Path) -> None:
    rows = [{"sequence_id": design_id(l), "kind": "paired",
             "heavy_sequence": vh_with_loop(l), "light_sequence": TZ_VL,
             "design_description": f"CDRH3 {l} ({'/'.join(subs(l)) or 'wild-type'})"}
            for l in loops]
    pd.DataFrame(rows).to_csv(path, index=False)


def _run_shard(loops: list[str], track: str, tmp: Path, k: int) -> pd.DataFrame:
    table, report = tmp / f"designs_{k}.csv", tmp / f"report_{k}.csv"
    _table(loops, table)
    cmd = [sys.executable, "-c", WRAPPER, str(SCRIPTS), f"check_{track}", str(table), "-o", str(report)]
    if track == "optimisation":
        cmd += ["--provided", str(PROVIDED)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if not report.exists():
        raise RuntimeError(f"shard {k} produced no report (exit {p.returncode}):\n{p.stderr[-1500:]}")
    (tmp / f"stdout_{k}.txt").write_text(p.stdout, encoding="utf-8")
    return pd.read_csv(report)


def validate(loops: list[str], track: str = "optimisation", shards: int = 12,
             keep_dir: Path | None = None) -> pd.DataFrame:
    """Real checker verdicts for each loop. Columns follow the checker's report."""
    loops = list(dict.fromkeys(loops))
    if not loops:
        return pd.DataFrame()
    shards = max(1, min(shards, len(loops)))
    chunks = [loops[i::shards] for i in range(shards)]
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        with ThreadPoolExecutor(shards) as ex:
            parts = list(ex.map(lambda a: _run_shard(a[1], track, tmp, a[0]), enumerate(chunks)))
        if keep_dir:
            keep_dir.mkdir(parents=True, exist_ok=True)
            for f in tmp.glob("*_*.csv"):
                (keep_dir / f.name).write_bytes(f.read_bytes())
            for f in tmp.glob("stdout_*.txt"):
                (keep_dir / f.name).write_bytes(f.read_bytes())
    df = pd.concat(parts, ignore_index=True)
    df["loop"] = df.sequence_id.str.removeprefix("LEO-")
    return df


if __name__ == "__main__":
    import argparse
    import time

    ap = argparse.ArgumentParser()
    ap.add_argument("csv", type=Path, help="table with a cdrh3 column (10-mers)")
    ap.add_argument("--track", default="optimisation", choices=["optimisation", "denovo"])
    ap.add_argument("--shards", type=int, default=12)
    ap.add_argument("-o", type=Path, required=True)
    a = ap.parse_args()
    loops = pd.read_csv(a.csv).cdrh3.tolist()
    t0 = time.time()
    r = validate(loops, a.track, a.shards)
    r.to_csv(a.o, index=False)
    print(f"{len(r)} designs, {int(r.passes.sum())} PASS ({a.track}) in {time.time() - t0:.0f}s -> {a.o}")
    print(r.track_assigned.value_counts().to_string())


def validate_rows(rows: list[dict], track: str = "optimisation", shards: int = 4) -> pd.DataFrame:
    """Real checker verdicts for arbitrary designs: dicts with sequence_id, kind, heavy_sequence,
    light_sequence, design_description. Used for designs that differ from trastuzumab outside CDRH3."""
    shards = max(1, min(shards, len(rows)))
    chunks = [rows[i::shards] for i in range(shards)]
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)

        def run(a):
            k, part = a
            table, report = tmp / f"t_{k}.csv", tmp / f"r_{k}.csv"
            pd.DataFrame(part).to_csv(table, index=False)
            cmd = [sys.executable, "-c", WRAPPER, str(SCRIPTS), f"check_{track}", str(table), "-o", str(report)]
            if track == "optimisation":
                cmd += ["--provided", str(PROVIDED)]
            p = subprocess.run(cmd, capture_output=True, text=True)
            if not report.exists():
                raise RuntimeError(f"shard {k} produced no report:\n{p.stderr[-1200:]}")
            return pd.read_csv(report)

        with ThreadPoolExecutor(shards) as ex:
            return pd.concat(list(ex.map(run, enumerate(chunks))), ignore_index=True)
