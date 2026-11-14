#!/usr/bin/env python3
"""Measure what abguard actually enforces, instead of assuming it.

Three things decide how many substitutions a legal optimisation design needs,
and all three are abguard's to define, not ours:

  1. which residues abguard calls the CDRH3 (the 10-mer we design, or a wider
     window including the SR prefix and the Y suffix),
  2. what identity it computes against trastuzumab's CDRH3,
  3. how it counts the distance to the nearest reference CDRH3.

The ledger records a disagreement on exactly this: AGENTS.md reasons from an
11-mer and the pre-filter guide from a 13-mer, giving different substitution
floors. So this probe walks a ladder of 0 to 6 substitutions in the designed
window and prints what abguard says about each.

Run from the task folder. Writes nothing; prints a table.
"""
from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from her2lib import TZ_LOOP, design_rows  # noqa: E402

HERE = Path(__file__).resolve().parent
CHECKERS = Path("/home/hackathon/wt/leo-design-loop/NICCOLO/scripts")
SPR = HERE / "data" / "trastuzumab_her2_spr_all_768.csv"

# A ladder of substitutions placed on the never-contacting positions first
# (106, 107, 108 -> indices 7, 8, 9), then outward to contacting ones.
# Each rung adds one substitution to the rung above it.
LADDER = [
    ("0 subs (wild type)",        TZ_LOOP),
    ("1 sub  (A106L)",            "WGGDGFYLMD"),
    ("2 subs (+M107F)",           "WGGDGFYLFD"),
    ("3 subs (+D108V)",           "WGGDGFYLFV"),
    ("4 subs (+D102S)",           "WGGSGFYLFV"),
    ("5 subs (+F104P)",           "WGGSGPYLFV"),
    ("6 subs (+G103S)",           "WGGSSPYLFV"),
]


def main() -> int:
    rows = list(design_rows([loop for _, loop in LADDER], prefix="PROBE"))
    labels = {r["sequence_id"]: lab for r, (lab, _) in zip(rows, LADDER)}
    loops = {r["sequence_id"]: loop for r, (_, loop) in zip(rows, LADDER)}

    with tempfile.TemporaryDirectory() as td:
        table = Path(td) / "probe.csv"
        with table.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

        report = Path(td) / "report.csv"
        cmd = [sys.executable, str(CHECKERS / "check_optimisation.py"), str(table),
               "-o", str(report), "--provided", str(SPR)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        print("--- checker stdout ---")
        print(proc.stdout.strip() or "(none)")
        if proc.returncode not in (0, 1):
            print("--- stderr ---")
            print(proc.stderr.strip()[-3000:])
            return proc.returncode
        if proc.stderr.strip():
            print("--- stderr (non-fatal) ---")
            print(proc.stderr.strip()[-1500:])

        out = list(csv.DictReader(report.open(encoding="utf-8")))

    if not out:
        print("no report rows")
        return 1

    print("\n--- columns abguard returned ---")
    print(", ".join(out[0].keys()))

    print("\n--- the ladder ---")
    hdr = f"{'design':<22} {'loop':<11} {'abguard cdrh3':<16} {'ident%':>7} {'minDist':>7} {'nearest':<16} {'pass':>5}"
    print(hdr)
    print("-" * len(hdr))
    for r in out:
        sid = r["sequence_id"]
        ident = r.get("cdrh3_identity_to_trastuzumab", "")
        try:
            ident = f"{float(ident) * 100:.1f}"
        except ValueError:
            ident = "?"
        print(f"{labels[sid]:<22} {loops[sid]:<11} {r.get('cdrh3', ''):<16} {ident:>7} "
              f"{r.get('cdrh3_min_distance', ''):>7} {r.get('cdrh3_nearest', ''):<16} "
              f"{r.get('passes', ''):>5}")
        fails = r.get("failures", "")
        if fails and fails not in ("[]", ""):
            print(f"{'':<22} -> {fails}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
