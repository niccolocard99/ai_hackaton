#!/usr/bin/env python3
"""Build the BASE VALIDATION table (AGENTS.md) from one or more RFD4 design runs.

Reads runs/<run>/design/generation/inference_metrics/per_example_metrics.csv for each run and writes

  designs.csv          sequence_id, kind, heavy_sequence, light_sequence, design_description
                       (the table check_denovo.py, check_optimisation.py and the generation ledger take;
                       a VHH goes whole in heavy_sequence, light_sequence stays empty)
  designs_metrics.csv  sequence_id plus the RFD4 metrics that matter for triage

`designed_sequences` in the metrics file is the binder chain with the retained framework: for these
designs, chain H. It ends ...QVTVS, exactly as generated (the organisers' framework drops the last S);
it is not edited, so the validated sequence is the designed one.

    ./make_designs_table.py RUN [RUN ...] [-o designs.csv] [--metrics designs_metrics.csv]

Pure standard library.
"""
import argparse
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
METRICS = "design/generation/inference_metrics/per_example_metrics.csv"

KEEP = {
    "hotspot_obedience": "hotspot_obedience.hotspot_obedience",
    "sequence_adherence": "sequence_adherence.sequence_adherence",
    "clash_frac_1p5A": "clash_metrics.frac_atoms_clashing_1p5A",
    "hbond_pairs": "interchain_interactions.hbond_atom_pair_count_coi",
    "saltbridge_pairs": "interchain_interactions.saltbridge_residue_pair_count_coi",
    "hydrophobic_pairs": "interchain_interactions.hydrophobic_atom_pair_count_coi",
    "pistack_pairs": "interchain_interactions.pistack_residue_pair_count_coi",
    "pre_oracle_pass": "pre_oracle_pass",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="run names under runs/ (e.g. vhh-fixedfw-n100-s2)")
    ap.add_argument("-o", "--out", default="designs.csv")
    ap.add_argument("--metrics", default="designs_metrics.csv")
    a = ap.parse_args()

    table, metrics, seen = [], [], {}
    for run in a.runs:
        path = HERE / "runs" / run / METRICS
        if not path.is_file():
            sys.exit(f"no metrics file for run {run}: {path}")
        for r in csv.DictReader(open(path, encoding="utf8")):
            chain, _, heavy = r["designed_sequences"].partition(":")
            if chain != "H" or not heavy:
                sys.exit(f"{r['example_id']}: expected 'H:<sequence>' in designed_sequences, got {r['designed_sequences'][:40]!r}")
            variant = "freefw" if "freefw" in r["task"] else "fixedfw"
            sid = f"{run}__{r['example_id']}"
            if heavy in seen:
                print(f"note: {sid} has the same sequence as {seen[heavy]}", file=sys.stderr)
            seen.setdefault(heavy, sid)
            table.append({
                "sequence_id": sid, "kind": "vhh", "heavy_sequence": heavy, "light_sequence": "",
                "design_description": (f"{variant} VHH on the trastuzumab epitope; hotspots D593 NZ, F573 CD1, I591 CG2, "
                                      f"P572 CD; RFD4 antibody-2026-10-30; run {run}; seed {int(float(r['seed']))}")})
            metrics.append({"sequence_id": sid, "variant": variant, "run": run,
                            **{k: r[v] for k, v in KEEP.items()}})

    with open(a.out, "w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, fieldnames=list(table[0]))
        w.writeheader()
        w.writerows(table)
    with open(a.metrics, "w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, fieldnames=list(metrics[0]))
        w.writeheader()
        w.writerows(metrics)
    print(f"wrote {len(table)} designs -> {a.out}, metrics -> {a.metrics}")


if __name__ == "__main__":
    main()
