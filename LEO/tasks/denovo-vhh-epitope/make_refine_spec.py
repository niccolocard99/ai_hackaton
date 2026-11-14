#!/usr/bin/env python3
"""Build an RFD4-Proteina refinement spec: close variants of chosen generated designs, by partial diffusion.

Why: in batch 1 only 8% of valid designs reached the team's epitope bar (epitope_cov >= 0.6 and
contacts_in_ep >= 0.7), and epitope coverage was the one axis with a usable signal (EXP-011). Rather than
re-rolling from scratch at that rate, re-run part of the sampling trajectory on the designs that already
reach the bar. Partial diffusion keeps the docked pose and varies loops and sequence, so the children stay on
the same epitope. The parent's annotations (fixed HER2, the four hotspots, the framework) travel with its
CIF, so nothing about the epitope is re-specified.

    ./make_refine_spec.py results/batch1-refine-parents.csv spec/design_refine_strong.json

Reads sequence_id (RUN__EXAMPLE) from the parents CSV, copies each runs/RUN/design/.../EXAMPLE_generated.cif
to spec/refine/EXAMPLE.cif, and writes one design per parent (input_as_annotated_cif true, so the parent's
conditioning is kept). Run it with run_design.sh in refine mode:

    PARTIAL=0.4 ./run_design.sh refine_strong 16 4

Per apheris/rfd4/04-refine.md: keep the parent's weights (antibody-2026-10-30) and recycle schedule, a
partial-trajectory fraction of 0.4, hotspot resampling off (never combine it with partial diffusion).
Smaller fractions stay closer to the parent. Evaluate children, do not assume they improve on the parent.
Pure standard library; run it where runs/ is (the VM).
"""
import csv
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    parents_csv, out = Path(sys.argv[1]), Path(sys.argv[2])
    rows = list(csv.DictReader(open(HERE / parents_csv, encoding="utf8")))
    refine_dir = (HERE / out).parent / "refine"
    refine_dir.mkdir(parents=True, exist_ok=True)
    spec, seen = {}, set()
    for r in rows:
        run, _, example = r["sequence_id"].partition("__")
        key = "refine_" + example.replace("her2d4_vhh_", "")
        if key in seen:
            sys.exit(f"two parents map to the same design id {key}; example ids repeat across runs")
        seen.add(key)
        src = HERE / "runs" / run / "design" / "generation" / "inference_samples" / example / f"{example}_generated.cif"
        if not src.is_file():
            sys.exit(f"parent CIF not found: {src}")
        shutil.copy2(src, refine_dir / f"{example}.cif")
        spec[key] = {"input": f"refine/{example}.cif", "input_as_annotated_cif": True}
        tail = r["heavy_sequence"][-8:]
        print(f"{key:<22} parent {r['sequence_id']}  CDRH3 {len(r['cdrh3'])} aa  cov {r['epitope_cov']}  ends ...{tail}")
    (HERE / out).write_text(json.dumps(spec, indent=2) + "\n", encoding="utf8")
    print(f"\n{len(spec)} parents -> {out}; copies in {refine_dir}")
    print("check the 'ends' column: refinement propagates an unwanted terminal tail if a parent has one")


if __name__ == "__main__":
    main()
