#!/usr/bin/env python3
"""Make a design-output tree that holds only chosen designs, so refold-prepare and refold-score can refold
exactly those.

refold-prepare takes only --top-n or --top-fraction, best-ranked by RFD4's own criteria, with no way to name
designs. It insists on a design-output directory: a .apheris-foundry marker, config.json, result.json,
generation/inference_metrics/per_example_metrics.csv and the generated CIFs. This builds that from the
shortlisted designs of one or more runs:

    ./select_for_refold.py results/batch1-shortlist.csv runs/shortlist1/in

For each sequence_id (RUN__EXAMPLE) in the shortlist it copies the design's CIF, takes its row from
runs/RUN/design/generation/inference_metrics/per_example_metrics.csv (paths rewritten to the new tree),
and writes a result.json listing exactly those designs, in shortlist order. Then pass OUT as --input and set
--top-n to the number of designs (refold.sh does this: DESIGN_DIR=OUT ./refold.sh NAME N MODEL).

Example ids (her2d4_vhh_fixedfw_0, ...) repeat across runs, so two shortlisted designs with the same id
cannot share one tree; the script stops if that happens. Pure standard library; run it where runs/ is (the VM).
"""
import csv
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DESIGN = Path("design")
PATH_COLS = ("designed_structure_path", "output_base_dir")


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    shortlist, out = Path(sys.argv[1]), Path(sys.argv[2])
    wanted = [r["sequence_id"] for r in csv.DictReader(open(shortlist, encoding="utf8"))]
    examples = [sid.partition("__")[2] for sid in wanted]
    dup = {e for e in examples if examples.count(e) > 1}
    if dup:
        sys.exit(f"example ids repeat across runs and would collide in one tree: {sorted(dup)}; split the shortlist")

    # container path of the new tree: the task folder is mounted at /work
    rel = out.relative_to(HERE) if out.is_absolute() else out
    base = f"/work/{rel.as_posix()}/generation/inference_samples"

    runs, header, rows, designs = {}, None, [], []
    for sid in wanted:
        run, _, example = sid.partition("__")
        if run not in runs:
            d = HERE / "runs" / run / DESIGN
            reader = csv.DictReader(open(d / "generation" / "inference_metrics" / "per_example_metrics.csv", encoding="utf8"))
            runs[run] = (d, {r["example_id"]: r for r in reader})
            if header is None:
                header = reader.fieldnames
            elif header != reader.fieldnames:
                sys.exit(f"run {run} has a different metrics header; refold its designs separately")
        d, table = runs[run]
        if example not in table:
            sys.exit(f"{sid}: example {example} not found in run {run}")
        src = d / "generation" / "inference_samples" / example / f"{example}_generated.cif"
        dst = HERE / out / "generation" / "inference_samples" / example / src.name
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        row = dict(table[example])
        for c in PATH_COLS:
            if c in row:
                row[c] = f"{base}/{example}/{src.name}" if c == "designed_structure_path" else base
        rows.append(row)
        designs.append(f"generation/inference_samples/{example}/{src.name}")

    first = next(iter(runs.values()))[0]
    mdir = HERE / out / "generation" / "inference_metrics"
    mdir.mkdir(parents=True, exist_ok=True)
    with open(mdir / "per_example_metrics.csv", "w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        w.writerows(rows)
    (HERE / out / ".apheris-foundry").write_text("", encoding="utf8")
    shutil.copy2(first / "config.json", HERE / out / "config.json")
    result = json.loads((first / "result.json").read_text(encoding="utf8"))
    result["result"]["designs"] = designs
    (HERE / out / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf8")
    print(f"{len(rows)} designs from {len(runs)} run(s) -> {HERE / out}")
    print(f"next: DESIGN_DIR={out} ./refold.sh NAME {len(rows)} of3ob")


if __name__ == "__main__":
    main()
