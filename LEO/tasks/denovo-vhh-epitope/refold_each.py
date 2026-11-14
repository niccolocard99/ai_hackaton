#!/usr/bin/env python3
"""Refold designs one at a time, then collect the scores. Run ON THE VM.

refold-prepare refuses a batch whose designs do not share the same conditioned residues ("refold each
target separately"). Designs with different CDR lengths shift the framework residue indices, so they do not
share them. For each design in the list this builds a one-design tree (select_for_refold.py), then runs
refold.sh on it (prepare, fold, score), and finally merges every refold_metrics.csv into one summary.

<<<<<<< HEAD
    FORCE=1 ./refold_each.py results/batch1-shortlist.csv shortlist1 of3ob [N]

LIST is a CSV with a sequence_id column (RUN__EXAMPLE), e.g. the shortlist from score_designs.py. NAME is the
output folder under runs/. MODEL is of3ob or protenix-v2 (no Boltz: team decision). N limits the run to the first N designs.
=======
    FORCE=1 ./refold_each.py results/batch1-shortlist.csv shortlist1 boltz2 [N]

LIST is a CSV with a sequence_id column (RUN__EXAMPLE), e.g. the shortlist from score_designs.py. NAME is the
output folder under runs/. MODEL is boltz2, of3ob or protenix-v2. N limits the run to the first N designs.
>>>>>>> origin/main
Output: runs/NAME/refold_summary.csv. Each design takes about one docker start per step, so budget minutes
per design, not seconds. FORCE=1 is passed through to the GPU guard in refold.sh.

Use the scores as a sanity check on pose and hotspot contacts, not to rank (NICCOLO/eda/AGENT_GUIDE.md
section 6). Compare scores only within one folding model.
"""
import csv
import os
import subprocess
import sys
<<<<<<< HEAD
import time
=======
>>>>>>> origin/main
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    list_csv, name, model = sys.argv[1:4]
    n = int(sys.argv[4]) if len(sys.argv) > 4 else None
    reader = csv.DictReader(open(HERE / list_csv, encoding="utf8"))
    rows = list(reader)[:n]
    summary, failed = [], []
    for i, r in enumerate(rows, 1):
        sid = r["sequence_id"]
        d = Path("runs") / name / sid
        (HERE / d).mkdir(parents=True, exist_ok=True)
        one = HERE / d / "one.csv"
        with open(one, "w", newline="", encoding="utf8") as f:
            w = csv.DictWriter(f, fieldnames=reader.fieldnames)
            w.writeheader()
            w.writerow(r)
        print(f"\n##### {i}/{len(rows)}  {sid}", flush=True)
        subprocess.run([sys.executable, str(HERE / "select_for_refold.py"), str(one), str(d / "in")], check=True, cwd=HERE)
<<<<<<< HEAD
        for _ in range(180):  # refold.sh exits 3 only from its GPU guard (another oslo-hackathon container is running)
            res = subprocess.run([str(HERE / "refold.sh"), f"{name}/{sid}", "1", model],
                                 env={**os.environ, "DESIGN_DIR": str(d / "in")}, cwd=HERE)
            if res.returncode != 3:
                break
            print("GPU busy, retrying in 20 s", flush=True)
            time.sleep(20)
=======
        res = subprocess.run([str(HERE / "refold.sh"), f"{name}/{sid}", "1", model],
                             env={**os.environ, "DESIGN_DIR": str(d / "in")}, cwd=HERE)
>>>>>>> origin/main
        scored = HERE / d / f"scored-{model}" / "refold_metrics.csv"
        if res.returncode != 0 or not scored.is_file():
            failed.append(sid)
            print(f"FAILED: {sid} (exit {res.returncode})", flush=True)
            continue
        for row in csv.DictReader(open(scored, encoding="utf8")):
            summary.append({"sequence_id": sid, **row})
    if summary:
        out = HERE / "runs" / name / "refold_summary.csv"
        cols = list(dict.fromkeys(k for row in summary for k in row))
        with open(out, "w", newline="", encoding="utf8") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(summary)
        print(f"\nrefolded {len(summary)} of {len(rows)} designs -> {out}")
    if failed:
        print("failed:", ", ".join(failed))
    print("REFOLD_EACH_DONE")


if __name__ == "__main__":
    main()
