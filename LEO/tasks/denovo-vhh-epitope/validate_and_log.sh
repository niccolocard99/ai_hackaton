#!/usr/bin/env bash
# BASE VALIDATION and generation-ledger logging for finished design runs (AGENTS.md, both mandatory).
# Run ON THE VM, after run_design.sh.
#
#   ./validate_and_log.sh TAG RUN [RUN ...]
#
# TAG names the output files (reports/TAG-*.csv); RUN is a folder under runs/. Steps:
#   1. make_designs_table.py   the BASE VALIDATION table (+ RFD4 metrics for triage)
#   2. ledger check --new-only a sequence already in the ledger is not re-validated: reuse its verdict
#   3. check_denovo.py         the assigned track must be de_novo and every shared and de novo check must pass
#      check_optimisation.py   run for the record: a de novo design is expected to FAIL it, never reclassified
#   4. ledger add, per run     logs each design with base_validation, track_assigned and cdrh3
#
# abguard and immunum are on PATH on the VM (~/.local/bin; a non-interactive shell may need the export below).
# The local track checkers measure identity to trastuzumab only: a PASS is necessary, not sufficient.
set -euo pipefail

TAG=${1:?usage: validate_and_log.sh TAG RUN [RUN ...]}
shift
[ "$#" -ge 1 ] || { echo "give at least one RUN" >&2; exit 2; }

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$HERE"
export PATH="$HOME/.local/bin:$PATH"
LEDGER=${LEDGER:-"$HERE/../generation-ledger/ledger.py"}
CHECKERS=${CHECKERS:-"$HERE/../../../NICCOLO/scripts"}
PROVIDED=${PROVIDED:-"$HERE/../../../apheris/hackathon-data/her2-trastuzumab/trastuzumab_her2_spr_all_768.csv"}
[ -f "$LEDGER" ] || { echo "generation ledger not found at $LEDGER (set LEDGER=...)" >&2; exit 2; }
[ -f "$CHECKERS/check_denovo.py" ] || { echo "track checkers not found at $CHECKERS (set CHECKERS=...)" >&2; exit 2; }
mkdir -p reports

echo "=== 1. table ==="
python3 make_designs_table.py "$@" -o "reports/$TAG-designs.csv" --metrics "reports/$TAG-metrics.csv"

echo "=== 2. ledger: which sequences are new? ==="
set +e
python3 "$LEDGER" check "reports/$TAG-designs.csv" --new-only "reports/$TAG-fresh.csv"
rc=$?
set -e
[ "$rc" -eq 0 ] || [ "$rc" -eq 3 ] || exit "$rc"
total=$(($(wc -l < "reports/$TAG-designs.csv") - 1))
fresh=$(($(wc -l < "reports/$TAG-fresh.csv") - 1))
echo "$fresh of $total designs are new to the ledger"
[ "$fresh" -eq "$total" ] || echo "WARNING: $((total - fresh)) already logged; their earlier verdicts apply and they are skipped below" >&2

echo "=== 3. base validation (de novo track) ==="
python3 "$CHECKERS/check_denovo.py" "reports/$TAG-fresh.csv" -o "reports/$TAG-denovo.csv" 2>&1 | grep -v "DEBUG" | tail -n 3 || true
echo "=== 3b. optimisation checker, for the record ==="
python3 "$CHECKERS/check_optimisation.py" "reports/$TAG-fresh.csv" -o "reports/$TAG-optimisation.csv" --provided "$PROVIDED" 2>&1 | grep -v "DEBUG" | tail -n 1 || true

echo "=== 4. log to the ledger, per run ==="
python3 - "$TAG" "$@" <<'EOF'
import csv, sys
tag, runs = sys.argv[1], sys.argv[2:]
designs = {r["sequence_id"]: r for r in csv.DictReader(open(f"reports/{tag}-fresh.csv", encoding="utf8"))}
report = {r["sequence_id"]: r for r in csv.DictReader(open(f"reports/{tag}-denovo.csv", encoding="utf8"))}
for run in runs:
    rows = []
    for sid, d in designs.items():
        if sid.split("__")[0] != run:
            continue
        r = report[sid]
        rows.append({**d, "cdrh3": r["cdrh3"], "track_assigned": r["track_assigned"],
                     "base_validation": "PASS" if r["passes"].lower() == "true" else "FAIL"})
    if rows:
        with open(f"reports/{tag}-{run}-ledger.csv", "w", newline="", encoding="utf8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"{run}: {len(rows)} designs, {sum(r['base_validation'] == 'PASS' for r in rows)} PASS")
EOF
for run in "$@"; do
  f="reports/$TAG-$run-ledger.csv"
  [ -f "$f" ] && python3 "$LEDGER" add "$f" --run-id "$run"
done
echo "done: reports/$TAG-denovo.csv, reports/$TAG-metrics.csv"
