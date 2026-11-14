#!/usr/bin/env bash
# Design VHHs against the trastuzumab epitope on HER2 domain IV (RFD4-Proteina, antibody weights).
# Run ON THE HACKATHON VM (needs Docker and the H100). Not runnable on a laptop.
#
#   ./run_design.sh <fixedfw|freefw> <num-replicates> <seed> [run-name]
#
# Order of operations (AGENTS.md, GENERATION LEDGER, mandatory):
#   1. check_spec.py        the spec must resolve against its structures
#   2. ledger check-run     exit 3 = this exact run was already launched: stop and reuse its output
#   3. GPU guard            one H100 is shared; refuse while another oslo-hackathon container runs
#   4. ledger add-run       registered BEFORE the job starts, so a teammate sees it
#   5. docker run           the design job
#   6. ledger set-status    done / failed / empty
#
# Smoke test first (N=2), read the output, then scale. Run it in the background and poll the log:
#   nohup ./run_design.sh fixedfw 2 1 > /dev/null 2>&1 &
# Outputs go to runs/<run-name>/ (git-ignored). Override the ledger location with LEDGER=/path/ledger.py.
set -euo pipefail

VARIANT=${1:?usage: run_design.sh <fixedfw|freefw> <num-replicates> <seed> [run-name]}
N=${2:?num-replicates}
SEED=${3:?seed}
NAME=${4:-vhh-${VARIANT}-n${N}-s${SEED}}

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$HERE"
. ./common.sh

SPEC=spec/design_${VARIANT}.json
[ -f "$SPEC" ] || { echo "no such spec: $SPEC" >&2; exit 2; }
[ ! -e "runs/$NAME" ] || { echo "runs/$NAME already exists; choose another run name" >&2; exit 2; }

LEDGER=${LEDGER:-"$HERE/../generation-ledger/ledger.py"}
[ -f "$LEDGER" ] || { echo "generation ledger not found at $LEDGER; AGENTS.md makes it mandatory before every run (set LEDGER=...)" >&2; exit 2; }

<<<<<<< HEAD
# PARTIAL=0.4 switches to refinement: re-run that fraction of the sampling trajectory on the parent CIFs listed
# in the spec (make_refine_spec.py), instead of designing from scratch. The parents carry their own conditions.
PARTIAL=${PARTIAL:-}
if [ -n "$PARTIAL" ]; then
  INPUTS=(--input "$SPEC"); for f in spec/refine/*.cif; do INPUTS+=(--input "$f"); done
  EXTRA=(--partial-trajectory-fraction "$PARTIAL")
  PARTIAL_SETTING=(--setting "partial-trajectory-fraction=$PARTIAL")
else
  INPUTS=(--input "$SPEC" --input spec/target_her2_d4.cif --input spec/framework_vhh.cif)
  EXTRA=()
  PARTIAL_SETTING=()
fi

# Every flag that changes what is generated goes to the ledger, spelled as in the command without dashes.
# --input is passed per file: the ledger only follows `input:` inside .yaml specs, and ours is JSON.
LEDGER_FLAGS=(--method design --model "$DESIGN_IMG" "${INPUTS[@]}"
  --seed "$SEED" --setting "num-replicates=$N" --setting "batch-size=1"
  --setting "recycle-strategy=high/front_loaded" "${PARTIAL_SETTING[@]}")

# 1. Fail on a malformed spec before touching the GPU or the ledger. (A refine spec has no contig or conditions
#    to check: the parents carry them, and RFD4 prevalidates the inputs itself.)
if [ -z "$PARTIAL" ]; then
  python3 check_spec.py "$SPEC" > /dev/null || { echo "check_spec.py failed for $SPEC" >&2; exit 2; }
fi
=======
# Every flag that changes what is generated goes to the ledger, spelled as in the command without dashes.
# --input is passed per file: the ledger only follows `input:` inside .yaml specs, and ours is JSON.
LEDGER_FLAGS=(--method design --model "$DESIGN_IMG"
  --input "$SPEC" --input spec/target_her2_d4.cif --input spec/framework_vhh.cif
  --seed "$SEED" --setting "num-replicates=$N" --setting "batch-size=1"
  --setting "recycle-strategy=high/front_loaded")

# 1. Fail on a malformed spec before touching the GPU or the ledger.
python3 check_spec.py "$SPEC" > /dev/null || { echo "check_spec.py failed for $SPEC" >&2; exit 2; }
>>>>>>> origin/main

# 2. Already run? (exit 3). Capture the code: set -e would otherwise hide the ledger's own message.
set +e
python3 "$LEDGER" check-run "${LEDGER_FLAGS[@]}"
rc=$?
set -e
if [ "$rc" -eq 3 ]; then echo "not launching: reuse the output listed above" >&2; exit 3; fi
[ "$rc" -eq 0 ] || exit "$rc"

# 3. Shared GPU.
gpu_guard

# 4. Register before launching.
mkdir -p "runs/$NAME"
python3 "$LEDGER" add-run "${LEDGER_FLAGS[@]}" --run-id "$NAME" --output "runs/$NAME"

# 5. Pass the spec FILE as --input: a parent folder silently ignores contig and conditions.
# batch-size 1 for the longest input avoids running out of GPU memory; high/front_loaded is for antibodies.
set +e
$RUN "$DESIGN_IMG" module-run design \
  --input "$SPEC" --output "runs/$NAME/design" \
<<<<<<< HEAD
  --num-replicates "$N" --seed "$SEED" --batch-size 1 --recycle-strategy high/front_loaded "${EXTRA[@]}" \
=======
  --num-replicates "$N" --seed "$SEED" --batch-size 1 --recycle-strategy high/front_loaded \
>>>>>>> origin/main
  2>&1 | tee "runs/$NAME/design.log"
rc=${PIPESTATUS[0]}
set -e

# 6. Record the outcome.
status=failed
if [ "$rc" -eq 0 ]; then
  n_out=$(ls runs/"$NAME"/design/generation/inference_samples/*/*_generated.cif 2>/dev/null | wc -l)
  if [ "$n_out" -gt 0 ]; then status=done; else status=empty; fi
fi
python3 "$LEDGER" set-status --run-id "$NAME" --status "$status" --notes "docker exit $rc"

echo
echo "run $NAME: $status"
echo "designs : runs/$NAME/design/generation/inference_samples/*/*_generated.cif"
echo "metrics : runs/$NAME/design/generation/inference_metrics/per_example_metrics.csv"
echo "next    : build the BASE VALIDATION table, then 'ledger.py add designs.csv --run-id $NAME'"
[ "$status" = done ]
