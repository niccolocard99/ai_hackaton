#!/usr/bin/env bash
# Refold the best designs of a run with a template-aware folding model and score them.
# Run ON THE HACKATHON VM after run_design.sh. Commands follow 18-vm.md, "Check designs by refolding".
#
#   ./refold.sh <run-name> <top-n> <of3ob|protenix-v2>
#   DESIGN_DIR=runs/shortlist1/in ./refold.sh shortlist1 12 of3ob    # refold exactly chosen designs
#       (build DESIGN_DIR with select_for_refold.py: refold-prepare itself only takes --top-n)
#
# NO BOLTZ: the team decided on 14 Nov not to use it anywhere, so the Boltz-2 option was removed.
# Use the scores as a sanity check, not a ranking gate: NICCOLO/eda/AGENT_GUIDE.md section 6 measured
# ipTM, min_ipae and pLDDT at chance for separating binders from non-binders (optimisation-track
# variants; there is no labelled de novo data). Check that the refolded complex reproduces the designed
# pose and the four hotspot contacts. Compare scores only within one folding model.
# Before refolding a table of designs, run the ledger's `check designs.csv --new-only fresh.csv` and
# refold only what is new (AGENTS.md, GENERATION LEDGER). This is not base validation: every design must
# still pass NICCOLO/scripts/check_denovo.py before it is ranked or submitted (see NOTES.md).
set -euo pipefail

NAME=${1:?usage: refold.sh <run-name> <top-n> <of3ob|protenix-v2>}
TOPN=${2:?top-n}
FM=${3:?folding model}

case "$FM" in
  of3ob)       FOLD_IMG_TAG=apheris-openfold3-v1.3.0-of3ob;    MODEL=openfold3-of3ob ;;
  protenix-v2) FOLD_IMG_TAG=apheris-protenix-v1.3.0-v2;        MODEL=protenix-v2 ;;
  *) echo "unknown folding model: $FM (of3ob or protenix-v2: template-aware; Boltz is excluded by team decision)" >&2; exit 2 ;;
esac

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$HERE"
. ./common.sh
DESIGN_DIR=${DESIGN_DIR:-"runs/$NAME/design"}
[ -d "$DESIGN_DIR" ] || { echo "no design output at $DESIGN_DIR" >&2; exit 2; }
gpu_guard

echo "=== 1/3 refold-prepare (top $TOPN) ==="
$RUN "$DESIGN_IMG" module-run refold-prepare \
  --input "$DESIGN_DIR" --output "runs/$NAME/prepared" --top-n "$TOPN"

echo "=== 2/3 fold with $MODEL ==="
FOLD_IMG=$Q:$FOLD_IMG_TAG
CMD=$(docker inspect -f '{{join .Config.Cmd " "}}' "$FOLD_IMG")
TEMPLATE=()
if [ -f "runs/$NAME/prepared/holo/template.cif" ]; then
  TEMPLATE=(--template-paths "runs/$NAME/prepared/holo/template.cif")
fi
# shellcheck disable=SC2086  # $CMD is the image's own entry command, intentionally word-split
$RUN "$FOLD_IMG" $CMD predict \
  --input "runs/$NAME/prepared/holo/queries.json" --output "runs/$NAME/$FM" \
  "${TEMPLATE[@]}" --write-detailed-scores true \
  --num-diffusion-samples 1 --seeds 42 --query-failures skip

echo "=== 3/3 refold-score ==="
$RUN "$DESIGN_IMG" module-run refold-score \
  --input "$DESIGN_DIR" --output "runs/$NAME/scored-$FM" \
  --predictions "runs/$NAME/$FM" --model "$MODEL"

echo
echo "scores: runs/$NAME/scored-$FM/refold_metrics.csv  (one row per design)"
