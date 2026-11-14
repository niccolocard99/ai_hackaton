#!/usr/bin/env bash
# Round 2 of the de novo VHH work: RFD4-Proteina ONLY (Leo's directive, 14 Nov: no Boltz, no other folding
# model, nothing but RFD4-Proteina design runs), strictly one GPU job at a time. Each step waits for the card,
# and run_design.sh keeps its own GPU guard. Run ON THE VM, detached:
#
#   setsid nohup ./queue_round2.sh > runs/round2.log 2>&1 < /dev/null &
#
# What changed after batch 1 was scored, and why (NOTES.md section 9):
#   - epitope coverage is the one axis with a usable signal (EXP-011), and only 8% of valid designs reach the
#     team's bar, so
#   1. REFINE the parents that already reach it, by partial diffusion (RFD4, antibody weights, fraction 0.4):
#      the docked pose, and with it the epitope, stays put; loops and sequence vary
#   2. a FRESH batch with CDR3 constrained to 17-21 (batch 1: short CDR3s carried HIGH estimated track risk at
#      the same epitope quality), a different seed, same hotspots
#
# The Protenix v2 refold of the batch-1 shortlist that EXP-011 recommends is NOT queued: only RFD4-Proteina is
# to be used. refold.sh and refold_each.py stay in the folder, unused.
#
# Needs LEDGER (and LEDGER_OWNER) in the environment, as run_design.sh does. A failed step does not stop the
# next one. Check the log, then score with validate_and_log.sh and epitope_check.py.
set -uo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$HERE"

echo "$(date -u +%H:%M:%S) step 1: refine strong-epitope parents (partial diffusion 0.4, 16 children each)"
PARTIAL=0.4 ./wait_and_run.sh refine_strong 16 4 || echo "step 1 ended with exit $?"

echo "$(date -u +%H:%M:%S) step 2: fresh batch, CDR3 17-21, 100 designs, seed 3"
./wait_and_run.sh fixedfw_long 100 3 || echo "step 2 ended with exit $?"

echo "$(date -u +%H:%M:%S) ROUND2_DONE"
