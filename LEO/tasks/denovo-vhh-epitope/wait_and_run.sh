#!/usr/bin/env bash
# Queue a design run behind other RFD4 jobs: wait until no oslo-hackathon container holds the GPU, then run
# run_design.sh. The team rule is one RFD4 job on the H100 at a time (AGENTS.md; reaffirmed 14 Nov).
#
#   ./wait_and_run.sh <variant> <num-replicates> <seed> [max-wait-minutes]      (default wait: 120)
#
# run_design.sh still does its own ledger check-run / add-run and its own GPU guard, so a job that slips in
# between the last check here and the launch makes it stop rather than overlap. Never sets FORCE.
# Start it detached:  setsid nohup ./wait_and_run.sh fixedfw_long 100 3 > runs/queue.log 2>&1 < /dev/null &
set -euo pipefail

VARIANT=${1:?usage: wait_and_run.sh <variant> <num-replicates> <seed> [max-wait-minutes]}
N=${2:?num-replicates}
SEED=${3:?seed}
MAX_MIN=${4:-120}

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$HERE"

waited=0
while docker ps --format '{{.Image}}' | grep -q 'oslo-hackathon'; do
  if [ "$waited" -ge $((MAX_MIN * 60)) ]; then
    echo "gave up after $MAX_MIN minutes: the GPU is still held by another job" >&2
    exit 3
  fi
  [ $((waited % 120)) -ne 0 ] || echo "$(date -u +%H:%M:%S) waiting for: $(docker ps --format '{{.Names}}' | tr '\n' ' ')"
  sleep 10
  waited=$((waited + 10))
done
echo "$(date -u +%H:%M:%S) GPU free after ${waited}s; launching $VARIANT n=$N seed=$SEED"
exec ./run_design.sh "$VARIANT" "$N" "$SEED"
