#!/usr/bin/env bash
# Pilot: de novo VHH, caplacizumab framework, HER2 domain IV, hotspot set B,
# framework sequence FREE (bbonly mode).
set -euo pipefail
cd /home/hackathon/hackathon-team-3/apheris/rfd4
IMG=quay.io/apheris/oslo-hackathon:apheris-rfd4-proteina-v0.3.0-antibody-2026-10-30
docker run --rm --gpus all --shm-size 16g \
  --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -v /etc/passwd:/etc/passwd:ro -v "$PWD:/work" -w /work \
  "$IMG" module-run design \
  --input runs/denovo-pilot/in/pilot.json \
  --output runs/denovo-pilot/bbonly \
  --num-replicates 2 --seed 1 --recycle-strategy high/front_loaded
