# Shared by run_design.sh and refold.sh. Source it; do not execute it.
# Defines Q, DESIGN_IMG, RUN, and gpu_guard. Everything here follows apheris/rfd4/18-vm.md.

Q=quay.io/apheris/oslo-hackathon
DESIGN_IMG=$Q:apheris-rfd4-proteina-v0.3.0-antibody-2026-10-30

# Options that go in front of every model command (18-vm.md, "Run a design"). Mounts the current
# directory as /work, so scripts cd to the task folder first and all outputs stay under runs/.
# --cpus 3 and the thread caps: the VM is CPU-saturated (team request, 14 Nov: stay under about 3 cores).
RUN="docker run --rm --gpus all --shm-size 16g --cpus 3 --user $(id -u):$(id -g) -e HOME=/tmp \
  -e OMP_NUM_THREADS=2 -e MKL_NUM_THREADS=2 -e OPENBLAS_NUM_THREADS=2 \
  -v /etc/passwd:/etc/passwd:ro -v $PWD:/work -w /work"

# One H100 80 GB is shared by the whole team: do not start a heavy job while another is running.
gpu_guard() {
  echo "--- GPU / containers before starting ---"
  nvidia-smi --query-gpu=name,memory.used,memory.total --format=csv,noheader || true
  docker ps --format '{{.Names}}\t{{.Image}}\t{{.RunningFor}}'
  if [ "${FORCE:-0}" != 1 ] && docker ps --format '{{.Image}}' | grep -q 'oslo-hackathon'; then
    echo "another oslo-hackathon container is running; wait for it, or re-run with FORCE=1" >&2
    exit 3
  fi
}
