"""What the machine looks like right now: GPU, heavy containers, disk, load, JapanFold credit.

Every probe degrades to `None`/`False` when its tool is missing or slow, so the loop never stops
for want of telemetry. The JapanFold key is read from the environment and never written anywhere;
only the numeric balance is recorded.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

JAPANFOLD_BASE = "https://api.japanfold.aiand.com"
HEAVY_IMAGE_HINTS = ("apheris", "boltz", "protenix", "openfold", "esmfold", "rfdiffusion", "alphafold", "chai")


def _run(cmd: list, timeout_s: float = 10) -> str | None:
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s)
        return p.stdout if p.returncode == 0 else None
    except Exception:
        return None


def parse_gpu(stdout: str | None) -> dict | None:
    """`nvidia-smi --query-gpu=...` csv -> {util_pct, mem_used_mb, mem_total_mb} for the first GPU."""
    if not stdout or not stdout.strip():
        return None
    line = stdout.strip().splitlines()[0]
    try:
        util, used, total = (int(round(float(x))) for x in line.split(","))
    except ValueError:
        return None
    return {"util_pct": util, "mem_used_mb": used, "mem_total_mb": total}


def parse_docker_ps(stdout: str | None) -> dict:
    """`docker ps --format '{{.Image}}'` -> container images running now, plus heavy-image detection."""
    images = [ln.strip() for ln in (stdout or "").splitlines() if ln.strip()]
    heavy = sorted({img for img in images if any(h in img.lower() for h in HEAVY_IMAGE_HINTS)})
    return {"containers": images, "heavy_images": heavy, "heavy_job_running": bool(heavy)}


def japanfold_balance(timeout_s: float = 8) -> float | None:
    """GET /v1/credits through the JapanFold key, if one is configured. Returns dollars, or None."""
    key = os.environ.get("JAPANFOLD_API_KEY")
    if not key:
        return None
    try:
        req = urllib.request.Request(f"{JAPANFOLD_BASE}/v1/credits",
                                     headers={"Authorization": f"Bearer {key}"})
        with urllib.request.urlopen(req, timeout=timeout_s) as r:
            body = json.loads(r.read().decode())
        for k in ("balance", "credits", "amount"):
            if isinstance(body, dict) and isinstance(body.get(k), (int, float)):
                return float(body[k])
        if isinstance(body, dict) and isinstance(body.get("data"), dict):
            for k in ("balance", "credits", "amount"):
                if isinstance(body["data"].get(k), (int, float)):
                    return float(body["data"][k])
        return None
    except Exception:
        return None


def probe(out_dir: Path | None = None) -> dict:
    """One snapshot. Cheap (sub-second except the optional JapanFold call); safe to run every round."""
    gpu_raw = None
    if shutil.which("nvidia-smi"):
        gpu_raw = _run(["nvidia-smi", "--query-gpu=utilization.gpu,memory.used,memory.total",
                        "--format=csv,noheader,nounits"])
    docker_raw = _run(["docker", "ps", "--format", "{{.Image}}"]) if shutil.which("docker") else None
    free_gb = None
    try:
        import shutil as sh
        usage = sh.disk_usage(str(out_dir or Path.cwd()))
        free_gb = round(usage.free / 1e9, 1)
    except Exception:
        pass
    try:
        load1 = round(os.getloadavg()[0], 2)
    except Exception:
        load1 = None
    return {"utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "gpu": parse_gpu(gpu_raw), **parse_docker_ps(docker_raw),
            "disk_free_gb": free_gb, "load1": load1,
            "japanfold_balance_usd": japanfold_balance(),
            "amass_key_present": bool(os.environ.get("AMASS_API_KEY"))}


def seconds_until(deadline) -> float | None:
    if deadline is None:
        return None
    return max(0.0, deadline.timestamp() - time.time())
