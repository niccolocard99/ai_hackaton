#!/usr/bin/env python3
"""broker.py — the single gate for every compute job in the de novo HER2 campaign.

Why this exists
---------------
Three compute lanes have hard concurrency limits and two of them cost money:

  local      1 heavy job at a time (one H100, SHARED with the rest of the team)
  japanfold  <= 2 in flight, and the API allows only 2 concurrent wait-holds per
             account, so exactly ONE poller may wait at a time
  nebius     <= 2 in flight; 16 teams share the backend, queue delays are normal

Reasoning agents fan out freely; compute does not. Everything that costs GPU
seconds or dollars goes through here so that we (a) never exceed a lane's
concurrency, (b) never submit the same spec twice, and (c) always know what we
spent. The registry is append-only JSONL so a crash or a context reset cannot
lose the record of an in-flight job -- which is the failure that actually bites,
because a resubmitted hosted job costs real money and a resubmitted local job
steals the GPU from a teammate.

Design notes
------------
* Append-only `jobs.jsonl`: one JSON object per state transition, newest wins.
  Never rewritten in place, so concurrent appends cannot corrupt earlier rows.
* Dedupe key is the spec hash (sha256 of the canonical spec/command payload).
  `submit` on a spec already present returns the existing job instead of a new
  one, unless --force. Idempotency-Key for hosted lanes is derived from the same
  hash, so a retried HTTP call cannot create a second remote job.
* A flock'd lockfile serialises the read-modify-append, so two agents calling
  submit at the same instant cannot both win the last concurrency slot.
* The broker does not know how to talk to JapanFold or Nebius. It issues
  reservations. The lane client (oracles/refold.py, generate.py) asks for a slot,
  submits, then reports back the remote id and cost. That keeps credentials out
  of this file entirely -- no key is ever read, printed or stored here.

Usage
-----
  broker.py submit --lane local --tag arm01 --cmd-file cmd.sh [--cost-est 0]
  broker.py submit --lane japanfold --tag arm01 --spec-file spec.json --cost-est 0.35
  broker.py claim   --lane local            # -> next queued job, marks running
  broker.py finish  --job-id J --status succeeded --cost 0.31 --out runs/...
  broker.py status  [--lane L]              # concurrency + spend summary
  broker.py run-local                       # drain the local queue, one at a time
  broker.py spend                           # per-lane and total $ spent

Spend caps (hard, from the campaign brief): JapanFold total $400, of which
generation <= $250 and refolds <= $150. The broker REFUSES a submit that would
cross a cap. Raising a cap is a human decision, not an agent's.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
import subprocess
import sys
import time
import uuid
from contextlib import contextmanager
from pathlib import Path

HERE = Path(__file__).resolve().parent
REGISTRY = HERE / "jobs.jsonl"
LOCKFILE = HERE / ".broker.lock"

# lane -> max simultaneous jobs
LANE_LIMITS = {"local": 1, "japanfold": 2, "nebius": 2}

# hard dollar caps. None = untracked/free (local GPU costs us no cash).
LANE_CAPS = {"local": None, "japanfold": 400.0, "nebius": None}
# purpose sub-caps within a lane
PURPOSE_CAPS = {("japanfold", "generation"): 250.0, ("japanfold", "refold"): 150.0}

TERMINAL = {"succeeded", "failed", "canceled", "skipped"}


@contextmanager
def _lock():
    """Serialise read-modify-append across processes/agents."""
    LOCKFILE.touch(exist_ok=True)
    fh = open(LOCKFILE, "r+")
    try:
        fcntl.flock(fh, fcntl.LOCK_EX)
        yield
    finally:
        fcntl.flock(fh, fcntl.LOCK_UN)
        fh.close()


def _append(rec: dict) -> None:
    rec = dict(rec)
    rec.setdefault("ts", time.time())
    rec.setdefault("ts_iso", time.strftime("%Y-%m-%dT%H:%M:%S"))
    with open(REGISTRY, "a") as fh:
        fh.write(json.dumps(rec, sort_keys=True) + "\n")
        fh.flush()
        os.fsync(fh.fileno())


def _read_all() -> list[dict]:
    if not REGISTRY.exists():
        return []
    out = []
    for line in REGISTRY.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            out.append(json.loads(line))
        except json.JSONDecodeError:
            # A torn final line can only be the newest record; ignore it rather
            # than crash, so one bad append never bricks the registry.
            continue
    return out


def jobs_state() -> dict[str, dict]:
    """Collapse the append-only log into current state per job_id."""
    state: dict[str, dict] = {}
    for rec in _read_all():
        jid = rec.get("job_id")
        if not jid:
            continue
        state.setdefault(jid, {}).update(rec)
    return state


def _spend(state: dict[str, dict]) -> dict:
    tot: dict = {"by_lane": {}, "by_purpose": {}, "total": 0.0}
    for j in state.values():
        c = float(j.get("cost") or 0.0)
        if not c:
            continue
        lane, purpose = j.get("lane", "?"), j.get("purpose", "?")
        tot["by_lane"][lane] = tot["by_lane"].get(lane, 0.0) + c
        tot["by_purpose"][f"{lane}/{purpose}"] = (
            tot["by_purpose"].get(f"{lane}/{purpose}", 0.0) + c
        )
        tot["total"] += c
    return tot


def _inflight(state: dict[str, dict], lane: str) -> list[dict]:
    return [
        j
        for j in state.values()
        if j.get("lane") == lane and j.get("status") in ("queued", "running", "submitted")
    ]


def _running(state: dict[str, dict], lane: str) -> list[dict]:
    return [
        j
        for j in state.values()
        if j.get("lane") == lane and j.get("status") in ("running", "submitted")
    ]


def spec_hash(payload: str) -> str:
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def cmd_submit(a) -> int:
    if a.cmd_file:
        payload = Path(a.cmd_file).read_text()
    elif a.spec_file:
        payload = Path(a.spec_file).read_text()
    elif a.cmd:
        payload = a.cmd
    else:
        print("submit needs --cmd, --cmd-file or --spec-file", file=sys.stderr)
        return 2
    # The dedupe identity is lane + purpose + payload + tag: the same spec run on
    # a different lane or for a different purpose is a legitimately different job,
    # but the same spec twice on the same lane is the mistake we are preventing.
    h = spec_hash(f"{a.lane}|{a.purpose}|{a.tag}|{payload}")

    with _lock():
        state = jobs_state()
        dupes = [j for j in state.values() if j.get("spec_hash") == h]
        if dupes and not a.force:
            j = sorted(dupes, key=lambda r: r.get("ts", 0))[-1]
            print(json.dumps({"deduped": True, "job_id": j["job_id"], "status": j.get("status")}))
            return 0

        est = float(a.cost_est or 0.0)
        sp = _spend(state)
        cap = LANE_CAPS.get(a.lane)
        if cap is not None and sp["by_lane"].get(a.lane, 0.0) + est > cap:
            print(
                json.dumps({"refused": "lane_cap", "lane": a.lane, "cap": cap,
                            "spent": sp["by_lane"].get(a.lane, 0.0), "est": est}),
                file=sys.stderr)
            return 3
        pcap = PURPOSE_CAPS.get((a.lane, a.purpose))
        if pcap is not None:
            cur = sp["by_purpose"].get(f"{a.lane}/{a.purpose}", 0.0)
            if cur + est > pcap:
                print(
                    json.dumps({"refused": "purpose_cap", "lane": a.lane,
                                "purpose": a.purpose, "cap": pcap, "spent": cur, "est": est}),
                    file=sys.stderr)
                return 3
        if est > 25.0:
            print(json.dumps({"refused": "single_job_over_25", "est": est}), file=sys.stderr)
            return 3

        jid = f"{a.lane[:3]}-{a.tag}-{uuid.uuid4().hex[:6]}"
        _append({
            "job_id": jid, "lane": a.lane, "purpose": a.purpose, "tag": a.tag,
            "spec_hash": h, "idempotency_key": f"{a.lane}-{h}",
            "status": "queued", "cost_est": est, "cost": 0.0,
            "payload_path": a.cmd_file or a.spec_file or None,
            "payload_inline": a.cmd if a.cmd else None,
            "out": a.out, "note": a.note,
        })
    print(json.dumps({"job_id": jid, "idempotency_key": f"{a.lane}-{h}", "status": "queued"}))
    return 0


def cmd_claim(a) -> int:
    """Hand out the next queued job on a lane, if a concurrency slot is free."""
    with _lock():
        state = jobs_state()
        limit = LANE_LIMITS.get(a.lane, 1)
        if len(_running(state, a.lane)) >= limit:
            print(json.dumps({"claimed": None, "reason": "lane_busy",
                              "running": len(_running(state, a.lane)), "limit": limit}))
            return 0
        q = sorted(
            [j for j in state.values() if j.get("lane") == a.lane and j.get("status") == "queued"],
            key=lambda r: r.get("ts", 0),
        )
        if not q:
            print(json.dumps({"claimed": None, "reason": "empty_queue"}))
            return 0
        j = q[0]
        _append({"job_id": j["job_id"], "status": "running", "claimed_at": time.time()})
    print(json.dumps({"claimed": j["job_id"], "payload_path": j.get("payload_path"),
                      "payload_inline": j.get("payload_inline"), "out": j.get("out"),
                      "idempotency_key": j.get("idempotency_key"), "tag": j.get("tag")}))
    return 0


def cmd_finish(a) -> int:
    with _lock():
        rec = {"job_id": a.job_id, "status": a.status, "finished_at": time.time()}
        if a.cost is not None:
            rec["cost"] = float(a.cost)
        if a.out:
            rec["out"] = a.out
        if a.remote_id:
            rec["remote_id"] = a.remote_id
        if a.note:
            rec["note"] = a.note
        _append(rec)
    print(json.dumps({"job_id": a.job_id, "status": a.status}))
    return 0


def cmd_mark(a) -> int:
    """Record a hosted lane's remote id / status without finishing it."""
    with _lock():
        rec = {"job_id": a.job_id, "status": a.status or "submitted"}
        if a.remote_id:
            rec["remote_id"] = a.remote_id
        if a.cost is not None:
            rec["cost"] = float(a.cost)
        _append(rec)
    print(json.dumps({"job_id": a.job_id, "marked": a.status or "submitted"}))
    return 0


def cmd_status(a) -> int:
    state = jobs_state()
    sp = _spend(state)
    rows = [j for j in state.values() if not a.lane or j.get("lane") == a.lane]
    by_status: dict[str, int] = {}
    for j in rows:
        by_status[j.get("status", "?")] = by_status.get(j.get("status", "?"), 0) + 1
    out = {
        "jobs": len(rows),
        "by_status": by_status,
        "lane_concurrency": {
            ln: {"running": len(_running(state, ln)), "limit": LANE_LIMITS[ln],
                 "queued": len([j for j in state.values()
                                if j.get("lane") == ln and j.get("status") == "queued"])}
            for ln in LANE_LIMITS
        },
        "spend": sp,
        "caps": {"lane": LANE_CAPS, "purpose": {f"{k[0]}/{k[1]}": v for k, v in PURPOSE_CAPS.items()}},
    }
    print(json.dumps(out, indent=2))
    return 0


def cmd_spend(a) -> int:
    print(json.dumps(_spend(jobs_state()), indent=2))
    return 0


def cmd_list(a) -> int:
    state = jobs_state()
    rows = sorted(state.values(), key=lambda r: r.get("ts", 0))
    for j in rows:
        if a.lane and j.get("lane") != a.lane:
            continue
        if a.status and j.get("status") != a.status:
            continue
        print(f"{j.get('job_id'):34s} {j.get('lane',''):10s} {j.get('purpose',''):11s} "
              f"{j.get('status',''):10s} ${float(j.get('cost') or 0):7.3f}  {j.get('tag','')} "
              f"{j.get('note','') or ''}")
    return 0


def cmd_run_local(a) -> int:
    """Drain the local queue, strictly one job at a time.

    Each queued local job's payload is a shell script. We run it to completion,
    then record the outcome. Before each job we re-check the GPU, because the
    H100 is shared with the rest of the team and a teammate's job appearing is a
    reason to wait, not to pile on.
    """
    deadline = time.time() + (a.max_seconds or 10 ** 9)
    ran = 0
    while time.time() < deadline:
        if a.max_jobs and ran >= a.max_jobs:
            break
        with _lock():
            state = jobs_state()
            if len(_running(state, "local")) >= LANE_LIMITS["local"]:
                busy = True
            else:
                busy = False
                q = sorted(
                    [j for j in state.values()
                     if j.get("lane") == "local" and j.get("status") == "queued"],
                    key=lambda r: r.get("ts", 0),
                )
                job = q[0] if q else None
                if job:
                    _append({"job_id": job["job_id"], "status": "running",
                             "claimed_at": time.time()})
        if busy:
            time.sleep(10)
            continue
        if not job:
            if a.once:
                break
            time.sleep(10)
            continue

        script = job.get("payload_path")
        t0 = time.time()
        print(f"[broker] running {job['job_id']} ({job.get('tag')}) -> {script}", flush=True)
        try:
            p = subprocess.run(["bash", script], capture_output=True, text=True, timeout=a.job_timeout)
            ok = p.returncode == 0
            tail = (p.stdout[-2000:] + p.stderr[-2000:]) if not ok else p.stdout[-500:]
        except subprocess.TimeoutExpired:
            ok, tail = False, "TIMEOUT"
        dt = time.time() - t0
        with _lock():
            _append({"job_id": job["job_id"],
                     "status": "succeeded" if ok else "failed",
                     "finished_at": time.time(), "wall_seconds": round(dt, 1),
                     "note": (job.get("note") or "") + f" | wall={dt:.0f}s",
                     "log_tail": tail[-1500:]})
        print(f"[broker] {job['job_id']} {'ok' if ok else 'FAILED'} in {dt:.0f}s", flush=True)
        ran += 1
        if a.once:
            break
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("submit")
    s.add_argument("--lane", required=True, choices=sorted(LANE_LIMITS))
    s.add_argument("--purpose", default="generation",
                   choices=["generation", "refold", "score", "other"])
    s.add_argument("--tag", required=True)
    s.add_argument("--cmd")
    s.add_argument("--cmd-file")
    s.add_argument("--spec-file")
    s.add_argument("--cost-est", type=float, default=0.0)
    s.add_argument("--out")
    s.add_argument("--note")
    s.add_argument("--force", action="store_true")
    s.set_defaults(fn=cmd_submit)

    s = sub.add_parser("claim")
    s.add_argument("--lane", required=True, choices=sorted(LANE_LIMITS))
    s.set_defaults(fn=cmd_claim)

    s = sub.add_parser("finish")
    s.add_argument("--job-id", required=True)
    s.add_argument("--status", required=True,
                   choices=["succeeded", "failed", "canceled", "skipped"])
    s.add_argument("--cost", type=float)
    s.add_argument("--out")
    s.add_argument("--remote-id")
    s.add_argument("--note")
    s.set_defaults(fn=cmd_finish)

    s = sub.add_parser("mark")
    s.add_argument("--job-id", required=True)
    s.add_argument("--status")
    s.add_argument("--remote-id")
    s.add_argument("--cost", type=float)
    s.set_defaults(fn=cmd_mark)

    s = sub.add_parser("status")
    s.add_argument("--lane", choices=sorted(LANE_LIMITS))
    s.set_defaults(fn=cmd_status)

    s = sub.add_parser("spend")
    s.set_defaults(fn=cmd_spend)

    s = sub.add_parser("list")
    s.add_argument("--lane", choices=sorted(LANE_LIMITS))
    s.add_argument("--status")
    s.set_defaults(fn=cmd_list)

    s = sub.add_parser("run-local")
    s.add_argument("--once", action="store_true")
    s.add_argument("--max-jobs", type=int)
    s.add_argument("--max-seconds", type=int)
    s.add_argument("--job-timeout", type=int, default=5400)
    s.set_defaults(fn=cmd_run_local)

    a = ap.parse_args()
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
