#!/usr/bin/env python3
"""GPU worker: folds designs with Boltz-2 whenever the shared H100 is free.

It is separate from the generation loop on purpose. Waiting for a card a teammate holds must not
stall CPU-side generation, and the loop must not be able to hold the card. The worker reads
out/pool.csv (the loop writes it atomically), builds a queue with portfolio members first, folds
in small batches so others can slot in between, and appends to out/boltz_results.csv, which the
loop merges. Create out/STOP_GPU to stop it.

Team rule (root AGENTS.md): one heavy job at a time on the shared H100. This worker only starts a
job when nvidia-smi shows the card essentially idle, so it never runs beside someone else's job.
"""
from __future__ import annotations

import sys
import time
import traceback
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, "/home/hackathon/wt/leo-design-loop/NICCOLO/eda")
import gpu_stage as gs   # noqa: E402
import run_loop as rl    # noqa: E402

OUT = rl.OUT
RES = OUT / "boltz_results.csv"
RES_X = OUT / "boltz_results_exact.csv"      # folds of exact final sequences, keyed by design id
Q_X = OUT / "fold_queue.csv"                  # written by final_portfolio.py: id, heavy_sequence, light_sequence
BATCH = 10


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def queue() -> list[str]:
    try:
        pool = rl.load_pool()
    except Exception:
        return []                      # the loop may be mid-write; try again next poll
    if pool.empty or "abg_pass" not in pool.columns:
        return []
    pool = rl.ensure_cols(pool)
    ok = pool[pool.abg_pass == True]   # noqa: E712
    done = set(pd.read_csv(RES).cdrh3) if RES.exists() else set()
    port = rl.build_portfolio(pool)
    first = port.cdrh3.tolist() if len(port) else []
    rest = ok.assign(rk=ok.groupby("arm")["R"].rank(method="first")).sort_values(["rk", "R"]).cdrh3.tolist()
    seen, q = set(), []
    for l in first + rest:
        if l not in done and l not in seen:
            seen.add(l)
            q.append(l)
    return q


def exact_queue() -> list:
    if not Q_X.exists():
        return []
    q = pd.read_csv(Q_X)
    done = set(pd.read_csv(RES_X).id) if RES_X.exists() else set()
    return [r for r in q.itertuples() if r.id not in done]


def main() -> None:
    gs.wait_for_gpu.__defaults__ = (45, 10, 3000)    # short patience per attempt; we simply poll again
    noise = None
    last_note = 0.0
    log("GPU worker started; folding only when the card is free")
    while not (OUT / "STOP_GPU").exists():
        xq = exact_queue()
        if xq:                                   # the final list goes first: exact sequences, add-ons included
            if not gs.wait_for_gpu():
                if time.time() - last_note > 300:
                    log(f"card busy ({gs.gpu_used_mb():.0f} MiB held by another job); {len(xq)} final designs queued")
                    last_note = time.time()
                continue
            try:
                if noise is None:
                    noise = gs.wt_noise_floor()
                    if not noise:
                        continue
                res = gs.run_boltz_seqs({r.id: (r.heavy_sequence, r.light_sequence) for r in xq[:BATCH]},
                                        f"bzx_{int(time.time())}", chunk=BATCH)
            except Exception:
                log("exact-sequence fold batch failed:\n" + traceback.format_exc()[-900:])
                time.sleep(30)
                continue
            if res is not None and len(res):
                floor = max(2.0, 2.0 * noise.get("h3_rmsd_pose", 0.0))
                res["pose_flag"] = res.h3_rmsd_pose > floor
                res.to_csv(RES_X, mode="a", header=not RES_X.exists(), index=False)
                log(f"folded {len(res)} exact final designs; flagged {int(res.pose_flag.sum())}")
            continue
        q = queue()
        if not q:
            time.sleep(30)
            continue
        if not gs.wait_for_gpu():
            if time.time() - last_note > 300:
                log(f"card busy ({gs.gpu_used_mb():.0f} MiB held by another job); {len(q)} designs queued")
                last_note = time.time()
            continue
        batch = q[:BATCH]
        try:
            if noise is None:
                noise = gs.wt_noise_floor()
                if not noise:
                    continue
            res = gs.run_boltz(batch, f"bzw_{int(time.time())}", chunk=BATCH)
        except Exception:
            log("fold batch failed:\n" + traceback.format_exc()[-900:])
            time.sleep(30)
            continue
        if res is None or not len(res):
            continue
        floor = max(2.0, 2.0 * noise.get("h3_rmsd_pose", 0.0))
        res["pose_flag"] = res.h3_rmsd_pose > floor
        res.to_csv(RES, mode="a", header=not RES.exists(), index=False)
        log(f"folded {len(res)}/{len(batch)}; flagged {int(res.pose_flag.sum())} (H3 more than {floor:.2f} A "
            f"from the control); about {max(0, len(q) - len(res))} left in the queue")
    log("STOP_GPU found, worker exiting")


if __name__ == "__main__":
    main()
