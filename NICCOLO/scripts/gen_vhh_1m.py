#!/usr/bin/env python
"""Scale ANDERS/vhh_variants.py to ~1M de novo-track VHH candidates.

Reuses ANDERS' four strategies (recombine / extend / mutate / graft) and its abguard
pre-screen unchanged; adds chunked parallel generation, global dedupe and streamed
output. Run with the abguard interpreter on the VM:

  ~/.venvs/abguard/bin/python NICCOLO/scripts/gen_vhh_1m.py RANKED.csv -o OUTDIR \
      --target 1000000 --workers 16 --seed 11

Every kept row cleared the internal margins (identity <= 0.62, CDRH3 >= 4 edits from
the reference set and trastuzumab) per abguard. Scores are inherited parent priors,
used for ordering only. A random sample is re-checked with the official check_denovo.py.
"""
from __future__ import annotations

import argparse
import csv
import math
import multiprocessing as mp
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "ANDERS"))
import vhh_variants as V  # noqa: E402  (re-execs into the abguard interpreter if needed)

G: dict = {}  # worker-shared state, set before fork


def _init(ranked_path, parents_top, identity_max, cdrh3_min):
    ranked, _, _ = V.parse_ranked(Path(ranked_path), None, None)
    V.parent_cdrs(ranked)
    top = [p for p in ranked if p.get("usable")][:parents_top]
    probe = V.prescreen([dict(r, sequence_id=f"probe_{r['sequence_id']}") for r in top],
                        identity_max, cdrh3_min)
    ok = {r["sequence_id"].removeprefix("probe_") for r in probe}
    parents = [p for p in top if p["sequence_id"] in ok] or top
    G.update(parents=parents, pssm=V.build_pssm(parents),
             identity_max=identity_max, cdrh3_min=cdrh3_min)


def _task(job: tuple[int, int, dict]) -> tuple[int, int, list[dict]]:
    tid, seed, quota = job
    rng = random.Random(seed)
    parents, pssm = G["parents"], G["pssm"]
    n = len(parents)
    cand: list[dict] = []
    if quota.get("extend"):
        cand += V.strat_extend(parents, rng, math.ceil(quota["extend"] / n), f"e{tid}")
    if quota.get("mutate"):
        cand += V.strat_mutate(parents, pssm, rng, math.ceil(quota["mutate"] / n), f"m{tid}")
    if quota.get("recombine"):
        cand += V.strat_recombine(parents, rng, quota["recombine"], f"r{tid}")
    if quota.get("graft"):
        cand += V.strat_graft(parents, rng, quota["graft"], f"g{tid}")
    seen, clean = set(), []
    for c in cand:
        s = c["heavy_sequence"]
        if (not s or s in seen or len(s) > 250 or s.count("C") % 2
                or any(ch not in "ACDEFGHIKLMNPQRSTVWY" for ch in s)):
            continue
        seen.add(s)
        clean.append(c)
    kept = V.prescreen(clean, G["identity_max"], G["cdrh3_min"])
    for r in kept:
        r.pop("seg", None)
    return tid, len(clean), kept


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("ranked")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--target", type=int, default=1_000_000)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--chunk", type=int, default=2000, help="candidates generated per batch/task")
    ap.add_argument("--parents-top", type=int, default=300)
    ap.add_argument("--identity-max", type=float, default=V.INTERNAL_IDENTITY_MAX)
    ap.add_argument("--cdrh3-dist-min", type=int, default=V.INTERNAL_CDRH3_DIST_MIN)
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    ranked, score_col, _ = V.parse_ranked(Path(a.ranked), None, None)
    V.parent_cdrs(ranked)
    top = [p for p in ranked if p.get("usable")][: a.parents_top]
    probe = V.prescreen([dict(r, sequence_id=f"probe_{r['sequence_id']}") for r in top],
                        a.identity_max, a.cdrh3_dist_min)
    ok = {r["sequence_id"].removeprefix("probe_") for r in probe}
    parents = [p for p in top if p["sequence_id"] in ok] or top
    print(f"parents: {len(parents)}/{len(top)} clear the gate with margins "
          f"(score col: {score_col})", flush=True)
    G.update(parents=parents, pssm=V.build_pssm(parents),
             identity_max=a.identity_max, cdrh3_min=a.cdrh3_dist_min)

    ch = a.chunk
    quota = {"extend": int(ch * .40), "mutate": int(ch * .40), "recombine": ch - 2 * int(ch * .40)}
    bdir = out / "batches"
    bdir.mkdir(exist_ok=True)
    seen: set[str] = set()
    n_kept = n_batches = clean_total = tid = 0
    t0 = time.time()
    SUB = ["sequence_id", "kind", "heavy_sequence", "light_sequence", "design_description"]
    META = ["sequence_id", "strategy", "parents", "parent_score", "identity_to_trastuzumab",
            "cdrh3_len", "cdrh3", "cdrh3_min_distance", "cdrh3_dist_trastuzumab",
            "predicted_track"]

    def job():
        nonlocal tid
        q = dict(quota)
        if tid == 0:
            q["graft"] = 5
        j = (tid, a.seed * 1_000_003 + tid, q)
        tid += 1
        return j

    def write_batch(rows):
        nonlocal n_kept, n_batches
        rows = rows[: a.target - n_kept]
        if not rows:
            return
        n_batches += 1
        tag = f"{n_batches:05d}"
        # write to .tmp then rename, so a batch file is never seen half-written
        with (bdir / f"batch_{tag}.csv.tmp").open("w", newline="") as fh, \
                (bdir / f"batch_{tag}_meta.csv.tmp").open("w", newline="") as mh:
            w, m = csv.writer(fh), csv.writer(mh)
            w.writerow(SUB)
            m.writerow(META)
            for i, r in enumerate(rows):
                sid = f"vhh1m_{tag}_{i:04d}"
                w.writerow([sid, V.VHH_KIND, r["heavy_sequence"], "",
                            f"{r['strategy']} variant; de novo track"])
                m.writerow([sid, r["strategy"], r["parents"], r["parent_scores"],
                            f"{r['identity_to_trastuzumab']:.4f}", r["cdrh3_len"], r["cdrh3"],
                            r["cdrh3_min_distance"], r["cdrh3_dist_trastuzumab"],
                            r["novelty_class"]])
        (bdir / f"batch_{tag}.csv.tmp").rename(bdir / f"batch_{tag}.csv")
        (bdir / f"batch_{tag}_meta.csv.tmp").rename(bdir / f"batch_{tag}_meta.csv")
        n_kept += len(rows)

    ctx = mp.get_context("spawn")  # fork deadlocks: polars threads + fork
    with ctx.Pool(a.workers, initializer=_init,
                  initargs=(str(Path(a.ranked).resolve()), a.parents_top,
                            a.identity_max, a.cdrh3_dist_min)) as pool:
        pending = [pool.apply_async(_task, (job(),)) for _ in range(a.workers * 2)]
        while pending:
            for r in list(pending):
                if not r.ready():
                    continue
                pending.remove(r)
                _t, n_clean, rows = r.get()
                clean_total += n_clean
                fresh = []
                for x in rows:
                    if x["heavy_sequence"] not in seen:
                        seen.add(x["heavy_sequence"])
                        fresh.append(x)
                write_batch(fresh)
                # keep the pipeline full until the target is covered
                if n_kept + len(pending) * ch * 0.9 < a.target:
                    pending.append(pool.apply_async(_task, (job(),)))
                print(f"[{time.time() - t0:6.0f}s] batches={n_batches} "
                      f"kept={n_kept}/{a.target} screened={clean_total}", flush=True)
            time.sleep(0.2)
    print(f"DONE: {n_kept} unique candidates in {n_batches} batches -> {bdir}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
