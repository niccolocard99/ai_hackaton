#!/usr/bin/env python3
"""The loop, run as a background process on the VM.

  generate -> gates -> classifier -> ledger check -> real abguard -> H100 -> pool -> portfolio

Each round is one generation run in the GENERATION LEDGER (method prefilter, as the rule
asks for a screen that scores millions): checked before launch, registered before it
starts, its shortlist logged after. State lives in out/, so the loop resumes where it
stopped. Create out/STOP to finish the current round and exit.

What the loop can and cannot do. Both the classifier and the ranker are additive over
positions, so their optimum is reachable in one shot (Jahn pipeline LOOP.md says the same
and an independent implementation lands on the same motif). Extra rounds therefore do not
make the surrogate score better; they widen coverage, push far more candidates through the
real abguard, and bring in the scorers the surrogate lacks: epitope contact (8 structures),
Tm (NanoMelt) and structure plausibility (Boltz-2), all on the H100.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, "/home/hackathon/wt/leo-design-loop/NICCOLO/eda")

import designloop as dl            # noqa: E402
import gpu_stage as gs             # noqa: E402
import validate_abguard as vab     # noqa: E402
from her2lib import TZ_VL, vh_with_loop   # noqa: E402
from novelty import encode         # noqa: E402

OUT = Path(os.environ.get("DL_OUT", HERE / "out"))          # tests point these at scratch folders
LED = Path(os.environ.get("DL_LEDGER", HERE / "ledger"))
dl.OUT = OUT
POOL = OUT / "pool.csv"
SPR = HERE / "data" / "trastuzumab_her2_spr_all_768.csv"
FACS = HERE / "ref" / "tz_her2_facs_3classes_524346.csv"
EDA = Path("/home/hackathon/wt/leo-design-loop/NICCOLO/eda")
JAHN_RUNS = Path("/home/hackathon/wt/jahn/JAHN/pipeline/runs")
MODEL = "prefilter-niccolo-eda-with-leo-design-loop-v1"   # no slash: the ledger keeps the last path segment

ABG_COLS = ["abg_pass", "abg_track", "abg_cdrh3", "abg_min_dist", "abg_ident", "abg_failures"]
NM_COLS = ["tm_pred", "dtm_vs_wt"]
FOLD_COLS = ["iptm", "ptm", "plddt", "gpde", "ab_her2_iptm", "vh_rmsd", "vl_rmsd",
             "h3_rmsd_pose", "h3_rmsd_shape", "pose_flag", "model_cif"]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# --------------------------------------------------------------------------- the ledger
def led(*args: str) -> tuple[int, str]:
    p = subprocess.run([sys.executable, str(LED / "ledger.py"), *args], capture_output=True, text=True)
    return p.returncode, (p.stdout + p.stderr).strip()


def run_flags(cfg: dict, boot_seed: int, parents_hash: str = "") -> list[str]:
    """Every file the round reads and every setting that changes what it generates."""
    f = ["--method", "prefilter", "--model", MODEL, "--seed", str(boot_seed)]
    for path in (SPR, FACS, EDA / "prefilter.py", HERE / "designloop.py", HERE / "models.py",
                 HERE / "novelty.py", HERE / "her2lib.py"):
        f += ["--input", str(path)]
    st = {k: (",".join(map(str, v)) if isinstance(v, (list, tuple)) else v)
          for k, v in cfg.items() if k != "extra_arms"}
    st.update(min_sub=dl.MIN_SUB, d_max=dl.D_MAX, p_thr=dl.P_THR, k_lcb=dl.K_LCB,
              min_support=dl.MIN_SUPPORT)
    for name, loops in (cfg.get("extra_arms") or {}).items():
        st[f"extra_{name}"] = hashlib.sha256(",".join(sorted(loops)).encode()).hexdigest()[:10]
    if parents_hash:
        st["refine_parents"] = parents_hash
    for k, v in sorted(st.items()):
        f += ["--setting", f"{k}={v}"]
    return f


def validation_table(rows: pd.DataFrame) -> pd.DataFrame:
    """The BASE VALIDATION table the ledger takes; verdicts only where abguard has run."""
    out = []
    for _, r in rows.iterrows():
        row = dict(sequence_id=vab.design_id(r.cdrh3), kind="paired",
                   heavy_sequence=vh_with_loop(r.cdrh3), light_sequence=TZ_VL,
                   design_description=f"CDRH3 {r.cdrh3} ({r.subs}); arm {r.arm}; round {int(r['round'])}",
                   cdrh3="", track_assigned="", base_validation="")
        if "abg_pass" in r and pd.notna(r.abg_pass):
            row.update(cdrh3=r.abg_cdrh3, track_assigned=r.abg_track,
                       base_validation="PASS" if bool(r.abg_pass) else "FAIL")
        out.append(row)
    return pd.DataFrame(out)


# --------------------------------------------------------------------------- pool
def load_pool() -> pd.DataFrame:
    return pd.read_csv(POOL) if POOL.exists() else pd.DataFrame()


def ensure_cols(pool: pd.DataFrame) -> pd.DataFrame:
    for c in ABG_COLS + NM_COLS + FOLD_COLS:
        if c not in pool.columns:
            pool[c] = np.nan
    return pool


def merge_shortlist(pool: pd.DataFrame, sl: pd.DataFrame, rnd: int) -> tuple[pd.DataFrame, int]:
    """Add new designs; for ones already known keep their results and union the arms."""
    sl = sl.copy()
    sl["first_round"] = rnd
    if pool.empty:
        return ensure_cols(sl), len(sl)
    known = set(pool.cdrh3)
    new = sl[~sl.cdrh3.isin(known)]
    old = sl[sl.cdrh3.isin(known)].set_index("cdrh3").arm
    pool = pool.set_index("cdrh3", drop=False)
    for l, a in old.items():
        merged = "+".join(sorted(set(str(pool.at[l, "arm"]).split("+")) | set(a.split("+"))))
        pool.at[l, "arm"] = merged
    pool = pd.concat([pool.reset_index(drop=True), new], ignore_index=True)
    return ensure_cols(pool), len(new)


def put(pool: pd.DataFrame, df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    """Write df[cols] (keyed by cdrh3) into the pool rows with that cdrh3."""
    pool = pool.set_index("cdrh3", drop=False)
    for c in cols:
        if c not in pool.columns:
            pool[c] = np.nan
        pool[c] = pool[c].astype(object)
    d = df.set_index("cdrh3")
    for c in cols:
        pool.loc[d.index, c] = d[c].astype(object).values
    return pool.reset_index(drop=True)


def save_pool(pool: pd.DataFrame) -> None:
    tmp = POOL.with_suffix(".tmp")
    pool.to_csv(tmp, index=False)
    os.replace(tmp, POOL)           # atomic: the GPU worker reads this file while the loop runs


def merge_gpu(pool: pd.DataFrame) -> pd.DataFrame:
    """Pick up whatever the GPU worker has folded so far."""
    f = OUT / "boltz_results.csv"
    if pool.empty or not f.exists():
        return pool
    r = pd.read_csv(f).drop_duplicates("cdrh3", keep="last")
    r = r[r.cdrh3.isin(set(pool.cdrh3))]
    return put(pool, r, FOLD_COLS) if len(r) else pool


def abguard_cache() -> pd.DataFrame:
    """Real abguard verdicts from earlier rounds: a design already validated is not validated again."""
    parts = []
    for f in OUT.glob("round_*/abguard*/report_*.csv"):
        try:
            parts.append(pd.read_csv(f))
        except Exception:
            pass
    if not parts:
        return pd.DataFrame()
    d = pd.concat(parts, ignore_index=True)
    d["loop"] = d.sequence_id.str.removeprefix("LEO-")
    return d.drop_duplicates("loop", keep="last")


# --------------------------------------------------------------------------- stages
def import_abpipe() -> list[str]:
    """10-mer windows that Jahn's abpipe loop selected, read-only. An independent implementation."""
    loops = set()
    for f in sorted(JAHN_RUNS.glob("*/round_*/selected.csv")):
        try:
            d = pd.read_csv(f, usecols=["info__cdrh3"])
        except Exception:
            continue
        for c in d["info__cdrh3"].dropna():
            w = c[2:12]
            if len(c) == 13 and c.startswith("SR") and c.endswith("Y") and set(w) <= set(dl.AA):
                loops.add(w)
    return sorted(loops)


def stage_abguard(pool, loops, rdir, ctx):
    if not loops:
        return pool
    t0 = time.time()
    cached = abguard_cache()
    have = cached[cached["loop"].isin(loops)] if len(cached) else cached
    todo = [l for l in loops if l not in set(have["loop"])] if len(have) else list(loops)
    parts = [have] if len(have) else []
    if len(have):
        log(f"abguard: reusing {len(have)} earlier real verdicts, validating {len(todo)} new")
    if todo:
        parts.append(vab.validate(todo, "optimisation", shards=14, keep_dir=rdir / f"abguard_{int(time.time())}"))
    rep = pd.concat(parts, ignore_index=True)
    rep = rep.rename(columns={"passes": "abg_pass", "track_assigned": "abg_track", "cdrh3": "abg_cdrh3",
                              "cdrh3_min_distance": "abg_min_dist",
                              "cdrh3_identity_to_trastuzumab": "abg_ident", "failures": "abg_failures"})
    rep["cdrh3"] = rep["loop"]
    off = rep.abg_track != "novel"        # abguard handed it to the other track: rejected, never reclassified
    rep.loc[off, "abg_pass"] = False
    fast = ctx.novelty.min_distance_capped(encode(rep.cdrh3.tolist()))
    theirs = np.minimum(rep.abg_min_dist.astype(int).to_numpy(), 3)
    log(f"abguard (real, optimisation track): {len(rep)} designs in {time.time() - t0:.0f}s, "
        f"{int(rep.abg_pass.sum())} PASS, {int(off.sum())} assigned off-track; "
        f"fast novelty index agrees on {int((fast == theirs).sum())}/{len(rep)}")
    if (fast != theirs).any():
        bad = rep[fast != theirs][["cdrh3", "abg_min_dist"]].assign(fast=fast[fast != theirs])
        bad.to_csv(OUT / "index_disagreements.csv", mode="a", header=not (OUT / "index_disagreements.csv").exists(),
                   index=False)
        log("  DISAGREEMENT, first rows:\n" + bad.head(5).to_string(index=False))
    return put(pool, rep, ABG_COLS)


def nm_cache() -> dict:
    """Tm from earlier NanoMelt runs, keyed by CDRH3 window, so a restart never recomputes them."""
    out = {}
    for fa in sorted(gs.GPU.glob("nm_*/in.fa")):
        csvf = fa.parent / "nanomelt.csv"
        if not csvf.exists():
            continue
        vh, name = {}, None
        for line in fa.read_text().splitlines():
            if line.startswith(">"):
                name = line[1:].strip()
            elif name:
                vh[name] = line.strip()
        d = pd.read_csv(csvf)
        tm = dict(zip(d["ID"].astype(str), d["NanoMelt Tm (C)"].astype(float)))
        if "WT" not in tm:
            continue
        for k, seq in vh.items():
            if k != "WT" and k in tm:
                out[seq[gs.LOOP_START: gs.LOOP_START + 10]] = (tm[k], tm[k] - tm["WT"])
    return out


def stage_gpu(pool, rnd):
    """NanoMelt only: seconds, and cached. Boltz-2 runs in gpu_worker.py whenever the card is free."""
    ok = pool.abg_pass == True   # noqa: E712
    need = pool[ok & pool.tm_pred.isna()].cdrh3.tolist()
    if not need:
        return pool
    cache = nm_cache()
    got = [l for l in need if l in cache]
    if got:
        df = pd.DataFrame({"cdrh3": got, "tm_pred": [cache[l][0] for l in got],
                           "dtm_vs_wt": [cache[l][1] for l in got]})
        pool = put(pool, df, NM_COLS)
        log(f"NanoMelt: reused {len(got)} earlier predictions")
    need = [l for l in need if l not in cache]
    if need:
        try:
            pool = put(pool, gs.run_nanomelt(need, f"nm_r{rnd:02d}_{int(time.time())}"), NM_COLS)
        except Exception as e:
            log(f"NanoMelt failed: {e}")
    return pool


# --------------------------------------------------------------------------- one round
def generate(rnd, cfg, ctx, pool, boot_seed):
    """One generation run under the ledger protocol: check-run, add-run, run, set-status."""
    rdir = OUT / f"round_{rnd:02d}"
    refine_pool, parents_hash = None, ""
    if cfg.get("refine_top") and not pool.empty:
        refine_pool = pool[pool.abg_pass == True]   # noqa: E712
        top = refine_pool.sort_values("R").head(cfg["refine_top"]).cdrh3.tolist()
        parents_hash = hashlib.sha256(",".join(sorted(top)).encode()).hexdigest()[:10]
    flags = run_flags(cfg, boot_seed, parents_hash)
    rid = f"leo-dl-r{rnd:02d}-{hashlib.sha256(' '.join(flags).encode()).hexdigest()[:6]}"
    code, _ = led("check-run", "--ledger-dir", str(LED), *flags)
    log(f"ledger check-run: {'SEEN (exit 3)' if code == 3 else 'NEW'}  run_id {rid}")
    notes = f"round {rnd} '{cfg['name']}': enumerate arms, gate, bootstrap-score, diverse shortlist"
    if code == 3:
        if (rdir / "shortlist.csv").exists():
            log("this exact run is already in the ledger and its output exists: reusing it, not regenerating")
            return pd.read_csv(rdir / "shortlist.csv"), rid
        flags += ["--setting", "rerun=2"]
        rid += "-rerun2"
        notes += " | deliberate repeat: the logged run's output folder is missing"
    code, txt = led("add-run", "--ledger-dir", str(LED), *flags, "--run-id", rid,
                    "--output", str(rdir), "--notes", notes)
    if code != 0:
        raise RuntimeError(f"ledger add-run failed: {txt}")
    try:
        sl = dl.run_round(rnd, cfg, ctx, refine_pool)
    except Exception:
        led("set-status", "--ledger-dir", str(LED), "--run-id", rid, "--status", "failed",
            "--notes", "exception, see loop.log")
        raise
    led("set-status", "--ledger-dir", str(LED), "--run-id", rid, "--status", "done" if len(sl) else "empty")
    return sl, rid


def do_round(rnd, cfg, ctx, pool, boot_seed, fold_n):
    t0 = time.time()
    out = generate(rnd, cfg, ctx, pool, boot_seed)
    sl, rid = out
    if sl is None or len(sl) == 0:
        log("round produced no shortlist")
        return pool, 0
    rdir = OUT / f"round_{rnd:02d}"
    pool, n_new = merge_shortlist(pool, sl, rnd)
    mine = pool[pool.cdrh3.isin(sl.cdrh3)]
    # the rule: check the ledger before validating, refolding or ranking any table of designs
    pre = rdir / "base_validation_pre.csv"
    validation_table(mine).to_csv(pre, index=False)
    _, txt = led("check", str(pre), "--ledger-dir", str(LED), "--new-only", str(rdir / "fresh.csv"))
    log("ledger check: " + next((l for l in txt.splitlines() if " designs:" in l), "?"))
    todo = mine[mine.abg_pass.isna()].cdrh3.tolist()
    pool = stage_abguard(pool, todo, rdir, ctx)
    post = rdir / "base_validation.csv"
    validation_table(pool[pool.cdrh3.isin(sl.cdrh3)]).to_csv(post, index=False)
    code, txt = led("add", str(post), "--ledger-dir", str(LED), "--run-id", rid)
    log("ledger add: " + (txt.splitlines()[-1] if txt else "?") + ("" if code == 0 else f"  (exit {code})"))
    pool = merge_gpu(stage_gpu(pool, rnd))
    save_pool(pool)
    log(f"round {rnd} done in {time.time() - t0:.0f}s: shortlist {len(sl)}, new to the pool {n_new}, "
        f"pool {len(pool)}, abguard PASS in pool {int((pool.abg_pass == True).sum())}")   # noqa: E712
    return pool, n_new


# --------------------------------------------------------------------------- portfolio
def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def build_portfolio(pool: pd.DataFrame, n: int = 20, min_div: int = 3) -> pd.DataFrame:
    """A ranked list of up to n that hedges across hypotheses, since only the top 4 are guaranteed a test.

    The four lead slots are four different bets, not four copies of the best score: the
    surrogate has one optimum, so the first four by score would be near-duplicates and
    would all fail together. After them the list fills by lower-bound predicted KD.
    """
    if pool.empty:
        return pd.DataFrame()
    el = pool[(pool.abg_pass == True) & (pool.support_min >= dl.MIN_SUPPORT)        # noqa: E712
              & (pool.p_binder_p10 >= dl.P_THR)].copy()
    if el.empty:
        return el
    el["gpu_ok"] = el.iptm.notna() & (el.pose_flag == False)                        # noqa: E712
    el["tm_ok"] = el.dtm_vs_wt.notna() & (el.dtm_vs_wt.astype(float) > -4.1)
    el["in_abpipe"] = el.arm.astype(str).str.contains("abpipe")
    el = el.sort_values("R")
    slots = [
        ("affinity push: best lower-bound predicted KD", el),
        ("epitope conservative: contact risk <= 0.75", el[el.epi_risk <= 0.75]),
        ("independent consensus: also chosen by Jahn abpipe", el[el.in_abpipe & (el.epi_risk <= 1.5)]),
        ("structure and Tm clean: folds in pose, no Tm penalty", el[el.gpu_ok & el.tm_ok]),
    ]
    chosen: list[str] = []
    label: list[str] = []

    def room(l: str) -> bool:
        return all(hamming(l, c) >= min_div for c in chosen)

    for name, cand in slots:
        for _, r in cand.iterrows():
            if r.cdrh3 not in chosen and room(r.cdrh3):
                chosen.append(r.cdrh3)
                label.append(name)
                break
    for _, r in el.iterrows():
        if len(chosen) >= n:
            break
        if r.cdrh3 not in chosen and room(r.cdrh3):
            chosen.append(r.cdrh3)
            label.append("fill: next best lower bound")
    out = el.set_index("cdrh3").loc[chosen].reset_index()
    out.insert(0, "slot", label)
    out.insert(0, "rank", range(1, len(out) + 1))
    out["heavy_sequence"] = [vh_with_loop(l) for l in out.cdrh3]
    out["light_sequence"] = TZ_VL
    return out


def stability(el: pd.DataFrame, ctx, seeds=(101, 102, 103), top: int = 60) -> dict:
    """How often a design stays in the top `top` when the training data are resampled afresh."""
    codes = encode(el.cdrh3.tolist()).astype(np.int16)
    hits = pd.Series(0.0, index=el.cdrh3.values)
    for s in seeds:
        e = dl.ensemble_predict(dl.bootstrap_tables(ctx.spr, B=100, seed=s), codes)
        R = e["kd_mean"] + dl.K_LCB * e["kd_sd"]
        hits.loc[el.cdrh3.values[np.argsort(R)[:top]]] += 1
    return (hits / len(seeds)).to_dict()


def write_md(port: pd.DataFrame, ctx, extra: dict) -> None:
    tz = dl.ensemble_predict(ctx.boot, dl.TZ_CODES[None, :].astype(np.int16))
    stamp = time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime())
    L = ["# Portfolio from the design loop", "",
         f"Written {stamp}. Every row passed the real abguard optimisation checker (BASE VALIDATION) "
         "and is assigned `novel`, which is the optimisation track.", "",
         "**What the numbers are.** `pred` is the ranker bootstrap mean of log10 KD (nM) plus one standard "
         "deviation, so lower is better and uncertain designs are penalised. It is a ranking, not a KD: the "
         f"same ranker puts trastuzumab at {float(tz['kd_mean'][0]):+.2f} +/- {float(tz['kd_sd'][0]):.2f}, "
         "against a measured log10(0.152) = -0.82. The surrogate is additive, validated only near measured "
         "variants, and the top of a search over tens of millions of options is where its error is largest. "
         "It says nothing about thermostability or epitope preservation.", "",
         "| # | slot | CDRH3 | substitutions | pred | p10 binder | epi risk | nearest measured (d, KD nM) "
         "| dTm C | H3 A | P(beats tz) | top-60 freq |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in port.iterrows():
        dtm = "" if pd.isna(r.dtm_vs_wt) else f"{float(r.dtm_vs_wt):+.1f}"
        h3 = "" if pd.isna(r.h3_rmsd_pose) else f"{float(r.h3_rmsd_pose):.1f}"
        if r.pose_flag == True:   # noqa: E712
            h3 += "!"
        L.append(f"| {r['rank']} | {r.slot.split(':')[0]} | `{r.cdrh3}` | {r.subs} | {r.R:+.2f} | "
                 f"{r.p_binder_p10:.2f} | {r.epi_risk:.2f} | {r.nearest_binder} ({int(r.d_meas)}, "
                 f"{r.nearest_binder_kd_nM:.2g}) | {dtm} | {h3} | {r.p_beats_tz:.0%} | {extra.get(r.cdrh3, float('nan')):.0%} |")
    (OUT / "PORTFOLIO.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def finalize(pool: pd.DataFrame, ctx, fold_all: bool = True, gpu_wait_s: int = 900) -> pd.DataFrame:
    """Every portfolio member gets a Tm; folds are collected from the worker; both tracks are checked."""
    pool = merge_gpu(stage_gpu(pool, 99))
    port = build_portfolio(pool)
    if port.empty:
        log("finalize: empty portfolio")
        return pool
    t0 = time.time()
    while fold_all and time.time() - t0 < gpu_wait_s:
        missing = port[port.iptm.isna()].cdrh3.tolist()
        if not missing:
            break
        log(f"finalize: waiting for the GPU worker to fold {len(missing)} of {len(port)} portfolio designs "
            f"({int(time.time() - t0)}s of {gpu_wait_s}s)")
        time.sleep(60)
        pool = merge_gpu(pool)
        port = build_portfolio(pool)
    save_pool(pool)
    log("finalize: BASE VALIDATION, both track scripts (the other track is expected to FAIL)")
    loops = port.cdrh3.tolist()
    opt = vab.validate(loops, "optimisation", shards=4, keep_dir=OUT / "final_optimisation")
    den = vab.validate(loops, "denovo", shards=4, keep_dir=OUT / "final_denovo")
    bv = opt[["loop", "track_assigned", "passes"]].rename(columns={"passes": "optimisation_pass"}).merge(
        den[["loop", "passes"]].rename(columns={"passes": "denovo_pass"}), on="loop")
    bv.to_csv(OUT / "portfolio_base_validation.csv", index=False)
    log(f"  optimisation PASS {int(bv.optimisation_pass.sum())}/{len(bv)}; "
        f"de novo PASS {int(bv.denovo_pass.sum())}/{len(bv)}")
    el = pool[(pool.abg_pass == True) & (pool.support_min >= dl.MIN_SUPPORT)].copy()   # noqa: E712
    st = stability(el, ctx)
    port["stability"] = port.cdrh3.map(st)
    port.to_csv(OUT / "portfolio.csv", index=False)
    write_md(port, ctx, st)
    log(f"finalize: portfolio of {len(port)} written to out/portfolio.csv and out/PORTFOLIO.md")
    return pool


# --------------------------------------------------------------------------- schedule and main
SCHEDULE = [
    dict(name="broad", product_arms=["conservative", "balanced", "affinity"], m_balanced=9, m_affinity=7,
         keep=800_000, shell_seeds=120, shell_radii=[3], m_shell=7, keep_per_seed=6000, per_arm=150,
         min_div=2, import_abpipe=True),
    dict(name="wide", product_arms=["balanced", "affinity"], m_balanced=10, m_affinity=8, keep=1_000_000,
         shell_seeds=187, shell_radii=[3, 4], m_shell=7, keep_per_seed=5000, per_arm=150, min_div=2),
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--until", default="", help="HH:MM UTC: do not start a round after this")
    ap.add_argument("--max-rounds", type=int, default=10)
    ap.add_argument("--fold-per-round", type=int, default=24)
    ap.add_argument("--boot-seed", type=int, default=0)
    ap.add_argument("--gpu-wait", type=int, default=15, help="minutes finalize waits for the GPU worker")
    ap.add_argument("--no-boltz", action="store_true", help="skip Boltz-2 folds (test runs)")
    a = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    ctx = dl.Context(B=100, seed=a.boot_seed)
    pool = load_pool()
    abp = import_abpipe()
    log(f"imported {len(abp)} 10-mer windows chosen by Jahn abpipe (read-only) as an independent arm; "
        f"pool on disk: {len(pool)}")
    for rnd in range(1, a.max_rounds + 1):
        if (OUT / "STOP").exists():
            log("STOP file found, finishing")
            break
        if a.until and time.strftime("%H:%M", time.gmtime()) >= a.until:
            log(f"time budget reached ({a.until} UTC), finishing")
            break
        if rnd <= len(SCHEDULE):
            cfg = dict(SCHEDULE[rnd - 1])
        else:
            cfg = dict(name=f"refine-{rnd}", refine_top=60, per_arm=120, min_div=1)
        if cfg.pop("import_abpipe", False) and abp:
            cfg["extra_arms"] = {"abpipe": abp}
        try:
            pool, n_new = do_round(rnd, cfg, ctx, pool, a.boot_seed, 0 if a.no_boltz else a.fold_per_round)
        except Exception:
            log("round failed:\n" + traceback.format_exc())
            continue
        port = build_portfolio(pool)
        port.to_csv(OUT / "portfolio.csv", index=False)
        log(f"portfolio after round {rnd}: {len(port)} designs; lead {port.cdrh3.iloc[0] if len(port) else '-'}")
        if rnd > len(SCHEDULE) and n_new == 0:
            log("a refine round found nothing new: the search has converged")
            break
    pool = finalize(pool, ctx, fold_all=not a.no_boltz, gpu_wait_s=a.gpu_wait * 60)
    log("loop finished")


if __name__ == "__main__":
    main()
