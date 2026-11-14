#!/usr/bin/env python3
"""One generation round for the HER2 optimisation track.

    enumerate (several arms) -> hard gates -> score with uncertainty -> diverse selection

Everything here is CPU and numpy. The GPU stages (NanoMelt, Boltz-2) and the real
abguard run live in their own modules and consume the shortlist this writes.

The three gates, all measured (probe_abguard.py, feasible.py), none assumed:

  identity   >= 4 substitutions from trastuzumab: the 13-mer CDRH3 SRWGGDGFYAMDY
             must fall under 70% identity, and 4/13 -> 69.2% is the first rung
             that does.
  novelty    >= 3 mutations from every one of the 525,063 reference CDRH3s. This
             is the gate that bites: every measured SPR variant is itself in the
             reference tables, so none of them can be submitted.
  scope      the nearest measured binder is <= D_MAX away. The ranker is an
             interpolator: its skill falls from Spearman +0.67 to +0.43 as the
             nearest training sequence goes from 3 to 5 mutations away, and to
             -0.06 at 6 (NICCOLO/eda/AGENT_GUIDE.md section 3). Novelty already
             forces the distance to be >= 3, so the usable window is 3 to D_MAX.

Novelty is applied before ranking, as root AGENTS.md requires: the best-scoring
designs are exactly the ones the novelty filter rejects.

Selection uses a lower confidence bound over a cluster bootstrap, not the point
estimate, because the candidates are the argmax of a surrogate with validated
Spearman 0.63 taken over tens of millions of options (see models.py).
"""
from __future__ import annotations

import json
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, "/home/hackathon/wt/leo-design-loop/NICCOLO/eda")

from her2lib import CONTACT_FREQ, N_STRUCT, POSITIONS, TZ_LOOP, epitope_risk, subs  # noqa: E402
from models import TZ_CODES, bootstrap_tables, ensemble_predict, fit_point  # noqa: E402
from novelty import AA, NoveltyIndex, encode, load_reference, min_dist_to  # noqa: E402
from prefilter import load as pf_load  # noqa: E402

DATA, REF, OUT = HERE / "data", HERE / "ref", HERE / "out"
C_CODE = AA.index("C")
_IDX = np.arange(10)[None, :]
CONTACT = np.array([CONTACT_FREQ[p] for p in POSITIONS])  # per index 0..9

MIN_SUB = 4          # identity gate, see module docstring
D_MAX = 4            # scope gate: nearest measured binder within this many mutations
P_THR = 0.80         # classifier gate on the bootstrap 10th percentile of P(binder)
MIN_SUPPORT = 2      # every substituted residue must be carried by >= this many measured binders
K_LCB = 1.0          # selection score = mean + K_LCB * sd of predicted log10 KD


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def decode(codes: np.ndarray) -> list[str]:
    lut = np.frombuffer(AA.encode(), np.uint8)
    return [b.tobytes().decode() for b in lut[np.asarray(codes)]]


# --------------------------------------------------------------------------- context
class Context:
    """Everything that is expensive to build and shared by every round."""

    def __init__(self, B: int = 100, seed: int = 0):
        t0 = time.time()
        ref_codes, self.ref_counts = load_reference(DATA, REF)
        self.novelty = NoveltyIndex(ref_codes)
        log(f"novelty index: {self.ref_counts['unique_loops']:,} unique reference loops "
            f"({time.time() - t0:.0f}s)")

        self.spr, _ = pf_load()
        self.tab, self.pf = fit_point(self.spr)
        self.boot = bootstrap_tables(self.spr, B=B, seed=seed)
        log(f"bootstrap: {len(self.boot)} cluster-resampled refits")

        q = self.spr[self.spr.binding_result == "quantified_binding"].reset_index(drop=True)
        self.b_loops = q.cdrh3.tolist()
        self.b_codes = encode(self.b_loops).astype(np.int16)
        self.b_kd = q.kd_nM.to_numpy(float)
        self.n_binders = len(q)

        # evidence behind every (position, residue): how many measured variants carry it
        allc = encode(self.spr.cdrh3.tolist()).astype(np.int16)
        self.support = np.zeros((10, 20), int)
        self.support_b = np.zeros((10, 20), int)
        for p in range(10):
            np.add.at(self.support[p], allc[:, p], 1)
            np.add.at(self.support_b[p], self.b_codes[:, p], 1)

        # additive goodness table, higher = better, used only to enumerate. Selection
        # never uses it: it is replaced by the bootstrap lower bound below.
        pred_b = self.tab.pred_log10_kd(self.b_codes)
        sig_r = float(pred_b.std())
        sig_c = float(self.tab.logit(allc).std())
        self.G = -self.tab.Wr / sig_r + 0.5 * self.tab.Wc / sig_c
        self.G[:, C_CODE] = -1e9  # no cysteine: an odd one is unpaired and abguard rejects it
        self.tz_pred = float(self.tab.pred_log10_kd(TZ_CODES[None, :].astype(np.int16))[0])

    def top_alphabet(self, p: int, m: int, min_support: int = 2) -> np.ndarray:
        order = np.argsort(-self.G[p])
        keep = [a for a in order if a != C_CODE and self.support[p, a] >= min_support][:m]
        if TZ_CODES[p] not in keep:
            keep.append(int(TZ_CODES[p]))
        return np.array(sorted(set(int(a) for a in keep)), np.int16)

    def observed_alphabet(self, p: int, min_support: int = 3) -> np.ndarray:
        ks = [a for a in range(20) if a != C_CODE and self.support[p, a] >= min_support]
        if TZ_CODES[p] not in ks:
            ks.append(int(TZ_CODES[p]))
        return np.array(sorted(set(ks)), np.int16)


# --------------------------------------------------------------------------- arms
def product_topk(P, A, G, wt, min_sub, keep, block_elems=4_000_000):
    """Top-`keep` of the full Cartesian product over positions P, by additive score G.

    Positions outside P stay at trastuzumab. The product is split into two halves
    and scored by broadcasting, so a 24M-way space never exists in memory at once.
    Returns (codes, score, n_enumerated).
    """
    half = (len(P) + 1) // 2
    Lp, Rp = list(P[:half]), list(P[half:])

    def build(pos):
        if not pos:
            return np.zeros((1, 0), np.int16), np.zeros(1), np.zeros(1, np.int16)
        grids = np.meshgrid(*[A[p] for p in pos], indexing="ij")
        codes = np.stack([g.reshape(-1) for g in grids], 1).astype(np.int16)
        score = np.zeros(len(codes))
        nsub = np.zeros(len(codes), np.int16)
        for i, p in enumerate(pos):
            score += G[p, codes[:, i]]
            nsub += (codes[:, i] != wt[p]).astype(np.int16)
        return codes, score, nsub

    cl, sl, nl = build(Lp)
    cr, sr, nr = build(Rp)
    n_total = len(cl) * len(cr)
    rows = max(1, block_elems // len(cr))
    S, I, J, held = [], [], [], 0
    for i0 in range(0, len(cl), rows):
        i1 = min(len(cl), i0 + rows)
        tot = sl[i0:i1, None] + sr[None, :]
        tot = np.where((nl[i0:i1, None] + nr[None, :]) >= min_sub, tot, -np.inf)
        flat = tot.ravel()
        k = min(keep, flat.size)
        idx = np.argpartition(-flat, k - 1)[:k]
        idx = idx[np.isfinite(flat[idx])]
        S.append(flat[idx]); I.append(i0 + idx // len(cr)); J.append(idx % len(cr)); held += len(idx)
        if held > 3 * keep:   # running top-K, so a billion-way space never holds more than ~3K rows
            s_, i_, j_ = np.concatenate(S), np.concatenate(I), np.concatenate(J)
            sel = np.argpartition(-s_, keep - 1)[:keep]
            S, I, J, held = [s_[sel]], [i_[sel]], [j_[sel]], keep
    s, ii, jj = np.concatenate(S), np.concatenate(I), np.concatenate(J)
    if len(s) > keep:
        sel = np.argpartition(-s, keep - 1)[:keep]
        s, ii, jj = s[sel], ii[sel], jj[sel]
    codes = np.tile(wt.astype(np.int16), (len(s), 1))
    for i, p in enumerate(Lp):
        codes[:, p] = cl[ii, i]
    for i, p in enumerate(Rp):
        codes[:, p] = cr[jj, i]
    return codes, s, n_total


def shell_topk(seeds, radius, A, G, wt, min_sub, keep_per_seed):
    """All loops at Hamming exactly `radius` from each seed, over alphabet A, top by G.

    Returns (codes, n_enumerated).
    """
    out, n_enum = [], 0
    for s in seeds:
        blocks = []
        for combo in combinations(range(10), radius):
            alts = [A[p][A[p] != s[p]] for p in combo]
            if any(len(a) == 0 for a in alts):
                continue
            grids = np.meshgrid(*alts, indexing="ij")
            n = grids[0].size
            blk = np.tile(s, (n, 1))
            for p, g in zip(combo, grids):
                blk[:, p] = g.ravel()
            blocks.append(blk)
        if not blocks:
            continue
        cand = np.concatenate(blocks)
        n_enum += len(cand)
        cand = cand[(cand != wt).sum(1) >= min_sub]
        if len(cand) > keep_per_seed:
            sc = G[_IDX, cand].sum(1)
            cand = cand[np.argpartition(-sc, keep_per_seed - 1)[:keep_per_seed]]
        out.append(cand)
    return (np.concatenate(out) if out else np.zeros((0, 10), np.int16)), n_enum


def neighbours(codes, A):
    """All single-residue neighbours of each row over alphabet A (the refinement move)."""
    out = []
    for c in codes:
        for p in range(10):
            for a in A[p]:
                if a != c[p]:
                    n = c.copy(); n[p] = a
                    out.append(n)
    return np.array(out, np.int16) if out else np.zeros((0, 10), np.int16)


# --------------------------------------------------------------------------- gates
def gate(codes: np.ndarray, ctx: Context, d_max: int = D_MAX) -> tuple[np.ndarray, dict]:
    """identity -> novelty -> scope. Returns (surviving codes, d_meas for them, funnel)."""
    f = {"in": int(len(codes))}
    if not len(codes):
        return codes, np.zeros(0, np.int8), f
    codes = np.unique(codes, axis=0)
    f["unique"] = int(len(codes))
    codes = codes[(codes != TZ_CODES).sum(1) >= MIN_SUB]
    f["identity_ok"] = int(len(codes))
    codes = codes[(codes != C_CODE).all(1)]
    codes = codes[ctx.novelty.is_novel(codes)]
    f["novel"] = int(len(codes))
    if not len(codes):
        return codes, np.zeros(0, np.int8), f
    d = min_dist_to(codes, ctx.b_codes)
    ok = d <= d_max
    f["in_scope"] = int(ok.sum())
    return codes[ok], d[ok], f


# --------------------------------------------------------------------------- scoring
def score(codes: np.ndarray, d_meas: np.ndarray, ctx: Context) -> pd.DataFrame:
    e = ensemble_predict(ctx.boot, codes, ref=TZ_CODES.astype(np.int16))
    loops = decode(codes)
    # nearest measured binder (identity and KD), for context and for anchoring later
    near_i = np.empty(len(codes), int)
    for i0 in range(0, len(codes), 20000):
        blk = codes[i0 : i0 + 20000]
        dd = (blk[:, None, :] != ctx.b_codes[None, :, :]).sum(2)
        near_i[i0 : i0 + 20000] = dd.argmin(1)
    mutated = codes != TZ_CODES
    sup_b = np.where(mutated, ctx.support_b[_IDX, codes], 10**6).min(1)
    risk = (mutated * CONTACT[None, :]).sum(1) / N_STRUCT
    n_contact = (mutated & (CONTACT[None, :] > 0)).sum(1)
    df = pd.DataFrame({
        "cdrh3": loops,
        "n_sub": mutated.sum(1),
        "subs": ["/".join(subs(l)) for l in loops],
        "kd_mean": e["kd_mean"], "kd_sd": e["kd_sd"],
        "p_binder_p10": e["p_binder_p10"], "p_beats_tz": e["p_beats_ref"],
        "point_log10_kd": ctx.tab.pred_log10_kd(codes),
        "point_p_binder": ctx.tab.p_binder(codes),
        "d_meas": d_meas.astype(int),
        "nearest_binder": [ctx.b_loops[i] for i in near_i],
        "nearest_binder_kd_nM": ctx.b_kd[near_i],
        "epi_risk": risk, "n_contact_pos": n_contact,
        "support_min": sup_b,
    })
    df["R"] = df.kd_mean + K_LCB * df.kd_sd
    return df


def select_diverse(df: pd.DataFrame, n: int, min_div: int, sort_col: str = "R") -> pd.DataFrame:
    """Greedy: best first, skipping anything within `min_div - 1` mutations of a pick."""
    df = df.sort_values(sort_col).reset_index(drop=True)
    codes = encode(df.cdrh3.tolist()).astype(np.int16)
    chosen: list[int] = []
    for i in range(len(df)):
        if len(chosen) >= n:
            break
        if chosen and ((codes[chosen] != codes[i]).sum(1) < min_div).any():
            continue
        chosen.append(i)
    return df.iloc[chosen].reset_index(drop=True)


# --------------------------------------------------------------------------- one round
ARM_DEFS = {
    # name: (mutable index set, how the alphabet is chosen)
    "conservative": dict(positions=[1, 2, 5, 7, 8, 9], note="only positions contacting HER2 in <= 4 of 8 structures"),
    "balanced":     dict(positions=[0, 1, 2, 3, 5, 7, 8, 9], note="drops the two positions contacting in 7 of 8"),
    "affinity":     dict(positions=list(range(10)), note="no contact constraint"),
}


def run_round(rnd: int, cfg: dict, ctx: Context, prior_pool: pd.DataFrame | None) -> pd.DataFrame:
    rdir = OUT / f"round_{rnd:02d}"
    rdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    log(f"=== round {rnd}: {cfg['name']} ===")
    arms: dict[str, tuple[np.ndarray, int]] = {}  # arm -> (codes, n_enumerated)

    # --- product arms: full enumeration over a restricted position set
    for arm in cfg.get("product_arms", []):
        d = ARM_DEFS[arm]
        A = []
        for p in range(10):
            if p not in d["positions"]:
                A.append(np.array([TZ_CODES[p]], np.int16))
            elif arm == "conservative":
                A.append(ctx.observed_alphabet(p, min_support=3))
            else:
                A.append(ctx.top_alphabet(p, cfg["m_" + arm]))
        size = int(np.prod([len(a) for a in A], dtype=np.float64))
        codes, _, n_enum = product_topk(d["positions"], A, ctx.G, TZ_CODES, MIN_SUB, cfg["keep"])
        log(f"  arm {arm:<12} space {size:>13,}  alphabet sizes "
            f"{[len(a) for a in A]}  kept {len(codes):,}")
        arms[arm] = (codes, size)

    # --- seed-shell arm: recombination around the tightest measured binders
    if cfg.get("shell_seeds"):
        order = np.argsort(ctx.b_kd)[: cfg["shell_seeds"]]
        seeds = ctx.b_codes[order]
        for r in cfg["shell_radii"]:
            A = [ctx.top_alphabet(p, cfg["m_shell"]) if r >= 4 else ctx.observed_alphabet(p, 3)
                 for p in range(10)]
            codes, n_enum = shell_topk(seeds, r, A, ctx.G, TZ_CODES, MIN_SUB, cfg["keep_per_seed"])
            arms[f"shell{r}"] = (codes, n_enum)
            log(f"  arm shell{r:<7} seeds {len(seeds)}  enumerated {n_enum:>13,}  kept {len(codes):,}")

    # --- refinement arm: single-residue moves from the pool's best, the loop's feedback
    if cfg.get("refine_top") and prior_pool is not None and len(prior_pool):
        top = prior_pool.sort_values("R").head(cfg["refine_top"])
        A = [ctx.observed_alphabet(p, 3) for p in range(10)]
        codes = neighbours(encode(top.cdrh3.tolist()).astype(np.int16), A)
        arms["refine"] = (codes, len(codes))
        log(f"  arm refine       from top {len(top)} pool designs  moves {len(codes):,}")

    for name, loops in cfg.get("extra_arms", {}).items():
        codes = encode(list(loops)).astype(np.int16)
        arms[name] = (codes, len(codes))
        log(f"  arm {name:<12} imported {len(codes):,} windows")

    # --- gate and score each arm, then select per arm
    funnel, picked = {}, []
    total_enum = 0
    for arm, (codes, n_enum) in arms.items():
        total_enum += n_enum
        g_codes, d, f = gate(codes, ctx)
        f["enumerated"] = int(n_enum)
        funnel[arm] = f
        if not len(g_codes):
            log(f"  {arm:<12} -> 0 survivors  {f}")
            continue
        df = score(g_codes, d, ctx)
        df = df[(df.p_binder_p10 >= P_THR) & (df.support_min >= MIN_SUPPORT)]
        f["classifier_ok"] = int(len(df))
        sel = select_diverse(df, cfg["per_arm"], cfg.get("min_div", 2))
        sel["arm"] = arm
        f["selected"] = int(len(sel))
        picked.append(sel)
        best = sel.iloc[0]
        log(f"  {arm:<12} in {f['in']:>9,} -> novel {f['novel']:>8,} -> scope {f['in_scope']:>8,} "
            f"-> clf {f['classifier_ok']:>7,} -> picked {len(sel):>3}   best {best.cdrh3} "
            f"R {best.R:+.2f} ({best.kd_mean:+.2f}+/-{best.kd_sd:.2f}) epi {best.epi_risk:.2f}")

    if not picked:
        log("  no survivors this round")
        return pd.DataFrame()
    allp = pd.concat(picked, ignore_index=True)
    merged = (allp.sort_values("R").groupby("cdrh3", as_index=False)
              .agg(**{c: (c, "first") for c in allp.columns if c not in ("cdrh3", "arm")},
                   arm=("arm", lambda s: "+".join(sorted(set(s))))))
    merged["round"] = rnd
    merged = merged.sort_values("R").reset_index(drop=True)
    merged.to_csv(rdir / "shortlist.csv", index=False)
    (rdir / "funnel.json").write_text(json.dumps(
        dict(round=rnd, cfg={k: v for k, v in cfg.items() if k != "extra_arms"}, funnel=funnel, total_enumerated=int(total_enum),
             seconds=round(time.time() - t0, 1)), indent=1), encoding="utf-8")
    log(f"  round {rnd}: {total_enum:,} enumerated -> {len(merged)} shortlisted "
        f"in {time.time() - t0:.0f}s")
    return merged


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    ctx = Context()
    log(f"ranker places trastuzumab at log10 KD {ctx.tz_pred:+.3f}; each design below is "
        "judged on its lower bound, not its point estimate")
    cfg = dict(name="smoke", product_arms=["conservative", "balanced", "affinity"],
               m_balanced=6, m_affinity=5, keep=150_000, shell_seeds=30, shell_radii=[3],
               m_shell=6, keep_per_seed=4000, per_arm=60)
    df = run_round(0, cfg, ctx, None)
    print(df.head(15)[["cdrh3", "arm", "n_sub", "subs", "kd_mean", "kd_sd", "R",
                       "p_binder_p10", "d_meas", "epi_risk", "support_min"]].to_string())
