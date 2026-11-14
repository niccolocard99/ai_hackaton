#!/usr/bin/env python3
"""Generate trastuzumab CDRH3-window variants -> hard rules -> NICCOLO/eda prefilter -> diverse top list.

Optimisation track only (10-residue CDRH3 window on the trastuzumab VH/VL framework). Run with the abmodels venv:
    /home/hackathon/.venvs/abmodels/bin/python gen_score.py

Stages (counts are written to out/pool_stats.json):
  1. sample N_RAW variants with k in {4..7} substitutions; residues drawn from a position-specific enrichment
     table built from the FACS 'high' class (a proposal distribution only; the classifier is what scores).
  2. hard rules: k>=4 (IMGT CDRH3 identity < 70% needs >=4 edits on the 13-mer), no Cys (unpaired-Cys fail),
     no new N-glycosylation sequon, standard residues, length exactly 10.
  3. prefilter.score(): gate p_binder >= GATE_P, rank by pred_log10_kd ascending (a RANKING, not a KD).
  4. 'hole' check: >=3 Hamming from EVERY provided window (all 524,346 FACS + SPR windows).
  5. greedy diversification: selected windows differ pairwise in >= MIN_PAIR_HAMMING positions.
The RFD4 her2-novel seed-11 windows are pushed through the same stages as a labelled comparison source.
"""
import itertools, json, multiprocessing as mp, os, re, sys, time, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")

OPT = "/home/hackathon/nick-work/optim"
os.chdir(f"{OPT}/eda"); sys.path.insert(0, ".")
from prefilter import PreFilter, load, AA, TZ          # noqa: E402  (user's classifier, copied verbatim)

H = "/home/hackathon"
OUT = f"{OPT}/out"
N_RAW, K_CHOICES, K_P = 4_000_000, [4, 5, 6, 7], [.35, .35, .20, .10]
GATE_P = 0.80          # p_binder gate (the README sets none; recorded here and in pool_stats.json)
N_HOLE = 40_000        # best-by-predicted-KD candidates sent to the (expensive) hole check
N_SELECT = 200         # diversified list handed to the official checker
MIN_PAIR_HAMMING = 3
SEED, NPROC = 7, 8
T0 = time.time()
log = lambda *a: print(f"[{time.time() - T0:6.1f}s]", *a, flush=True)

IDXA = np.full(256, 255, np.uint8)
for i, c in enumerate(AA): IDXA[ord(c)] = i
W = 20 ** np.arange(9, -1, -1, dtype=np.int64)
A = {c: i for i, c in enumerate(AA)}


def enc(seqs):
    seqs = list(seqs)
    assert all(len(s) == 10 for s in seqs), "every window must be exactly 10 residues (score() does NOT check per sequence)"
    o = IDXA[np.frombuffer("".join(seqs).encode(), np.uint8).reshape(-1, 10)]
    assert (o != 255).all(), "non-standard residue"
    return o


dec = lambda S: np.array(list(AA), dtype="S1")[S].view("S10").ravel().astype(str)
WT = enc([TZ])[0]


def hard_rules(S):
    """True where a window passes: no Cys, no new N-X-S/T sequon (checked with the VH flanks SR..YWG)."""
    n = len(S)
    ext = np.concatenate([np.tile(np.array([A["S"], A["R"]], np.uint8), (n, 1)), S,
                          np.tile(np.array([A["Y"], A["W"], A["G"]], np.uint8), (n, 1))], 1)
    sequon = np.zeros(n, bool)
    for i in range(ext.shape[1] - 2):
        sequon |= (ext[:, i] == A["N"]) & (ext[:, i + 1] != A["P"]) & ((ext[:, i + 2] == A["S"]) | (ext[:, i + 2] == A["T"]))
    return (S != A["C"]).all(1) & ~sequon


# ---------------------------------------------------------------- data
spr, facs = load()
t = time.time(); pf = PreFilter().fit(spr, facs); log(f"classifier fitted in {time.time() - t:.1f}s on {len(spr)} SPR / {len(facs)} FACS rows")
stats = {"classifier_fit_rows": {"spr": len(spr), "facs": len(facs)}, "gate_p_binder": GATE_P}

spr_all = pd.read_csv("data/trastuzumab_her2_spr_all_768.csv")
facs_all = pd.concat([pd.read_csv(f"data/splits/{s}.csv") for s in ("train", "val", "test")], ignore_index=True)
provided = sorted(set(facs_all.seq) | set(spr_all.cdrh3[spr_all.cdrh3.str.len() == 10]))
Renc = enc(provided); Rcodes = np.unique(Renc.astype(np.int64) @ W)
log(f"provided windows (FACS train+val+test + SPR, 10-mers): {len(Rcodes):,}")
stats["provided_windows"] = int(len(Rcodes))

# ---------------------------------------------------------------- proposal distribution (FACS 'high' enrichment)
fe = enc(facs.seq.values); y = facs.label.values.astype(int); base = y.mean()
E = np.zeros((10, 20))
for p in range(10):
    a_all = np.bincount(fe[:, p], minlength=20).astype(float); a_hi = np.bincount(fe[y == 1, p], minlength=20).astype(float)
    E[p] = np.log(((a_hi + 2 * base) / (a_all + 2)) / base)
q = np.exp(E); q[np.arange(10), WT] = 0; q /= q.sum(1, keepdims=True)
cdf = np.cumsum(q, 1); cdf[:, -1] = 1.0

# ---------------------------------------------------------------- 1. sample
rng = np.random.default_rng(SEED)
ks = rng.choice(K_CHOICES, size=N_RAW, p=K_P)
rank = rng.random((N_RAW, 10)).argsort(1).argsort(1)
mut = rank < ks[:, None]; del rank
S = np.tile(WT, (N_RAW, 1))
for p in range(10):
    idx = np.flatnonzero(mut[:, p])
    S[idx, p] = np.searchsorted(cdf[p], rng.random(idx.size)).clip(0, 19)
del mut
codes = S.astype(np.int64) @ W
_, first = np.unique(codes, return_index=True); S = S[first]
stats["sampled"] = N_RAW; stats["unique_after_dedupe"] = int(len(S)); log(f"sampled {N_RAW:,} -> {len(S):,} unique")
nm = (S != WT).sum(1); assert nm.min() >= 4

# ---------------------------------------------------------------- 2. hard rules
ok = hard_rules(S); S = S[ok]; nm = nm[ok]
stats["after_hard_rules"] = int(len(S)); log(f"after hard rules (k>=4, no Cys, no new sequon): {len(S):,}")

# ---------------------------------------------------------------- 3. classifier
def score_all(Sx, chunk=400_000):
    parts = []
    for i in range(0, len(Sx), chunk):
        parts.append(pf.score(dec(Sx[i:i + chunk])))
    return pd.concat(parts, ignore_index=True)

t = time.time(); sc = score_all(S); log(f"scored {len(sc):,} windows in {time.time() - t:.1f}s")
sc["src"] = "insilico"
stats["p_binder_quantiles_all_k>=4"] = {str(f): round(float(np.quantile(sc.p_binder, f)), 3) for f in (.1, .5, .9, .99)}

# comparison source: RFD4 her2-novel seed-11 designs (CDRH3 B99-108 on the trastuzumab framework)
m = f"{H}/hackathon-team-3/apheris/rfd4/runs/pipe-her2novel-s11/design/generation/inference_metrics/per_example_metrics.csv"
rw = sorted({mt.group(1) for r in pd.read_csv(m).designed_sequences if (mt := re.search(r"YYCSR(.{10})YWGQG", r))})
Sr = enc(rw); keep_r = hard_rules(Sr) & ((Sr != WT).sum(1) >= 4)
sr = score_all(Sr[keep_r]); sr["src"] = "rfd4_her2novel"
stats["rfd4_her2novel"] = {"unique_windows": len(rw), "pass_hard_rules_k>=4": int(keep_r.sum()),
                           "p_binder>=gate": int((sr.p_binder >= GATE_P).sum())}
log("rfd4 her2-novel comparison source:", stats["rfd4_her2novel"])

pool = pd.concat([sc, sr], ignore_index=True)
pool["n_mut"] = pool.n_mut.astype(int)
g = pool[pool.p_binder >= GATE_P].sort_values("pred_log10_kd").reset_index(drop=True)
stats["pass_gate"] = int(len(g)); log(f"pass p_binder >= {GATE_P}: {len(g):,}")
top = g.head(N_HOLE).copy()
Stop = enc(top.cdrh3.values)

# ---------------------------------------------------------------- 4. hole check (Hamming <= 2 to any provided window)
def member(c):
    i = np.searchsorted(Rcodes, c.ravel()); i[i == len(Rcodes)] = 0
    return (Rcodes[i] == c.ravel()).reshape(c.shape).any(1)

def hole_flags(Sb):
    m_ = len(Sb); base_ = Sb.astype(np.int64) @ W
    d1 = (np.arange(20, dtype=np.int64)[None, None, :] - Sb[:, :, None].astype(np.int64)) * W[None, :, None]
    hit = member(base_[:, None] + d1.reshape(m_, -1))
    for p, q_ in itertools.combinations(range(10), 2):
        hit |= member(base_[:, None] + (d1[:, p, :, None] + d1[:, q_, None, :]).reshape(m_, -1))
    return hit

t = time.time()
batches = [Stop[i:i + 250] for i in range(0, len(Stop), 250)]
with mp.get_context("fork").Pool(NPROC) as pool_:
    hole = np.concatenate(pool_.map(hole_flags, batches))
top["within2_of_provided"] = hole
stats["hole_checked"] = int(len(top)); stats["fail_hole_(<=2_from_a_provided_window)"] = int(hole.sum())
log(f"hole check on {len(top):,} best candidates in {time.time() - t:.1f}s: {int(hole.sum()):,} within 2 of a provided window")
surv = top[~hole].reset_index(drop=True); Ssurv = enc(surv.cdrh3.values)
stats["survive_all_rules_and_gate"] = int(len(surv))

# ---------------------------------------------------------------- 5. diversify
from collections import Counter
MAX_SUB_SHARE = 0.25     # no single substitution (e.g. A106E) in more than 25% of the list at any prefix size
surv_subs = [{f"{TZ[p]}{99 + p}{w[p]}" for p in range(10) if w[p] != TZ[p]} for w in surv.cdrh3]
cnt, chosen = Counter(), []
for i in range(len(surv)):
    if chosen and (Ssurv[chosen] != Ssurv[i]).sum(1).min() < MIN_PAIR_HAMMING: continue
    cap = max(2, int(np.ceil(MAX_SUB_SHARE * (len(chosen) + 1))))
    if any(cnt[m] + 1 > cap for m in surv_subs[i]): continue
    chosen.append(i); cnt.update(surv_subs[i])
    if len(chosen) == N_SELECT: break
stats["max_substitution_share"] = MAX_SUB_SHARE
sel = surv.iloc[chosen].copy().reset_index(drop=True)
Ssel = Ssurv[chosen]
sel["rank"] = np.arange(1, len(sel) + 1)
sel["mutations"] = [",".join(f"{TZ[p]}{99 + p}{w[p]}" for p in range(10) if w[p] != TZ[p]) for w in sel.cdrh3]
sel["min_hamming_to_provided"] = [int((Renc != s).sum(1).min()) for s in Ssel]
sel = sel[["rank", "cdrh3", "n_mut", "mutations", "p_binder", "pred_log10_kd", "p_facs_high", "min_hamming_to_provided", "src"]]
sel = sel.rename(columns={"cdrh3": "window"})
sel.to_csv(f"{OUT}/selected_diverse_{N_SELECT}.csv", index=False)
g.head(5000).to_csv(f"{OUT}/pool_top5000_by_pred_kd_before_hole_check.csv", index=False)
stats["selected"] = int(len(sel)); stats["selected_by_source"] = sel.src.value_counts().to_dict()
stats["selected_n_mut"] = {int(k): int(v) for k, v in sel.n_mut.value_counts().sort_index().items()}
stats["selected_min_hamming_to_provided"] = {int(k): int(v) for k, v in sel.min_hamming_to_provided.value_counts().sort_index().items()}
json.dump(stats, open(f"{OUT}/pool_stats.json", "w"), indent=2)
log("done"); print(json.dumps(stats, indent=1))
print(sel.head(12).to_string(index=False))
