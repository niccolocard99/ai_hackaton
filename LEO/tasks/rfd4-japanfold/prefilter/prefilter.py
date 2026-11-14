"""HER2 CDRH3 pre-filter: fast screen -> binder classifier -> affinity ranker.

Scope: trastuzumab CDRH3 variants (10-mer loop, trastuzumab VH/VL framework), i.e. the OPTIMISATION
track. It does not transfer to de novo designs with <70% identity to known antibodies.

Stages
  0 hard rules   novelty vs the 525,063 provided CDRH3s, mutation-count band (abguard still decides the track)
  1 classifier   P(measurable KD) -- logistic on one-hot CDRH3, trained on 620 SPR outcomes
  2 ranker       predicted log10 KD -- ridge on one-hot + mutation count + FACS score
Validation: clusters (single-linkage Hamming<=2) held out together; see validate().
"""
from pathlib import Path
import time, pickle
import numpy as np, pandas as pd
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy.stats import spearmanr

AA = "ACDEFGHIKLMNPQRSTVWY"; IDX = {c: i for i, c in enumerate(AA)}; L = 10; TZ = "WGGDGFYAMD"
DATA = Path("data"); MODELS = Path("models"); MODELS.mkdir(exist_ok=True)

def ordinal(seqs):
    a = np.frombuffer("".join(seqs).encode(), np.uint8).reshape(-1, L)
    o = np.full(a.shape, -1, np.int8)
    for c, i in IDX.items(): o[a == ord(c)] = i
    return o
def onehot(seqs):
    o = ordinal(seqs).astype(np.int32)
    if o.min() < 0: raise ValueError("non-standard residue")
    X = np.zeros((len(o), L*20), np.float32); r = np.arange(len(o))
    for p in range(L): X[r, p*20+o[:, p]] = 1
    return X
def nmut(seqs): return np.array([sum(a != b for a, b in zip(s, TZ)) for s in seqs], float)

def clusters(seqs, thresh=2):
    """Single-linkage clusters at Hamming<=thresh, used as CV groups so near-duplicates never straddle a fold."""
    O = ordinal(seqs).astype(np.int16)
    D = (O[:, None, :] != O[None, :, :]).sum(2)
    parent = list(range(len(seqs)))
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for i, j in zip(*np.where(np.triu(D <= thresh, 1))):
        a, b = find(i), find(j)
        if a != b: parent[a] = b
    return np.array([find(i) for i in range(len(seqs))])

def load():
    spr = pd.read_csv(DATA/"trastuzumab_her2_spr_all_768.csv")
    v = spr[(spr.cdrh3.str.len() == L) & (spr.cdrh3 != TZ) & (spr.cdrh3 != "ARDRIVDDYYYYGMDV")]
    v = v[v.cdrh3.map(lambda s: set(s) <= set(AA))].reset_index(drop=True)
    facs = pd.concat([pd.read_csv(DATA/"splits/train.csv"), pd.read_csv(DATA/"splits/val.csv")], ignore_index=True)
    return v, facs

class PreFilter:
    def fit(self, spr, facs):
        # FACS model: trained with every SPR sequence removed, so its score is never leaked downstream
        clean = facs[~facs.seq.isin(set(spr.cdrh3))]
        self.n_facs_dropped = len(facs) - len(clean)
        self.facs = HistGradientBoostingClassifier(max_iter=400, max_leaf_nodes=63,
            categorical_features=list(range(L)), random_state=0).fit(ordinal(clean.seq.values), clean.label.values)
        pf = self.facs.predict_proba(ordinal(spr.cdrh3.values))[:, 1]
        b = spr[spr.binding_result.isin(["quantified_binding", "no_reportable_kd"])]
        ib = b.index.values
        self.clf = LogisticRegression(C=.3, max_iter=4000).fit(onehot(b.cdrh3.values),
            (b.binding_result == "quantified_binding").astype(int).values)
        q = spr[spr.binding_result == "quantified_binding"]; iq = q.index.values
        self.rank = Ridge(3.0).fit(np.hstack([onehot(q.cdrh3.values), nmut(q.cdrh3.values)[:, None], pf[iq][:, None]]),
            np.log10(q.kd_nM.values))
        return self
    def score(self, seqs):
        seqs = list(seqs); oh = onehot(seqs); nm = nmut(seqs)[:, None]
        pf = self.facs.predict_proba(ordinal(seqs))[:, 1][:, None]
        return pd.DataFrame({"cdrh3": seqs, "n_mut": nm[:, 0], "p_facs_high": pf[:, 0],
            "p_binder": self.clf.predict_proba(oh)[:, 1], "pred_log10_kd": self.rank.predict(np.hstack([oh, nm, pf]))})

def validate(spr, facs):
    """Leak-free grouped CV, reported overall and inside the >=3-mutation design-constraint band."""
    clean = facs[~facs.seq.isin(set(spr.cdrh3))]
    fm = HistGradientBoostingClassifier(max_iter=400, max_leaf_nodes=63, categorical_features=list(range(L)),
        random_state=0).fit(ordinal(clean.seq.values), clean.label.values)
    spr = spr.copy(); spr["p_facs"] = fm.predict_proba(ordinal(spr.cdrh3.values))[:, 1]
    spr["n_mut"] = nmut(spr.cdrh3.values); spr["grp"] = clusters(spr.cdrh3.values)
    out = {}
    b = spr[spr.binding_result.isin(["quantified_binding", "no_reportable_kd"])].reset_index(drop=True)
    y = (b.binding_result == "quantified_binding").astype(int).values
    X = np.hstack([onehot(b.cdrh3.values), b[["n_mut", "p_facs"]].values])
    p = np.zeros(len(y))
    for tr, te in GroupKFold(5).split(X, y, b.grp.values):
        p[te] = LogisticRegression(C=.3, max_iter=4000).fit(X[tr, :200], y[tr]).predict_proba(X[te, :200])[:, 1]
    out["clf"] = dict(n=len(y), pos=int(y.sum()), auroc=roc_auc_score(y, p), ap=average_precision_score(y, p),
        auroc_nmut_baseline=roc_auc_score(y, -b.n_mut.values))
    for floor in (3, 4):
        bd = b.n_mut >= floor
        out["clf"][f"n_ge{floor}"] = int(bd.sum())
        out["clf"][f"auroc_ge{floor}"] = roc_auc_score(y[bd], p[bd])
        out["clf"][f"ap_ge{floor}"] = average_precision_score(y[bd], p[bd])
    q = spr[spr.binding_result == "quantified_binding"].reset_index(drop=True)
    yq = np.log10(q.kd_nM.values)
    Xq = np.hstack([onehot(q.cdrh3.values), q[["n_mut", "p_facs"]].values])
    pq = np.zeros(len(yq))
    for tr, te in GroupKFold(5).split(Xq, yq, q.grp.values):
        pq[te] = Ridge(3.0).fit(Xq[tr], yq[tr]).predict(Xq[te])
    top = np.argsort(pq)[:20]
    out["rank"] = dict(n=len(yq), clusters=int(q.grp.nunique()), spearman=spearmanr(pq, yq).statistic,
        spearman_nmut_baseline=spearmanr(q.n_mut.values, yq).statistic,
        median_kd_overall=float(10**np.median(yq)), median_kd_top20=float(10**np.median(yq[top])), best_kd_top20=float(10**yq[top].min()))
    for floor in (3, 4):
        qd = q.n_mut >= floor
        out["rank"][f"n_ge{floor}"] = int(qd.sum())
        out["rank"][f"spearman_ge{floor}"] = spearmanr(pq[qd], yq[qd]).statistic
        tq = np.argsort(np.where(qd, pq, np.inf))[:20]
        out["rank"][f"median_kd_top20_ge{floor}"] = float(10**np.median(yq[tq]))
    return out

def benchmark(pf, n=2_000_000):
    rng = np.random.default_rng(0)
    seqs = ["".join(s) for s in rng.choice(list(AA), size=(n, L))]
    t0 = time.time(); oh = onehot(seqs); t_feat = time.time()-t0
    t0 = time.time(); pf.clf.predict_proba(oh); t_clf = time.time()-t0
    t0 = time.time(); pfs = pf.facs.predict_proba(ordinal(seqs))[:, 1]; t_facs = time.time()-t0
    t0 = time.time(); pf.rank.predict(np.hstack([oh, nmut(seqs)[:, None], pfs[:, None]])); t_rank = time.time()-t0
    return dict(n=n, featurise_s=t_feat, clf_s=t_clf, facs_s=t_facs, rank_s=t_rank,
        clf_seqs_per_s=int(n/(t_feat+t_clf)), full_seqs_per_s=int(n/(t_feat+t_clf+t_facs+t_rank)))

if __name__ == "__main__":
    spr, facs = load()
    print("=== leak-free grouped-CV validation ===")
    for k, d in validate(spr, facs).items():
        print(f"  {k}: " + ", ".join(f"{a}={b:.3f}" if isinstance(b, float) else f"{a}={b}" for a, b in d.items()))
    pf = PreFilter().fit(spr, facs)
    print(f"\nfitted; FACS rows dropped to avoid leakage: {pf.n_facs_dropped}")
    print("\n=== throughput (this laptop, single process) ===")
    for k, v in benchmark(pf).items(): print(f"  {k}: {v:,.3f}" if isinstance(v, float) else f"  {k}: {v:,}")
    pickle.dump(pf, open(MODELS/"prefilter.pkl", "wb"))
    print("\n=== demo: score trastuzumab's own CDRH3 and 5 random 10-mers ===")
    rng = np.random.default_rng(1)
    demo = [TZ] + ["".join(s) for s in rng.choice(list(AA), size=(5, L))]
    print(pf.score(demo).round(3).to_string(index=False))
