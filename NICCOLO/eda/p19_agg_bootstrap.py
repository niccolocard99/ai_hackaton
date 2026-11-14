"""Step 19: cluster bootstrap on the one finding that would change the config --
borda(pred_log10_kd + p_binder) vs pred_log10_kd alone, which is what niccolo_loop.yaml ships."""
import sys, pathlib
HERE = pathlib.Path(__file__).resolve().parent; REPO = HERE.parents[1]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(REPO / "JAHN/pipeline"))
import numpy as np, pandas as pd
from sklearn.linear_model import Ridge, LogisticRegression
from scipy.stats import spearmanr, rankdata
from prefilter import load, onehot, nmut, clusters, TZ

spr, _ = load()
q = spr[spr.binding_result == "quantified_binding"].reset_index(drop=True); q["n_mut"] = nmut(q.cdrh3.values)
y = np.log10(q.kd_nM.values); grp = clusters(q.cdrh3.tolist(), 2)
b = spr[spr.binding_result.isin(["quantified_binding", "no_reportable_kd"])].reset_index(drop=True)
b["n_mut"] = nmut(b.cdrh3.values); yb = (b.binding_result == "quantified_binding").astype(int).values
gb = clusters(b.cdrh3.tolist(), 2)
def oof(X, yy, g, kind, reps=8):
    acc = np.zeros(len(yy)); cnt = np.zeros(len(yy))
    for r in range(reps):
        rs = np.random.default_rng(r); u = np.unique(g); rs.shuffle(u)
        a = {x: i % 5 for i, x in enumerate(u)}; f = np.array([a[x] for x in g])
        for k in range(5):
            tr, te = np.where(f != k)[0], np.where(f == k)[0]
            if len(te) == 0: continue
            if kind == "clf":
                if len(np.unique(yy[tr])) < 2: continue
                acc[te] += LogisticRegression(C=.3, max_iter=4000).fit(X[tr], yy[tr]).predict_proba(X[te])[:, 1]
            else:
                acc[te] += Ridge(3.0).fit(X[tr], yy[tr]).predict(X[te])
            cnt[te] += 1
    return acc / np.maximum(cnt, 1)
pred = oof(np.hstack([onehot(q.cdrh3.values), q[["n_mut"]].values.astype(float)]), y, grp, "reg")
pb = pd.Series(oof(onehot(b.cdrh3.values), yb, gb, "clf"), index=b.cdrh3).reindex(q.cdrh3).to_numpy()

def borda(cols_asc):
    """mean of per-column ranks, 1 = best; each entry is (values, lower_is_better) as his _ranks does."""
    rs = [rankdata(v if low else -v) for v, low in cols_asc]
    return np.mean(rs, axis=0)

A = borda([(pred, True)])                       # what niccolo_loop.yaml ships
B = borda([(pred, True), (-pb, True)])          # adding p_binder (higher is better -> negate)
for lab, mask in [("all 179", np.ones(len(y), bool)), ("legal band >=4", (q.n_mut >= 4).to_numpy())]:
    yy, aa, bb, gg = y[mask], A[mask], B[mask], grp[mask]
    sA, sB = spearmanr(aa, yy).statistic, spearmanr(bb, yy).statistic
    tA, tB = 10**np.median(yy[np.argsort(aa)[:20]]), 10**np.median(yy[np.argsort(bb)[:20]])
    u = np.unique(gg); idx = {x: np.where(gg == x)[0] for x in u}; rs = np.random.default_rng(0)
    d, dt = [], []
    for _ in range(4000):
        p = np.concatenate([idx[x] for x in rs.choice(u, len(u), replace=True)])
        d.append(spearmanr(bb[p], yy[p]).statistic - spearmanr(aa[p], yy[p]).statistic)
        k = min(20, len(p))
        dt.append(np.median(yy[p][np.argsort(bb[p])[:k]]) - np.median(yy[p][np.argsort(aa[p])[:k]]))
    lo, hi = np.percentile(d, [2.5, 97.5]); tlo, thi = np.percentile(dt, [2.5, 97.5])
    v = "ADDS SIGNAL" if lo > 0 else ("hurts" if hi < 0 else "no difference (CI spans 0)")
    print(f"\n{lab} (n={mask.sum()})")
    print(f"  ranker alone (shipped config) spearman {sA:+.3f}  top20 median KD {tA:8.2f} nM")
    print(f"  + p_binder via borda          spearman {sB:+.3f}  top20 median KD {tB:8.2f} nM")
    print(f"  spearman difference {np.mean(d):+.3f}  95% CI [{lo:+.3f}, {hi:+.3f}]  -> {v}")
    print(f"  top20 median log10KD difference {np.mean(dt):+.3f}  95% CI [{tlo:+.3f}, {thi:+.3f}]"
          f"  -> {'better' if thi < 0 else ('worse' if tlo > 0 else 'no difference')}")
