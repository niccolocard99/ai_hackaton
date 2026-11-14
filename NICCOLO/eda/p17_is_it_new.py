"""Step 17: is the per-position result new information, or a restatement of what the ranker already learned?

Two tests:
  A  does the ranker's own coefficient structure already imply the same per-position directions?
  B  does adding a 'costly positions' counter improve the ranker out-of-fold? (cluster bootstrap)
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.compose import ColumnTransformer
from scipy.stats import spearmanr, mannwhitneyu
from prefilter import load, onehot, nmut, clusters
AA = "ACDEFGHIKLMNPQRSTVWY"; TZ = "WGGDGFYAMD"

spr, _ = load(); spr = spr.copy(); spr["n_mut"] = nmut(spr.cdrh3.values)
q = spr[spr.binding_result == "quantified_binding"].reset_index(drop=True)
y = np.log10(q.kd_nM.values)
L = q[q.n_mut >= 4].reset_index(drop=True); yL = np.log10(L.kd_nM.values)

print("=== A. does the ranker already know the per-position directions? ===")
# fit the production ranker on all 179 and read its implied cost of mutating each position
m = Ridge(3.0).fit(np.hstack([onehot(q.cdrh3.values), q[["n_mut"]].values.astype(float)]), y)
W = m.coef_[:200].reshape(10, 20)
imp = []
for p in range(10):
    wt = W[p, AA.index(TZ[p])]
    others = np.delete(W[p], AA.index(TZ[p]))
    imp.append(others.mean() - wt)          # >0 means the model thinks mutating here raises log KD (worse)
# independent Mann-Whitney effect from p16, recomputed here
mw = []
for p in range(10):
    a = L[L.cdrh3.str[p] != TZ[p]]; b = L[L.cdrh3.str[p] == TZ[p]]
    mw.append(np.log10(a.kd_nM).median() - np.log10(b.kd_nM).median() if len(a) >= 3 and len(b) >= 3 else np.nan)
t = pd.DataFrame({"pos": [p+99 for p in range(10)], "wt": list(TZ),
                  "ranker_implied_cost": imp, "measured_dlog": mw})
t["ranker_says"] = np.where(t.ranker_implied_cost > 0, "hurts", "helps")
t["data_says"] = np.where(t.measured_dlog > 0, "hurts", "helps")
t["agree"] = t.ranker_says == t.data_says
print(t.to_string(index=False, float_format=lambda x: f"{x:+.3f}"))
r = spearmanr(t.ranker_implied_cost, t.measured_dlog)
print(f"\n  Spearman(ranker implied cost, measured dlog) = {r.statistic:+.3f}  p={r.pvalue:.3f}")
print(f"  directional agreement: {t.agree.sum()}/10 positions")

print("\n=== B. does a 'costly positions' counter add anything out-of-fold? ===")
costly = [p for p in range(10) if t.measured_dlog[p] > 0]
print(f"  costly positions (measured): {[p+99 for p in costly]}")
def pipe(n_oh, dim):
    if n_oh >= dim: return Ridge(3.0)
    return make_pipeline(ColumnTransformer([("oh", "passthrough", list(range(n_oh))),
        ("num", StandardScaler(), list(range(n_oh, dim)))]), Ridge(3.0))
def folds(g, k, rs):
    u = np.unique(g); rs.shuffle(u); a = {x: i % k for i, x in enumerate(u)}
    f = np.array([a[x] for x in g]); return [(np.where(f != i)[0], np.where(f == i)[0]) for i in range(k)]
def oof(X, n_oh, g, yy, reps=12):
    acc = np.zeros(len(yy)); cnt = np.zeros(len(yy))
    for rr in range(reps):
        for tr, te in folds(g, 5, np.random.default_rng(rr)):
            if len(te) == 0: continue
            acc[te] += pipe(n_oh, X.shape[1]).fit(X[tr], yy[tr]).predict(X[te]); cnt[te] += 1
    return acc/np.maximum(cnt, 1)
g = clusters(q.cdrh3.tolist(), 2)
oh = onehot(q.cdrh3.values); nm = q[["n_mut"]].values.astype(float)
core = np.array([[sum(s[p] != TZ[p] for p in costly)] for s in q.cdrh3], float)
cur = oof(np.hstack([oh, nm]), oh.shape[1], g, y)
new = oof(np.hstack([oh, nm, core]), oh.shape[1], g, y)
only = oof(np.hstack([nm, core]), 0, g, y)
u = np.unique(g); idx = {x: np.where(g == x)[0] for x in u}; rs = np.random.default_rng(0)
d = []
for _ in range(3000):
    bb = np.concatenate([idx[x] for x in rs.choice(u, len(u), replace=True)])
    d.append(spearmanr(new[bb], y[bb]).statistic - spearmanr(cur[bb], y[bb]).statistic)
lo, hi = np.percentile(d, [2.5, 97.5])
print(f"  ranker (one-hot + n_mut)            OOF spearman {spearmanr(cur, y).statistic:+.3f}")
print(f"  + costly-position counter           OOF spearman {spearmanr(new, y).statistic:+.3f}")
print(f"  counter + n_mut only (no one-hot)    OOF spearman {spearmanr(only, y).statistic:+.3f}")
print(f"  difference {np.mean(d):+.3f}  95% CI [{lo:+.3f}, {hi:+.3f}]  -> "
      f"{'adds signal' if lo > 0 else ('hurts' if hi < 0 else 'no difference (CI spans 0)')}")
print("\n  reminder: the +0.505 quoted in LEDGER_CROSSCHECK was in-sample and univariate, on the 131")
print("  legal variants only. It is not comparable to an out-of-fold ranker score.")
