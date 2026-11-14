"""Step 20: did JAHN's selection fix pick better candidates, measured against real KD?

He reports the old select_diverse was starved by a cap that grew with the list, drawing 53 of 60 picks from
around the median of the pool. Tested here on the 131 legal variants whose KD is measured: rank them with
the pipeline's own borda(pred_log10_kd, p_binder) using OUT-OF-FOLD pre-filter scores, then run both
versions of select_diverse and compare the true affinity of what each returns.
"""
import sys, pathlib, importlib.util, math
from collections import Counter
HERE = pathlib.Path(__file__).resolve().parent; REPO = HERE.parents[1]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(REPO / "JAHN/pipeline"))
import numpy as np, pandas as pd, polars as pl
from sklearn.linear_model import Ridge, LogisticRegression
from scipy.stats import rankdata
from prefilter import load, onehot, nmut, clusters, TZ
from abpipe import _niccolo, novelty
from abpipe.base import Parent
from abpipe.iterate import select_diverse as select_new

spec = importlib.util.spec_from_file_location("iterate_old", "/tmp/iterate_old.py")
old_mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(old_mod)
select_old = old_mod.select_diverse

spr, _ = load(); spr = spr.copy(); spr["n_mut"] = nmut(spr.cdrh3.values)
b = spr[spr.binding_result.isin(["quantified_binding", "no_reportable_kd"])].reset_index(drop=True)
b["n_mut"] = nmut(b.cdrh3.values); yb = (b.binding_result == "quantified_binding").astype(int).values
gb = clusters(b.cdrh3.tolist(), 2)
q = spr[spr.binding_result == "quantified_binding"].reset_index(drop=True)
y_all = np.log10(q.kd_nM.values); gq = clusters(q.cdrh3.tolist(), 2)
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
pred = oof(np.hstack([onehot(q.cdrh3.values), q[["n_mut"]].values.astype(float)]), y_all, gq, "reg")
pb = pd.Series(oof(onehot(b.cdrh3.values), yb, gb, "clf"), index=b.cdrh3).reindex(q.cdrh3).to_numpy()

# the shipped ranking: borda over both pre-filter columns
score = (rankdata(pred) + rankdata(-pb)) / 2.0
keep = (q.n_mut >= 4).to_numpy()                      # the legal band: what we may actually submit
qq = q[keep].reset_index(drop=True); yq = y_all[keep]; sc = score[keep]
order = np.argsort(sc)
ranked_idx = order                                    # best first
raw = pd.read_csv("data/trastuzumab_her2_spr_all_768.csv")
row = raw[raw.cdrh3 == TZ].iloc[0]; H = row.vh_sequence; st = H.index(TZ)
parent = Parent(id="tz", heavy=H, light=row.vl_sequence, cdrh3="SR" + TZ + "Y", window=TZ, window_start=st)
ranked = pl.DataFrame({
    "sequence_id": [f"c{i}" for i in ranked_idx],
    "heavy_sequence": [H[:st] + qq.cdrh3[i] + H[st + 10:] for i in ranked_idx],
}).with_row_index("rank", offset=1)
print(f"ranked pool: {ranked.height} legal variants with measured KD, median {10**np.median(yq):.1f} nM")

def summarise(name, sel, n_req):
    ids = [int(s[1:]) for s in sel["sequence_id"].to_list()]
    yy = yq[ids]
    pool_ranks = [int(np.where(ranked_idx == i)[0][0]) + 1 for i in ids]
    wins = [qq.cdrh3[i] for i in ids]
    wt = TZ
    cnt = Counter((p, w[p]) for w in wins for p in range(10) if w[p] != wt[p])
    return dict(selector=name, asked=n_req, got=len(ids),
                median_kd=10 ** np.median(yy), best_kd=10 ** yy.min(),
                median_pool_rank=int(np.median(pool_ranks)), worst_pool_rank=max(pool_ranks),
                max_shared_sub=max(cnt.values()) if cnt else 0)
rows = []
for n_req in (10, 20, 30):
    rows.append(summarise("old (cap grows with list)", select_old(ranked, parent, n_req), n_req))
    rows.append(summarise("new (pool_frac + stable cap)", select_new(ranked, parent, n_req), n_req))
    rows.append(summarise("plain top-N by rank", ranked.head(n_req), n_req))
t = pd.DataFrame(rows)
print("\n=== what each selector actually returns (true KD, nM) ===")
print(t.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
print("\n  median_pool_rank: where in the ranking the picks come from; lower is better (1 = top of the pool)")
print("  max_shared_sub: most times any single substitution (e.g. A106E) appears in the returned list")

# ---------------------------------------------------------------- bootstrap the old-vs-new difference
print("\n=== cluster bootstrap on the difference (n_select=20) ===")
gql = gq[keep]
u = np.unique(gql); idx = {x: np.where(gql == x)[0] for x in u}
rs = np.random.default_rng(0)
d_old_new, got_new = [], []
for _ in range(400):
    pick = np.concatenate([idx[x] for x in rs.choice(u, len(u), replace=True)])
    o = np.argsort(sc[pick])                                  # rank this resample
    ids = pick[o]
    rk = pl.DataFrame({"sequence_id": [f"c{i}" for i in range(len(ids))],
                       "heavy_sequence": [H[:st] + qq.cdrh3[i] + H[st + 10:] for i in ids]}).with_row_index("rank", offset=1)
    so = select_old(rk, parent, 20); sn = select_new(rk, parent, 20)
    yo = np.array([yq[ids[int(s[1:])]] for s in so["sequence_id"].to_list()])
    yn = np.array([yq[ids[int(s[1:])]] for s in sn["sequence_id"].to_list()])
    if len(yo) == 0 or len(yn) == 0: continue
    d_old_new.append(np.median(yn) - np.median(yo))           # negative = new is tighter
    got_new.append(len(yn))
d = np.array(d_old_new); lo, hi = np.percentile(d, [2.5, 97.5])
print(f"  median log10 KD, new minus old: {d.mean():+.3f}  95% CI [{lo:+.3f}, {hi:+.3f}]")
print(f"  -> {'NEW IS BETTER (CI below 0)' if hi < 0 else ('old better' if lo > 0 else 'no difference (CI spans 0)')}")
print(f"  in linear terms: new selects a median {10**abs(d.mean()):.1f}x {'tighter' if d.mean()<0 else 'weaker'} than old")
print(f"  new returns {np.mean(got_new):.1f} of 20 on average (it shortens rather than pad)")
