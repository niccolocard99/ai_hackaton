"""Step 18: how effective is JAHN's scoring pipeline, measured against real KD?

Runs his scorers and aggregators on the 179 trastuzumab variants with a measured KD and asks whether the
composite ranking beats the pre-filter alone. The pre-filter columns are supplied OUT-OF-FOLD (grouped CV,
clusters held out) so the comparison is not rigged in its favour -- his scorers see no KD at all.
"""
import sys, pathlib
HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(REPO / "JAHN/pipeline"))
import numpy as np, pandas as pd, polars as pl
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.model_selection import GroupKFold
from scipy.stats import spearmanr
from prefilter import load, onehot, nmut, clusters, TZ
from abpipe.base import Context, Parent
from abpipe.scorers.parent_distance import ParentDistance
from abpipe.scorers.nearest_known_binder import NearestKnownBinder
from abpipe.aggregators.ranking import Borda, WeightedSum, Pareto

D = HERE / "data/abpipe_test"; C = D / "_eval_cache"; C.mkdir(parents=True, exist_ok=True)
spr_raw = pd.read_csv(D / "her2-trastuzumab/trastuzumab_her2_spr_all_768.csv")
row = spr_raw[spr_raw.cdrh3 == TZ].iloc[0]
H, Lc = row.vh_sequence, row.vl_sequence
st = H.index(TZ)
parent = Parent(id="trastuzumab", heavy=H, light=Lc, cdrh3="SR" + TZ + "Y", window=TZ, window_start=st)
ctx = Context(data_dir=D, repo_root=REPO, run_dir=C, cache_dir=C, seed=0, parent=parent)

spr, _ = load()
q = spr[spr.binding_result == "quantified_binding"].reset_index(drop=True)
q["n_mut"] = nmut(q.cdrh3.values)
y = np.log10(q.kd_nM.values)
grp = clusters(q.cdrh3.tolist(), 2)
print(f"evaluation set: {len(q)} variants with measured KD, {len(np.unique(grp))} sequence clusters")

# ---- out-of-fold pre-filter columns (his pipeline would get in-sample ones; we do not) ----
oh = onehot(q.cdrh3.values); nm = q[["n_mut"]].values.astype(float)
Xr = np.hstack([oh, nm])
b = spr[spr.binding_result.isin(["quantified_binding", "no_reportable_kd"])].reset_index(drop=True)
b["n_mut"] = nmut(b.cdrh3.values)
yb = (b.binding_result == "quantified_binding").astype(int).values
gb = clusters(b.cdrh3.tolist(), 2); ohb = onehot(b.cdrh3.values)
def oof_reg(X, yy, g, reps=8):
    acc = np.zeros(len(yy)); cnt = np.zeros(len(yy))
    for r in range(reps):
        rs = np.random.default_rng(r); u = np.unique(g); rs.shuffle(u)
        a = {x: i % 5 for i, x in enumerate(u)}; f = np.array([a[x] for x in g])
        for k in range(5):
            tr, te = np.where(f != k)[0], np.where(f == k)[0]
            if len(te) == 0: continue
            acc[te] += Ridge(3.0).fit(X[tr], yy[tr]).predict(X[te]); cnt[te] += 1
    return acc / np.maximum(cnt, 1)
def oof_clf(X, yy, g, reps=8):
    acc = np.zeros(len(yy)); cnt = np.zeros(len(yy))
    for r in range(reps):
        rs = np.random.default_rng(r); u = np.unique(g); rs.shuffle(u)
        a = {x: i % 5 for i, x in enumerate(u)}; f = np.array([a[x] for x in g])
        for k in range(5):
            tr, te = np.where(f != k)[0], np.where(f == k)[0]
            if len(te) == 0 or len(np.unique(yy[tr])) < 2: continue
            acc[te] += LogisticRegression(C=.3, max_iter=4000).fit(X[tr], yy[tr]).predict_proba(X[te])[:, 1]; cnt[te] += 1
    return acc / np.maximum(cnt, 1)
pred_kd_oof = oof_reg(Xr, y, grp)
pb_all = oof_clf(ohb, yb, gb)
pb_oof = pd.Series(pb_all, index=b.cdrh3).reindex(q.cdrh3).to_numpy()
print("out-of-fold pre-filter columns built (clusters held out)")

# ---- his scorers (they never see KD) ----
heavy = [H[:st] + w + H[st + 10:] for w in q.cdrh3]
df = pl.DataFrame({"sequence_id": [f"s{i}" for i in range(len(q))], "heavy_sequence": heavy,
                   "light_sequence": [Lc] * len(q)})
pd_out = ParentDistance().score(df, ctx)
print("parent_distance done")
nb_out = NearestKnownBinder(k=25).score(df, ctx)
print("nearest_known_binder done (524k FACS refs)")

scored = df.join(pd_out, on="sequence_id").join(nb_out, on="sequence_id").with_columns([
    pl.Series("score__p_binder", pb_oof), pl.Series("score__pred_log10_kd", pred_kd_oof)])
directions = {"score__parent_distance": "lower", "score__binder_distance": "lower",
              "score__neighbor_binder_frac": "higher", "score__p_binder": "higher",
              "score__pred_log10_kd": "lower"}

def report(name, series, lower_is_better):
    s = np.asarray(series, float)
    sign = 1.0 if lower_is_better else -1.0      # we want the score to increase with log KD
    r = spearmanr(sign * s, y)
    order = np.argsort(sign * s)
    top = order[:20]
    return dict(signal=name, spearman=r.statistic, p=r.pvalue,
                top20_median_kd=10 ** np.median(y[top]), top20_best_kd=10 ** y[top].min())

rows = [report(c, scored[c].to_numpy(), directions[c] == "lower") for c in directions]
print("\n=== individual signals vs measured log10 KD (n=%d) ===" % len(q))
t = pd.DataFrame(rows).sort_values("spearman", ascending=False)
print(t.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
print(f"  set median KD {10**np.median(y):.1f} nM; a random top-20 would median near that")

print("\n=== aggregators: does combining beat the pre-filter alone? ===")
combos = {
 "pre-filter ranker ALONE": ["score__pred_log10_kd"],
 "his 3 scorers, NO pre-filter": ["score__parent_distance", "score__binder_distance", "score__neighbor_binder_frac"],
 "everything": list(directions),
 "ranker + p_binder": ["score__pred_log10_kd", "score__p_binder"],
 "ranker + parent_distance": ["score__pred_log10_kd", "score__parent_distance"],
 "ranker + neighbor_binder_frac": ["score__pred_log10_kd", "score__neighbor_binder_frac"],
}
out = []
for label, cols in combos.items():
    for agg_name, agg in (("borda", Borda(columns=tuple(cols))), ("weighted_sum", WeightedSum(columns=tuple(cols))),
                          ("pareto", Pareto(columns=tuple(cols)))):
        if len(cols) == 1 and agg_name != "borda":
            continue
        res = agg.combine(scored, directions, ctx)
        idx = [int(s[1:]) for s in res["sequence_id"].to_list()]
        ranked_y = y[idx]
        out.append(dict(combo=label, aggregator=agg_name, n_cols=len(cols),
                        spearman=spearmanr(np.arange(len(ranked_y)), ranked_y).statistic,
                        top20_median_kd=10 ** np.median(ranked_y[:20]), top20_best_kd=10 ** ranked_y[:20].min()))
a = pd.DataFrame(out)
print(a.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

print("\n=== inside the legal band (>=4 substitutions), which is what we may submit ===")
mask = (q.n_mut >= 4).to_numpy()
sub = scored.filter(pl.Series(mask)); ysub = y[mask]
print(f"  n={mask.sum()}")
for label, cols in combos.items():
    res = Borda(columns=tuple(cols)).combine(sub, directions, ctx)
    idx = [int(s[1:]) for s in res["sequence_id"].to_list()]
    pos = {int(s[1:]): k for k, s in enumerate(res["sequence_id"].to_list())}
    yy = y[np.array(idx)]
    print(f"  {label:32s} spearman {spearmanr(np.arange(len(yy)), yy).statistic:+.3f}  top20 median KD {10**np.median(yy[:20]):8.2f} nM")
