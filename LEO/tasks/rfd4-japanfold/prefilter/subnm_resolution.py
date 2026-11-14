"""Does the pre-filter have resolution in the regime the hackathon is scored in?

EXP-007 reports grouped-CV Spearman +0.628 across the whole >=4 band, where
measured KD spans 0.023 nM to >100 nM. Rank skill over three orders of magnitude
says nothing about resolution INSIDE the sub-nanomolar cluster where
trastuzumab (0.152 nM) sits -- and that cluster is what we are scored against.

Two questions, both leak-free (grouped CV, clusters held out together, exactly
as validate() does it):
  1. Over the full range, does it rank? (expect yes, ~+0.63)
  2. Restricted to measured KD < 1 nM, does it still rank? (the real question)
Also reports the predicted spread, to show whether Ridge(3.0) has shrunk the
predictions into a band narrower than the differences we need to resolve.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold

sys.path.insert(0, str(Path(__file__).parent))
from prefilter import clusters, load, nmut, onehot, ordinal  # noqa: E402

spr, facs = load()
clean = facs[~facs.seq.isin(set(spr.cdrh3))]
fm = HistGradientBoostingClassifier(
    max_iter=400, max_leaf_nodes=63, categorical_features=list(range(10)), random_state=0
).fit(ordinal(clean.seq.values), clean.label.values)

spr = spr.copy()
spr["p_facs"] = fm.predict_proba(ordinal(spr.cdrh3.values))[:, 1]
spr["n_mut"] = nmut(spr.cdrh3.values)
spr["grp"] = clusters(spr.cdrh3.values)

q = spr[spr.binding_result == "quantified_binding"].reset_index(drop=True)
y = np.log10(q.kd_nM.values)
X = np.hstack([onehot(q.cdrh3.values), q[["n_mut", "p_facs"]].values])
pred = np.zeros(len(y))
for tr, te in GroupKFold(5).split(X, y, q.grp.values):
    pred[te] = Ridge(3.0).fit(X[tr], y[tr]).predict(X[te])
q["pred_kd_nM"] = 10**pred
q["pred_log10"] = pred

print(f"quantified variants: {len(q)}   clusters: {q.grp.nunique()}")
print(f"measured KD range: {q.kd_nM.min():.3f} - {q.kd_nM.max():.1f} nM")
print(f"PREDICTED KD range: {q.pred_kd_nM.min():.3f} - {q.pred_kd_nM.max():.1f} nM")
print(
    f"  measured spans {np.log10(q.kd_nM.max()/q.kd_nM.min()):.2f} log10 units; "
    f"predicted spans only {pred.max()-pred.min():.2f}"
)

print("\n=== rank skill by regime (leak-free grouped CV) ===")
print(f"{'regime':<26}{'n':>5}{'spearman':>11}{'median |err| (log10)':>23}")
for label, mask in [
    ("full range", np.ones(len(q), bool)),
    (">=4 subs (legal band)", (q.n_mut >= 4).values),
    ("measured < 10 nM", (q.kd_nM < 10).values),
    ("measured < 1 nM", (q.kd_nM < 1).values),
    ("<1 nM AND >=4 subs", ((q.kd_nM < 1) & (q.n_mut >= 4)).values),
]:
    n = int(mask.sum())
    if n < 4:
        print(f"{label:<26}{n:>5}{'(too few)':>11}")
        continue
    rho = spearmanr(pred[mask], y[mask]).statistic
    err = float(np.median(np.abs(pred[mask] - y[mask])))
    print(f"{label:<26}{n:>5}{rho:>11.3f}{err:>23.2f}")

sub = q[q.kd_nM < 1].sort_values("kd_nM")
print(f"\n=== every measured sub-nanomolar variant (n={len(sub)}) ===")
print(
    sub[["cdrh3", "n_mut", "kd_nM", "pred_kd_nM"]]
    .rename(columns={"kd_nM": "measured_nM", "pred_kd_nM": "predicted_nM"})
    .round(3)
    .to_string(index=False)
)
err_sub = np.median(np.abs(sub.pred_log10 - np.log10(sub.measured_nM if "measured_nM" in sub else sub.kd_nM)))
print(f"\nmedian absolute error in the sub-nM regime: {err_sub:.2f} log10 units "
      f"= {10**err_sub:.1f}x")
print("A predicted edge smaller than that factor is not a result.")
