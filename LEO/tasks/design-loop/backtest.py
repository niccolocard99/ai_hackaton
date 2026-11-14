#!/usr/bin/env python3
"""Does the selection method pick better binders than the plain ranker, on variants it never saw?

Held-out clusters (single-linkage Hamming <= 2, as NICCOLO/eda/prefilter.validate does), so every
test variant is >= 3 mutations from every training variant. That is exactly the situation of a
legal design, which novelty forces to be >= 3 from everything measured.

Predictors, all trained only on the training folds:
  ridge     the team's ranker, point estimate
  lcb       ridge bootstrap mean + 1 sd (what the loop selects on)
  anchored  for each training binder s within D_ANCH of the test variant x:
              y_s + (ridge(x) - ridge(s)),  weighted 1/d^2
            i.e. keep a real measurement and use the model only for the difference
  knn       the same neighbours' measured KD with no model at all
  random    the floor
Reported: Spearman against measured log10 KD by distance to the nearest training variant, and the
median measured KD of the 5 best-ranked variants per fold, which is the quantity that matters when
only a handful get tested.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, "/home/hackathon/wt/leo-design-loop/NICCOLO/eda")
from prefilter import TZ, clusters, load, nmut, onehot  # noqa: E402

D_ANCH, B, K_PICK, REPEATS = 5, 30, 5, 20


def feats(seqs):
    return np.hstack([onehot(seqs), nmut(seqs)[:, None]])


def main():
    spr, _ = load()
    q = spr[spr.binding_result == "quantified_binding"].reset_index(drop=True)
    seqs, y = q.cdrh3.values, np.log10(q.kd_nM.values)
    X = feats(seqs)
    O = np.array([[ord(c) for c in s] for s in seqs])
    grp = clusters(list(seqs))
    uniq = np.unique(grp)
    print(f"{len(q)} measured binders in {len(uniq)} clusters; KD {10**y.min():.3f}-{10**y.max():.0f} nM")
    rows, picks = [], []
    rng = np.random.default_rng(0)
    for rep in range(REPEATS):
        fold_of = dict(zip(uniq, rng.permutation(len(uniq)) % 5))
        fold = np.array([fold_of[g] for g in grp])
        for f in range(5):
            tr, te = np.where(fold != f)[0], np.where(fold == f)[0]
            m = Ridge(3.0).fit(X[tr], y[tr])
            yh_tr, yh_te = m.predict(X[tr]), m.predict(X[te])
            # bootstrap sd over cluster resamples of the training rows
            g_tr, ub = grp[tr], np.unique(grp[tr])
            preds = []
            for _ in range(B):
                ids = rng.choice(ub, len(ub), replace=True)
                idx = np.concatenate([np.where(g_tr == g)[0] for g in ids])
                preds.append(Ridge(3.0).fit(X[tr][idx], y[tr][idx]).predict(X[te]))
            lcb = np.mean(preds, 0) + np.std(preds, 0)
            d = (O[te][:, None, :] != O[tr][None, :, :]).sum(2)       # (test, train) Hamming
            near = d.min(1)
            anch, knn = yh_te.copy(), np.full(len(te), y[tr].mean())
            for i in range(len(te)):
                nb = np.where(d[i] <= D_ANCH)[0]
                if len(nb):
                    w = 1.0 / d[i, nb] ** 2
                    anch[i] = np.sum(w * (y[tr][nb] + yh_te[i] - yh_tr[nb])) / w.sum()
                    knn[i] = np.sum(w * y[tr][nb]) / w.sum()
            pr = {"ridge": yh_te, "lcb": lcb, "anchored": anch, "knn": knn, "random": rng.random(len(te))}
            for i, t in enumerate(te):
                rows.append(dict(rep=rep, fold=f, near=int(near[i]), nmut=int(nmut([seqs[t]])[0]),
                                 y=y[t], **{k: v[i] for k, v in pr.items()}))
            for k, v in pr.items():
                top = te[np.argsort(v)[:K_PICK]]
                picks.append(dict(rep=rep, fold=f, method=k, med_kd=float(np.median(10 ** y[top])),
                                  hit1=float(np.mean(10 ** y[top] < 1.0))))
    R, P = pd.DataFrame(rows), pd.DataFrame(picks)
    R.to_csv(HERE / "out/backtest_rows.csv", index=False)
    meth = ["ridge", "lcb", "anchored", "knn", "random"]
    print(f"\nSpearman(prediction, measured log10 KD), mean over {REPEATS} repeats of 5-fold grouped CV")
    print("  (lower predicted = tighter, so a good predictor has a POSITIVE correlation)")
    for label, sel in [("all held-out", R.near >= 3), ("nearest train = 3", R.near == 3),
                       ("nearest train = 4", R.near == 4), ("nearest train >= 5", R.near >= 5),
                       ("in legal band, >= 4 subs from trastuzumab", R.nmut >= 4)]:
        sub = R[sel]
        if len(sub) < 15:
            continue
        s = {k: np.mean([spearmanr(g[k], g.y)[0] for _, g in sub.groupby("rep") if len(g) > 5]) for k in meth}
        print(f"  {label:<44} n/rep~{len(sub) // REPEATS:>3}  " + "  ".join(f"{k} {s[k]:+.3f}" for k in meth))
    print(f"\nTop {K_PICK} per fold: median measured KD (nM) of the picks, and share under 1 nM")
    g = P.groupby("method")[["med_kd", "hit1"]].mean().loc[meth]
    print(g.round(3).to_string())
    base = R.y.pipe(lambda s: float(np.median(10 ** s)))
    print(f"\nfor scale: median KD of all held-out binders = {base:.1f} nM")


if __name__ == "__main__":
    main()
