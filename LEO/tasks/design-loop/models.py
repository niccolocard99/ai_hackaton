"""The team's CDRH3 pre-filter in closed form, plus a cluster bootstrap for uncertainty.

Why this exists
---------------
NICCOLO/eda/prefilter.py fits two linear models on one-hot CDRH3: a logistic
classifier (does it bind at all) and a ridge ranker (how tightly). Both are
additive over positions, so a design's score is a sum of per-(position, residue)
table entries. Writing them as tables does two things:

  1. scoring 60 million candidates is a table lookup, not 60 million feature
     matrices, so we can enumerate a design space instead of sampling it;
  2. refitting on bootstrap resamples is cheap, so every prediction gets an
     uncertainty. That matters because the designs we want are the *argmax* of
     a surrogate over tens of millions of candidates. The maximum of many noisy
     estimates is biased upward (the optimizer's curse): the ranker's
     validated skill is Spearman +0.63, and the top of a 60M-way search is
     exactly where that error is largest. A residue seen in 2 training variants
     can carry a large coefficient on almost no evidence, and a bootstrap
     exposes it because the coefficient swings between resamples.

The tables reproduce PreFilter.score exactly; check_against_prefilter() asserts it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/home/hackathon/wt/leo-design-loop/NICCOLO/eda")
from prefilter import AA, L, TZ, PreFilter, clusters, nmut, onehot  # noqa: E402
from sklearn.linear_model import LogisticRegression, Ridge  # noqa: E402

TZ_CODES = np.array([AA.index(c) for c in TZ], np.int16)
_IDX = np.arange(L)[None, :]


class Tables:
    """score = intercept + sum_p W[p, residue_p]. Holds one ranker and one classifier."""

    def __init__(self, Wr: np.ndarray, br: float, Wc: np.ndarray, bc: float):
        self.Wr, self.br, self.Wc, self.bc = Wr, br, Wc, bc

    @staticmethod
    def from_fitted(rank: Ridge, clf: LogisticRegression) -> "Tables":
        Wr = rank.coef_[:200].reshape(L, 20).astype(np.float64).copy()
        n_mut_coef = float(rank.coef_[200])
        # The ranker has one extra feature, the substitution count. It is a count
        # of positions differing from trastuzumab, so it folds into the table as
        # a constant added to every non-wild-type residue.
        for p in range(L):
            for a in range(20):
                if a != TZ_CODES[p]:
                    Wr[p, a] += n_mut_coef
        Wc = clf.coef_[0].reshape(L, 20).astype(np.float64)
        return Tables(Wr, float(rank.intercept_), Wc, float(clf.intercept_[0]))

    def pred_log10_kd(self, codes: np.ndarray) -> np.ndarray:
        return self.br + self.Wr[_IDX, codes].sum(1)

    def logit(self, codes: np.ndarray) -> np.ndarray:
        return self.bc + self.Wc[_IDX, codes].sum(1)

    def p_binder(self, codes: np.ndarray) -> np.ndarray:
        return 1.0 / (1.0 + np.exp(-self.logit(codes)))


def fit_point(spr) -> tuple[Tables, PreFilter]:
    pf = PreFilter().fit(spr)
    return Tables.from_fitted(pf.rank, pf.clf), pf


def _resample_clusters(groups: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Row indices of a bootstrap resample that draws whole clusters with replacement.

    Rows within Hamming 2 of each other are near-duplicates, so resampling rows
    independently would understate the uncertainty by a factor of 2-3 (the
    ledger's own bootstrap finding, NICCOLO/eda/AGENT_GUIDE.md section 4).
    """
    uniq = np.unique(groups)
    pick = rng.choice(uniq, size=len(uniq), replace=True)
    members = {g: np.flatnonzero(groups == g) for g in uniq}
    return np.concatenate([members[g] for g in pick])


def bootstrap_tables(spr, B: int = 100, seed: int = 0) -> list[Tables]:
    rng = np.random.default_rng(seed)
    both = spr[spr.binding_result.isin(["quantified_binding", "no_reportable_kd"])].reset_index(drop=True)
    quant = spr[spr.binding_result == "quantified_binding"].reset_index(drop=True)
    g_both = clusters(both.cdrh3.tolist())
    g_quant = clusters(quant.cdrh3.tolist())
    y_both = (both.binding_result == "quantified_binding").astype(int).values
    x_both = onehot(both.cdrh3.values)
    x_quant = np.hstack([onehot(quant.cdrh3.values), nmut(quant.cdrh3.values)[:, None]])
    y_quant = np.log10(quant.kd_nM.values)
    out = []
    for _ in range(B):
        ib = _resample_clusters(g_both, rng)
        iq = _resample_clusters(g_quant, rng)
        if len(np.unique(y_both[ib])) < 2:
            continue
        clf = LogisticRegression(C=0.3, max_iter=4000).fit(x_both[ib], y_both[ib])
        rank = Ridge(3.0).fit(x_quant[iq], y_quant[iq])
        out.append(Tables.from_fitted(rank, clf))
    return out


def ensemble_predict(tabs: list[Tables], codes: np.ndarray, chunk: int = 50000, ref=None):
    """Per-candidate bootstrap mean/sd of predicted log10 KD and of the binder logit."""
    n = len(codes)
    km = np.empty(n); ks = np.empty(n); lm = np.empty(n); l10 = np.empty(n); pb = np.full(n, np.nan)
    kref = None if ref is None else np.array([t.pred_log10_kd(ref[None, :])[0] for t in tabs])
    for i in range(0, n, chunk):
        c = codes[i : i + chunk]
        K = np.stack([t.pred_log10_kd(c) for t in tabs], 1)
        Lg = np.stack([t.logit(c) for t in tabs], 1)
        km[i : i + chunk] = K.mean(1)
        if kref is not None:
            pb[i : i + chunk] = (K < kref[None, :]).mean(1)
        ks[i : i + chunk] = K.std(1)
        lm[i : i + chunk] = Lg.mean(1)
        l10[i : i + chunk] = np.quantile(Lg, 0.10, axis=1)
    return dict(kd_mean=km, kd_sd=ks, logit_mean=lm, logit_p10=l10,
                p_binder_p10=1 / (1 + np.exp(-l10)), p_beats_ref=pb)


def check_against_prefilter(tab: Tables, pf: PreFilter, loops: list[str]) -> float:
    """Max abs difference between the closed form and PreFilter.score. Must be ~0."""
    codes = np.array([[AA.index(c) for c in s] for s in loops], np.int16)
    ref = pf.score(loops)
    d1 = np.abs(tab.pred_log10_kd(codes) - ref["pred_log10_kd"].to_numpy()).max()
    d2 = np.abs(tab.p_binder(codes) - ref["p_binder"].to_numpy()).max()
    return float(max(d1, d2))


if __name__ == "__main__":
    from prefilter import load

    spr, _ = load()
    tab, pf = fit_point(spr)
    probe = ["WGGDGFYAMD", "WGGSGPYLFD", "WGSQGFYAYD", "WAASGPYLFE", "FGVEGFYEMI"]
    print("closed form vs PreFilter.score, max abs diff:", check_against_prefilter(tab, pf, probe))
    import time

    t0 = time.time()
    tabs = bootstrap_tables(spr, B=100, seed=0)
    print(f"bootstrap: {len(tabs)} fits in {time.time() - t0:.1f}s")
    codes = np.array([[AA.index(c) for c in s] for s in probe], np.int16)
    e = ensemble_predict(tabs, codes)
    for s, m, sd, p in zip(probe, e["kd_mean"], e["kd_sd"], e["p_binder_p10"]):
        print(f"  {s}  log10KD {m:+.3f} +/- {sd:.3f}   p_binder(p10) {p:.3f}")
