#!/usr/bin/env python3
"""Re-score the 120 submitted CDRH3 windows with the CURRENT origin/main prefilter and compare with the scores used for selection.

gen_score.py selected with NICCOLO/eda/prefilter.py as of da498b8 (FACS score is a feature of the ranker). origin/main defaults to
use_facs=False. p_binder should not move; pred_log10_kd can. This also ranks the 120 by the team's frozen stage-1 criterion
(borda over pred_log10_kd and p_binder, NICCOLO/PIPELINE.md) and measures the distance to the nearest SPR training window, because the
ranker is documented as blind beyond ~4 substitutions from its training set and as extrapolating below log10 KD -1.02.

    /home/hackathon/.venvs/abmodels/bin/python rescore_current_prefilter.py     # needs eda_current/ = prefilter.py from origin/main + data symlink
"""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr
OPT = "/home/hackathon/nick-work/optim"
os.chdir(f"{OPT}/eda_current"); sys.path.insert(0, ".")
from prefilter import PreFilter, load            # noqa: E402  (origin/main version, use_facs=False by default)
pd.set_option("display.width", 220)

spr, _ = load(); pf = PreFilter().fit(spr)
v = pd.read_csv(f"{OPT}/out/final_variants_scored.csv")
s = pf.score(v.window.tolist())
v["p_binder_now"] = s.p_binder.values; v["pred_log10_kd_now"] = s.pred_log10_kd.values
print(f"current prefilter use_facs={pf.use_facs}; fitted on {len(spr)} SPR variants")
print(f"p_binder   max |now - then| = {np.abs(v.p_binder_now - v.p_binder).max():.2e}")
d = v.pred_log10_kd_now - v.pred_log10_kd
print(f"pred_log10_kd  Spearman(then, now) = {spearmanr(v.pred_log10_kd, v.pred_log10_kd_now).statistic:.3f}; shift now-then: "
      f"mean {d.mean():+.3f}, min {d.min():+.3f}, max {d.max():+.3f} (log10 nM)")
print(f"gate p_binder >= 0.80 under the current model: {(v.p_binder_now >= .8).sum()}/{len(v)}")
print(f"pred_log10_kd now: median {v.pred_log10_kd_now.median():.2f}, range {v.pred_log10_kd_now.min():.2f}..{v.pred_log10_kd_now.max():.2f}; "
      f"{(v.pred_log10_kd_now < -1.02).sum()} below -1.02 (tightest measured SPR KD 0.095 nM = extrapolation)")

# distance to the nearest SPR training window (Hamming), ranker's stated validity region
tr = np.array([list(w) for w in spr.cdrh3]); W = np.array([list(w) for w in v.window])
v["min_hamming_to_spr_train"] = [int((tr != w).sum(1).min()) for w in W]
print("min Hamming to the SPR training windows:", v.min_hamming_to_spr_train.value_counts().sort_index().to_dict())

# team's frozen stage-1 criterion: borda over both pre-filter columns (rank of pred_log10_kd ascending, p_binder descending)
r1 = v.pred_log10_kd_now.rank(method="average"); r2 = (-v.p_binder_now).rank(method="average")
v["borda_score_now"] = (r1 + r2) / 2; v["borda_rank_now"] = v.borda_score_now.rank(method="first").astype(int)
top = v.nsmallest(20, "borda_rank_now")
print(f"\nborda top-20 (current prefilter) overlaps my first 20 (old ranking) in {len(set(top['rank']) & set(range(1, 21)))} designs; "
      f"median rank of the borda top-20 in my list = {int(top['rank'].median())}")
cols = ["sequence_id", "rank", "window", "mutations", "n_mut", "p_binder", "pred_log10_kd", "p_binder_now", "pred_log10_kd_now",
        "min_hamming_to_spr_train", "borda_score_now", "borda_rank_now"]
v[cols].to_csv(f"{OPT}/out/final_variants_rescored_current_prefilter.csv", index=False)
print("wrote out/final_variants_rescored_current_prefilter.csv")
print(top[["rank", "window", "mutations", "p_binder_now", "pred_log10_kd_now", "min_hamming_to_spr_train", "borda_rank_now"]].round(3).to_string(index=False))
