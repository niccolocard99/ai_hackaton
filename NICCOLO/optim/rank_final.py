#!/usr/bin/env python3
"""Final ordering of the 120 variants: the classifier orders, RF3 promotes.

Evidence behind the rule: Jahn's 273 SPR-labelled RF3 folds (JAHN/pipeline/data/calibration_rf3.csv: same model, same HER2 construct,
wild type ipTM 0.594 in both runs):
  RF3 ipTM AUROC 0.82 for >=4 substitutions, and all of it sits in the tail: the highest non-binder of 136 is 0.622 and all 7 variants
  above 0.70 are binders (0.622-0.70 is empty). Below 0.62 it says nothing (44% binders vs a 50% base rate).
  RF3 pose metrics (epitope_recall) reach AUROC 0.63 only: not used. Boltz-2 ipTM is useless there (AUROC 0.53, non-binders reach 0.96): not run.
Rule: tier 1 = RF3 ipTM >= IPTM_TIER1; tier 2 = the rest. Inside a tier, order by the team's frozen stage-1 criterion
(borda over pred_log10_kd and p_binder of the CURRENT NICCOLO/eda/prefilter.py, see rescore_current_prefilter.py).

    /home/hackathon/.venvs/abmodels/bin/python rank_final.py
"""
import itertools
import numpy as np, pandas as pd
pd.set_option("display.width", 230); pd.set_option("display.max_columns", 30)
OPT = "/home/hackathon/nick-work/optim"
IPTM_TIER1 = 0.70

s = pd.read_csv(f"{OPT}/out/rf3_structure_scores.csv"); s = s[s.role == "design"][["rank", "name", "iptm", "iptm_AH", "iptm_AL", "ipae", "pae_cdrh3_ep", "ep_recall", "lrmsd_xtal", "has_clash"]]
r = pd.read_csv(f"{OPT}/out/final_variants_rescored_current_prefilter.csv")
f = pd.read_csv(f"{OPT}/out/final_variants.csv")
d = r.merge(s, on="rank", how="left"); assert len(d) == 120 and d.iptm.notna().all() and not d.has_clash.any()
d["tier"] = np.where(d.iptm >= IPTM_TIER1, 1, 2)
d = d.sort_values(["tier", "borda_rank_now"]).reset_index(drop=True); d.insert(0, "final_rank", range(1, len(d) + 1))
d.to_csv(f"{OPT}/out/final_variants_ranked.csv", index=False)
f.set_index("sequence_id").loc[d.sequence_id].reset_index()[["sequence_id", "kind", "heavy_sequence", "light_sequence", "design_description"]] \
    .to_csv(f"{OPT}/out/final_variants_ranked_abguard.csv", index=False)       # abguard input format, same order: head -n N+1 gives the top N

print(f"tier 1 (RF3 ipTM >= {IPTM_TIER1}): {int((d.tier == 1).sum())} designs; tier 2: {int((d.tier == 2).sum())}")
show = ["final_rank", "tier", "window", "mutations", "p_binder", "pred_log10_kd_now", "borda_rank_now", "iptm", "iptm_AH", "ipae", "ep_recall", "lrmsd_xtal"]
print(d[show].head(20).round(3).to_string(index=False))
for n in (10, 20):
    t = d.head(n); W = [list(w) for w in t.window]
    subs = pd.Series([m for ms in t.mutations for m in ms.split(",")]).value_counts()
    ham = min(sum(a != b for a, b in zip(x, y)) for x, y in itertools.combinations(t.window, 2))
    print(f"top {n}: tier-1 {int((t.tier == 1).sum())}, min pairwise Hamming {ham}, most repeated substitutions {subs.head(4).to_dict()}, "
          f"min p_binder {t.p_binder.min():.3f}, n_mut {t.n_mut.value_counts().sort_index().to_dict()}")
print("wrote out/final_variants_ranked.csv and out/final_variants_ranked_abguard.csv")
