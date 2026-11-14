"""Step 16: per-position effect on affinity, all 131 legal quantified variants, with FDR correction.
Tests EXP-008's claim that A106/M107/D102 are the positions winners use."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
TZ = "WGGDGFYAMD"
raw = pd.read_csv("data/trastuzumab_her2_spr_all_768.csv")
v = raw[(raw.cdrh3.str.len() == 10) & (raw.cdrh3 != TZ) & (raw.cdrh3 != "ARDRIVDDYYYYGMDV")].copy()
v["nmut"] = v.cdrh3.map(lambda s: sum(a != b for a, b in zip(s, TZ)))
q = v[v.binding_result == "quantified_binding"].copy()
L = q[q.nmut >= 4].copy(); L["y"] = np.log10(L.kd_nM)
print(f"n={len(L)} legal quantified variants, median KD {L.kd_nM.median():.1f} nM\n")
rows = []
for p in range(10):
    mut = L[L.cdrh3.str[p] != TZ[p]]; con = L[L.cdrh3.str[p] == TZ[p]]
    if len(mut) < 3 or len(con) < 3:
        rows.append(dict(pos=p+99, wt=TZ[p], n_mut=len(mut), n_con=len(con), dlog=np.nan, p=np.nan)); continue
    u = mannwhitneyu(mut.y, con.y)
    rows.append(dict(pos=p+99, wt=TZ[p], n_mut=len(mut), n_con=len(con),
                     med_mut=mut.kd_nM.median(), med_con=con.kd_nM.median(),
                     dlog=mut.y.median()-con.y.median(), p=u.pvalue))
r = pd.DataFrame(rows)
ok = r.p.notna()
# Benjamini-Hochberg
pv = r.loc[ok, "p"].values; order = np.argsort(pv); m = len(pv)
q_bh = np.empty(m); run = 1.0
for i in range(m-1, -1, -1):
    run = min(run, pv[order[i]]*m/(i+1)); q_bh[order[i]] = run
r.loc[ok, "q_BH"] = q_bh
r["effect"] = np.where(r.dlog > 0, "mutating HURTS", "mutating helps")
print("dlog = median log10 KD (mutated) - (conserved). Positive = mutating this position is WORSE.")
print(r.sort_values("p").to_string(index=False, float_format=lambda x: f"{x:.3f}"))
sig = r[(r.q_BH < 0.1) & ok]
print(f"\n  positions with BH q<0.10: {sorted(sig.pos.tolist()) if len(sig) else 'none'}")
print("\n=== is mutating A106/M107 more common among winners, or just common everywhere? ===")
for p, lab in [(7, "A106"), (8, "M107"), (3, "D102")]:
    print(f"  {lab}: mutated in {(L.cdrh3.str[p] != TZ[p]).mean():.0%} of ALL legal variants"
          f" | {(L.nsmallest(8,'kd_nM').cdrh3.str[p] != TZ[p]).mean():.0%} of the 8 tightest")
print("  -> if the two percentages are similar, 'mutated in 8 of 8 winners' carries no information.")
print("\n=== what DOES associate with tighter binding among legal variants? ===")
print(f"  n_mut vs log10 KD: spearman {spearmanr(L.nmut, L.y).statistic:+.3f} (p={spearmanr(L.nmut, L.y).pvalue:.4f})")
cons = [p for p in range(10) if r.loc[p, "dlog"] > 0]
L["n_core_hit"] = L.cdrh3.map(lambda s: sum(s[p] != TZ[p] for p in cons))
print(f"  substitutions at the 'mutating hurts' positions {[p+99 for p in cons]}: spearman {spearmanr(L.n_core_hit, L.y).statistic:+.3f} (p={spearmanr(L.n_core_hit, L.y).pvalue:.4f})")
