"""Step 15: cross-check LEO's ledger claims against the FULL SPR-768 table.
EXP-008 / TH-019 ("1 of 21 legal variants beats trastuzumab") came from the 60-variant fold-calibration
subset, which was chosen, not sampled. This answers IDEA-017 with all 179 quantified variants."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import numpy as np, pandas as pd
TZ = "WGGDGFYAMD"
raw = pd.read_csv("data/trastuzumab_her2_spr_all_768.csv")

tzr = raw[raw.cdrh3 == TZ]
print("=== the baseline we are scored against: trastuzumab's own measured KD ===")
print(f"  {len(tzr)} replicate rows | values nM: {sorted(tzr.kd_nM.dropna().round(3).tolist())}")
print(f"  min {tzr.kd_nM.min():.3f} | median {tzr.kd_nM.median():.3f} | max {tzr.kd_nM.max():.3f}")
print("  EXP-008 used a single row, 0.152 nM. The 8 replicates span 6x, so the count of variants that")
print("  'beat trastuzumab' depends on which one you pick -- shown below for each.")

v = raw[(raw.cdrh3.str.len() == 10) & (raw.cdrh3 != TZ) & (raw.cdrh3 != "ARDRIVDDYYYYGMDV")].copy()
v["nmut"] = v.cdrh3.map(lambda s: sum(a != b for a, b in zip(s, TZ)))
q = v[v.binding_result == "quantified_binding"].copy()
print(f"\n=== full table: {len(q)} quantified 10-mer variants (EXP-008 saw 30 of them) ===")
print(q.groupby("nmut").kd_nM.agg(["size", "min", "median"]).round(3).to_string())
legal = q[q.nmut >= 4]
print(f"\n  legal (>=4 substitutions): n={len(legal)}  [EXP-008 had 21]")
for lab, base in [("single row used by EXP-008", 0.152), ("median of 8 replicates", tzr.kd_nM.median()),
                  ("best replicate", tzr.kd_nM.min())]:
    w = legal[legal.kd_nM < base]
    print(f"    vs {lab:28s} ({base:.3f} nM): {len(w):2d} of {len(legal)} = {100*len(w)/len(legal):4.1f}%  best {legal.kd_nM.min():.3f} nM")
print(f"\n  EXP-008's existence proof holds, but the rate is a RATE now:")
w = legal[legal.kd_nM < tzr.kd_nM.median()]
print(f"    {len(w)} of {len(legal)} legal variants beat the replicate median; EXP-008 inferred 1 of 21 (4.8%)")
if len(w): print(w[["cdrh3", "nmut", "kd_nM"]].sort_values("kd_nM").to_string(index=False))

print("\n=== position convergence (EXP-008 used the 8 tightest of its subset) ===")
print("  Recomputed on ALL legal quantified variants, split by whether they beat the set median:")
tight = legal[legal.kd_nM <= legal.kd_nM.median()]
loose = legal[legal.kd_nM > legal.kd_nM.median()]
rows = []
for p in range(10):
    tm = (tight.cdrh3.str[p] != TZ[p]).mean(); lm = (loose.cdrh3.str[p] != TZ[p]).mean()
    rows.append(dict(pos=f"{p+99}", wt=TZ[p], mut_in_tight=f"{tm:.0%}", mut_in_loose=f"{lm:.0%}", diff=tm-lm))
pc = pd.DataFrame(rows).sort_values("diff", ascending=False)
print(pc.to_string(index=False, float_format=lambda x: f"{x:+.3f}"))
print(f"  (tight n={len(tight)}, loose n={len(loose)}; positive diff = mutated more often in the tighter half)")

print("\n=== does 'mutate A106/M107' actually predict tighter binding? (all legal variants) ===")
for p, lab in [(7, "A106"), (8, "M107"), (3, "D102"), (9, "D108"), (5, "F104")]:
    a = legal[legal.cdrh3.str[p] != TZ[p]].kd_nM
    b = legal[legal.cdrh3.str[p] == TZ[p]].kd_nM
    if len(a) < 3 or len(b) < 3: print(f"  {lab}: too few ({len(a)}/{len(b)})"); continue
    from scipy.stats import mannwhitneyu
    u = mannwhitneyu(np.log10(a), np.log10(b))
    print(f"  {lab}: mutated n={len(a)} median {a.median():8.2f} nM | conserved n={len(b)} median {b.median():8.2f} nM | p={u.pvalue:.3f}")
