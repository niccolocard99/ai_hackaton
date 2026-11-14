"""Step 21: score LEO's 20 final designs with the frozen stage-1 models.

Two things this must be honest about:
  1. The pre-filter scores the 10-residue CDRH3 window ONLY. It is blind to the framework and light-chain
     substitutions these designs carry ("addons": D28K on the VL, fwTm on the VH framework).
  2. The ranker interpolates. Its signal dies ~5 substitutions from its 179 training variants, so the
     minimum Hamming distance from each design to that training set is reported alongside every score.
"""
import sys, pathlib
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr
from prefilter import PreFilter, load, TZ, L

TZ_VH = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYC"
         "SRWGGDGFYAMDYWGQGTLVTVSS")
TZ_VL = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
WIN = TZ_VH.index(TZ)

d = pd.read_csv("data/leo_FINAL.csv")
print(f"{len(d)} designs from LEO's FINAL.csv")

# ---- what the pre-filter can and cannot see ----
d["vl_changed"] = d.light_sequence != TZ_VL
def fw_changes(h):
    if len(h) != len(TZ_VH): return None
    return sum(1 for i in range(len(h)) if not (WIN <= i < WIN + L) and h[i] != TZ_VH[i])
d["vh_fw_subs"] = d.heavy_sequence.map(fw_changes)
d["window"] = d.heavy_sequence.str[WIN:WIN + L]
assert (d.window == d.cdrh3).all(), "the CDRH3 column should match the window in the heavy chain"
print(f"  light chain changed: {d.vl_changed.sum()} of {len(d)}")
print(f"  heavy framework substitutions outside the window: {d.vh_fw_subs.sum()} across {(d.vh_fw_subs>0).sum()} designs")
print("  -> the pre-filter sees NONE of these; it reads the 10-residue window only")

# ---- score the window with the frozen models ----
spr, _ = load()
pf = PreFilter().fit(spr)
s = pf.score(d.window.tolist())
d = pd.concat([d, s.drop(columns=["cdrh3"]).add_prefix("pf_")], axis=1)

# ---- distance to the ranker's training set (the interpolation check) ----
train = spr[spr.binding_result == "quantified_binding"].cdrh3.tolist()
A = np.frombuffer("".join(train).encode(), np.uint8).reshape(-1, L)
B = np.frombuffer("".join(d.window).encode(), np.uint8).reshape(-1, L)
d["min_ham_to_train"] = (B[:, None, :] != A[None, :, :]).sum(2).min(1)

# ---- the frozen stage-1 ranking: borda over both pre-filter columns ----
d["pf_borda"] = (rankdata(d.pf_pred_log10_kd) + rankdata(-d.pf_p_binder)) / 2
d["pf_rank"] = rankdata(d.pf_borda, method="ordinal").astype(int)

cols = ["rank", "cdrh3", "addons", "n_sub", "kd_mean", "p_beats_tz", "pf_p_binder", "pf_pred_log10_kd",
        "pf_rank", "min_ham_to_train", "vl_changed", "vh_fw_subs"]
print("\n=== scored (LEO's rank vs the frozen pipeline's rank) ===")
print(d[cols].to_string(index=False, float_format=lambda x: f"{x:.3f}"))

print("\n=== do the two rankings agree? ===")
r = spearmanr(d["rank"], d.pf_rank)
print(f"  Spearman(LEO rank, pipeline rank) = {r.statistic:+.3f}  p={r.pvalue:.3f}")
r2 = spearmanr(d.kd_mean, d.pf_pred_log10_kd)
print(f"  Spearman(LEO kd_mean, pipeline pred_log10_kd) = {r2.statistic:+.3f}  p={r2.pvalue:.4f}")

print("\n=== blindness demonstration: designs sharing a window but differing in addons ===")
g = d.groupby("cdrh3").agg(n=("cdrh3", "size"), addons=("addons", lambda x: " | ".join(x)),
                           pf=("pf_pred_log10_kd", lambda x: round(x.iloc[0], 4)),
                           leo_spread=("kd_mean", lambda x: round(x.max() - x.min(), 3)))
dup = g[g.n > 1]
print(dup.to_string() if len(dup) else "  every design has a distinct window, so no direct demonstration here")

print("\n=== interpolation warning ===")
print(d.min_ham_to_train.value_counts().sort_index().rename("designs").to_string())
far = (d.min_ham_to_train >= 5).sum()
print(f"  designs >=5 substitutions from the nearest training variant: {far} of {len(d)}")
print("  the ranker's measured signal at that separation is Spearman -0.057, i.e. none.")
print("\n=== extrapolation warning (LOOP.md threshold) ===")
print(f"  predicted log10 KD below -1.02 (tighter than anything measured in the legal band): "
      f"{(d.pf_pred_log10_kd < -1.02).sum()} of {len(d)}")
print(f"  LEO's kd_mean below -1.02: {(d.kd_mean < -1.02).sum()} of {len(d)}")
d.to_csv("figs/leo_FINAL_scored.csv", index=False)
print("\nwritten: figs/leo_FINAL_scored.csv")
