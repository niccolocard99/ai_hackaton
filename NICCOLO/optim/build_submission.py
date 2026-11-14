#!/usr/bin/env python3
"""selected CDRH3 windows -> paired VH/VL designs -> OFFICIAL check_optimisation.py -> JapanFold payloads (prepared, NOT submitted).

    /home/hackathon/.venvs/abmodels/bin/python build_submission.py

Outputs (out/ and jobs/):
  out/designs_candidates.csv        abguard input format for every selected window
  out/check_fast.csv                check_optimisation.py (CDRH3 reference only)
  out/check_provided_*.csv          check_optimisation.py --provided <SPR table>  (whole VH+VL vs every provided antibody)
  out/final_variants.csv            abguard-format submission table: the first N_FINAL designs that pass BOTH checks
  out/final_variants_scored.csv     same + classifier scores + check columns
  jobs/{rf3,boltz2}_payload.json    POST /v1/predictions bodies (H/L/A chains, same layout as Jahn's calibration folds)
  jobs/targets_manifest.csv         one row per fold target (designs + controls)
"""
import concurrent.futures as cf, json, os, subprocess, sys, time
import pandas as pd

OPT = "/home/hackathon/nick-work/optim"; OUT = f"{OPT}/out"; JOBS = f"{OPT}/jobs"
REPO = "/home/hackathon/hackathon-team-3"
ABG = "/home/hackathon/.venvs/abguard/bin/python"
CHK = f"{REPO}/NICCOLO/scripts/check_optimisation.py"
SPR = f"{REPO}/apheris/hackathon-data/her2-trastuzumab/trastuzumab_her2_spr_all_768.csv"
# HER2 domain IV construct (residues 489-607, native contiguous sequence incl. the 581-590 loop); identical to Jahn's calibration folds
HER2 = ("CHQLCARGHCWGPGPTQCVNCSQFLRGQECVEECRVLQGLPREYVNARHCLPCHPECQPQNGSVTCFGPEADQCVACAHYKDPPFCVARCPS"
        "GVKPDLSYMPIWKFPDEEGACQPCPIN")
WT = "WGGDGFYAMD"
N_FINAL, N_PROVIDED_CHECK, NPROC = 120, 170, 8
T0 = time.time(); log = lambda *a: print(f"[{time.time() - T0:6.1f}s]", *a, flush=True)
AAS = set("ACDEFGHIKLMNPQRSTVWY")

# ------------------------------------------------------------------ designs
sel = pd.read_csv(f"{OUT}/selected_diverse_200.csv")
spr = pd.read_csv(SPR)
wt = spr[spr.cdrh3 == WT].iloc[0]; VH, VL = wt.vh_sequence, wt.vl_sequence
i0 = VH.find(WT); assert VH.count(WT) == 1 and i0 == 98, "unexpected trastuzumab VH layout"
assert len(VH) == 120 and len(VL) == 107 and len(HER2) == 119
vh_of = lambda w: VH[:i0] + w + VH[i0 + 10:]
sel["sequence_id"] = [f"opt{r:03d}_{w}" for r, w in zip(sel["rank"], sel.window)]
des = pd.DataFrame({
    "sequence_id": sel.sequence_id, "kind": "paired", "heavy_sequence": [vh_of(w) for w in sel.window], "light_sequence": VL,
    "design_description": [f"trastuzumab VH/VL, CDRH3 window {WT}->{w} ({m}); prefilter p_binder={p:.2f}, pred_log10_kd={k:.2f} (a ranking, not a KD)"
                           for w, m, p, k in zip(sel.window, sel.mutations, sel.p_binder, sel.pred_log10_kd)]})
des.to_csv(f"{OUT}/designs_candidates.csv", index=False)
log(f"{len(des)} candidate designs assembled (VH {len(VH)} aa, VL {len(VL)} aa)")

def run_check(csv, out, provided=False):
    cmd = [ABG, CHK, csv, "-o", out] + (["--provided", SPR] if provided else [])
    r = subprocess.run(cmd, capture_output=True, text=True)
    if not os.path.exists(out): raise SystemExit(f"checker failed:\n{r.stdout[-800:]}\n{r.stderr[-800:]}")
    return pd.read_csv(out)

# ------------------------------------------------------------------ official check 1 (fast): hard checks + CDRH3 reference
t = time.time(); c1 = run_check(f"{OUT}/designs_candidates.csv", f"{OUT}/check_fast.csv")
log(f"check_optimisation (CDRH3 reference): {int(c1.passes.sum())}/{len(c1)} pass in {time.time() - t:.0f}s")
if (~c1.passes).any():
    print(c1[~c1.passes].failures.str.slice(0, 110).value_counts().head(5).to_string())
ok1 = c1[c1.passes].sequence_id.tolist()[:N_PROVIDED_CHECK]

# ------------------------------------------------------------------ official check 2: --provided (whole VH+VL vs every provided antibody)
# pure-python Levenshtein vs ~700 provided antibodies per design is slow, so shard it
sub = des[des.sequence_id.isin(ok1)].reset_index(drop=True)
shards = [sub.iloc[i::NPROC] for i in range(NPROC)]
for i, s in enumerate(shards): s.to_csv(f"{OUT}/_shard{i}.csv", index=False)
t = time.time()
with cf.ThreadPoolExecutor(NPROC) as ex:
    parts = list(ex.map(lambda i: run_check(f"{OUT}/_shard{i}.csv", f"{OUT}/check_provided_{i}.csv", provided=True), range(NPROC)))
c2 = pd.concat(parts, ignore_index=True)
log(f"check_optimisation --provided: {int(c2.passes.sum())}/{len(c2)} pass in {time.time() - t:.0f}s")
if (~c2.passes).any(): print(c2[~c2.passes].failures.str.slice(0, 110).value_counts().head(5).to_string())
for i in range(NPROC): os.remove(f"{OUT}/_shard{i}.csv")
ok2 = set(c2[c2.passes].sequence_id)

final = des[des.sequence_id.isin(ok2)].head(N_FINAL).reset_index(drop=True)
final.to_csv(f"{OUT}/final_variants.csv", index=False)
chk = c2.set_index("sequence_id")[["identity_to_trastuzumab", "cdrh3", "cdrh3_min_distance", "cdrh3_identity_to_trastuzumab", "passes"]]
chk.columns = ["identity_to_trastuzumab", "imgt_cdrh3", "cdrh3_min_distance_to_reference", "cdrh3_identity_to_trastuzumab", "official_check_passes"]
fs = sel.set_index("sequence_id").loc[final.sequence_id].join(chk).reset_index()
fs.to_csv(f"{OUT}/final_variants_scored.csv", index=False)
log(f"final: {len(final)} designs pass the official optimisation checks (both modes)")

# ------------------------------------------------------------------ JapanFold payloads (prepared, NOT submitted)
spr["k"] = [sum(a != b for a, b in zip(s, WT)) if len(s) == 10 else -1 for s in spr.cdrh3]
by = spr[spr.k >= 4].groupby("cdrh3").agg(res=("binding_result", lambda x: set(x)), kd=("kd_nM", "median"), vh=("vh_sequence", "first"), vl=("vl_sequence", "first")).reset_index()
bind = by[by.res.map(lambda r: r == {"quantified_binding"})].sort_values("kd")
non = by[by.res.map(lambda r: r == {"no_reportable_kd"})].sample(3, random_state=0)
ctl = [("ctl_wt", WT, VH, VL, "trastuzumab wild type", wt.kd_nM)]
for tag, df in (("ctl_binder", bind.iloc[[0, len(bind) // 2, -1]]), ("ctl_nonbinder", non)):
    for j, r in enumerate(df.itertuples(), 1): ctl.append((f"{tag}{j}_{r.cdrh3}", r.cdrh3, r.vh, r.vl, tag.replace("ctl_", "SPR-measured "), getattr(r, "kd", None)))
rows = []
for r in fs.itertuples():
    rows.append(dict(name=f"v{r.rank:03d}_{r.window}", role="design", window=r.window, rank=r.rank, vh=vh_of(r.window), vl=VL,
                     p_binder=round(r.p_binder, 4), pred_log10_kd=round(r.pred_log10_kd, 4), note=r.mutations))
for n, w, vh, vl, note, kd in ctl:
    rows.append(dict(name=n, role="control", window=w, rank=None, vh=vh, vl=vl, p_binder=None, pred_log10_kd=None, note=f"{note}; measured KD nM={kd}"))
man = pd.DataFrame(rows)
for r in man.itertuples():
    assert set(r.vh + r.vl) <= AAS and len(r.vh) == 120 and len(r.vl) == 107, r.name
targets = [{"name": r.name, "content": f">H\n{r.vh}\n>L\n{r.vl}\n>A\n{HER2}"} for r in man.itertuples()]
assert len({t["name"] for t in targets}) == len(targets)
man.drop(columns=["vh", "vl"]).assign(total_residues=len(VH) + len(VL) + len(HER2)).to_csv(f"{JOBS}/targets_manifest.csv", index=False)
for model, params in (("rf3", {"write_pae": True}), ("boltz2", None)):
    p = {"model": model, "name": f"optim-variants-{model}", "targets": targets}
    if params: p["params"] = params
    json.dump(p, open(f"{JOBS}/{model}_payload.json", "w"), indent=1)
log(f"payloads written: {len(targets)} targets each ({int((man.role == 'design').sum())} designs + {int((man.role == 'control').sum())} controls), {len(VH) + len(VL) + len(HER2)} residues per complex")
