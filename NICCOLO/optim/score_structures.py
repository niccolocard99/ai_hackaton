#!/usr/bin/env python3
"""Score the RF3 folds of the optimisation variants against the 1N8Z crystal pose: epitope placement, Fv pose, interface PAE.

    /home/hackathon/.venvs/abmodels/bin/python -I score_structures.py STRUCT_DIR RESULTS_JSON MANIFEST_CSV VARIANTS_SCORED_CSV REF_CIF OUT_CSV

STRUCT_DIR is the unzipped JapanFold archive (an untrusted download: parsed as data only). Reference = her2_fv_1n8z.cif
(chains A=VL, B=VH, C+D = HER2 domain IV residues 489-580 and 591-607; the 581-590 loop is not resolved in the crystal).
Folds have chains H (VH, 120), L (VL, 107), A (HER2 489-607, 119) in that file order; CDRH3 window = VH residues 99-108.

Metrics (all per fold):
  n_her2_contact / ep_recall / ep_precision / n_offep   HER2 residues within 4.5 A (heavy atoms) of the Fv, vs the crystal epitope
  cdrh3_*                                              the same restricted to the 10-residue CDRH3 window
  lrmsd_xtal   Fv CA RMSD (CDRH3 window excluded) to the crystal after fitting HER2 CAs   (CAPRI-style ligand RMSD)
  drift_wt     the same against the RF3 fold of wild-type trastuzumab
  pae_fv_her2  mean PAE over all Fv<->HER2 token pairs; ipae = only the predicted contact residues; pae_cdrh3_ep = CDRH3 window <-> crystal epitope
"""
import json, sys
import numpy as np, pandas as pd
from Bio.PDB import MMCIFParser
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

STRUCT, RES, MAN, VAR, REF, OUT = sys.argv[1:7]
CUT, OFF = 4.5, 488                  # contact cutoff (A, as in LEO's interface_1n8z.py); fold HER2 residue k <-> crystal residue k+488
WIN = set(range(99, 109))            # CDRH3 window, VH residue numbers
LOOP = set(range(581, 591))          # HER2 loop not resolved in the crystal
T = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G", HIS="H", ILE="I", LEU="L", LYS="K", MET="M",
         PHE="F", PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V")
pd.set_option("display.width", 230); pd.set_option("display.max_columns", 40)


def parse(path):
    """chain id -> dict(n residue numbers, aa, xyz heavy atoms, an residue number per atom, ca {resnum: xyz}); dict order = file order"""
    out = {}
    for ch in MMCIFParser(QUIET=True).get_structure("x", path)[0]:
        n, aa, xyz, an, ca = [], [], [], [], {}
        for r in ch:
            if r.id[0] != " ": continue
            n.append(r.id[1]); aa.append(T.get(r.get_resname(), "X"))
            for a in r:
                if a.element != "H": xyz.append(a.coord); an.append(r.id[1])
            if "CA" in r: ca[r.id[1]] = np.array(r["CA"].coord, float)
        out[ch.id] = dict(n=np.array(n), aa="".join(aa), xyz=np.array(xyz, float), an=np.array(an), ca=ca)
    return out


def shift(d, k):
    return dict(d, n=d["n"] + k, an=d["an"] + k, ca={r + k: v for r, v in d["ca"].items()})


def cat(a, b):
    return dict(n=np.concatenate([a["n"], b["n"]]), aa=a["aa"] + b["aa"], xyz=np.vstack([a["xyz"], b["xyz"]]),
                an=np.concatenate([a["an"], b["an"]]), ca={**a["ca"], **b["ca"]})


def contact(a, b, cut=CUT):
    """(residue numbers of a, residue numbers of b) with any heavy-atom pair closer than cut"""
    nb = cKDTree(b["xyz"]).query_ball_point(a["xyz"], cut)
    ia = [i for i, js in enumerate(nb) if js]
    ib = sorted({j for js in nb for j in js})
    return set(a["an"][ia].tolist()), set(b["an"][ib].tolist())


def sub(d, nums):
    m = np.isin(d["an"], list(nums)); return dict(xyz=d["xyz"][m], an=d["an"][m])


def kabsch(P, Q):
    """R, t with P @ R.T + t ~ Q (least squares)"""
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    R = Vt.T @ np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))]) @ U.T
    return R, qc - R @ pc


def pose(m, ref):
    """fit m's HER2 CAs onto ref's, return (HER2 fit RMSD, Fv CA RMSD without the CDRH3 window)"""
    com = sorted(set(m["A"]["ca"]) & set(ref["A"]["ca"]))
    P = np.array([m["A"]["ca"][k] for k in com]); Q = np.array([ref["A"]["ca"][k] for k in com])
    R, t = kabsch(P, Q)
    fit = float(np.sqrt(((P @ R.T + t - Q) ** 2).sum(1).mean()))
    d = [((m[c]["ca"][k] @ R.T + t) - ref[c]["ca"][k]) ** 2 for c in "HL" for k in sorted(set(m[c]["ca"]) & set(ref[c]["ca"]))
         if not (c == "H" and k in WIN)]
    return fit, float(np.sqrt(np.array(d).sum(1).mean()))


def model(path):
    """fold -> {'H','L','A'} with HER2 renumbered to crystal numbering; also returns the raw file-order chain list"""
    m = parse(path); order = list(m)
    assert order == ["H", "L", "A"] and [len(m[c]["n"]) for c in order] == [120, 107, 119], (path, order)
    m["A"] = shift(m["A"], OFF)
    return m, order


# ---------------------------------------------------------------- crystal reference
rf = parse(REF)
her2_x = cat(rf["C"], rf["D"])
xtal = {"H": rf["B"], "L": rf["A"], "A": her2_x}
resolved = set(her2_x["ca"])
ep_h, _ = contact(her2_x, rf["B"]); ep_l, _ = contact(her2_x, rf["A"]); EP = ep_h | ep_l
EPc, _ = contact(her2_x, sub(rf["B"], WIN))
print(f"crystal {REF.split('/')[-1]}: HER2 resolved {len(resolved)} residues; epitope {len(EP)} residues "
      f"({len(ep_h)} via VH, {len(ep_l)} via VL); CDRH3-window epitope {len(EPc)}; epitope residues: {sorted(EP)}")

# ---------------------------------------------------------------- folds
rows = {r["id"]: r for r in json.load(open(RES))["rows"]}
wt, _ = model(f"{STRUCT}/ctl_wt.cif")
assert wt["H"]["aa"] == xtal["H"]["aa"] and wt["L"]["aa"] == xtal["L"]["aa"], "crystal Fv sequence differs from the WT fold"
out, margins = [], []
for name, r in rows.items():
    if r["status"] != "ok": out.append(dict(name=name, status=r["status"])); continue
    m, order = model(f"{STRUCT}/{name}.cif")
    A, H, L = m["A"], m["H"], m["L"]
    ph, fh = contact(A, H); pl, fl = contact(A, L); P = ph | pl
    pc, wc = contact(A, sub(H, WIN))
    fit_x, l_x = pose(m, xtal); _, l_wt = pose(m, wt)
    # PAE: token order = file order (H, L, A). Check it against the alternative (A, H, L): the right one has low within-chain blocks
    Pm = np.load(f"{STRUCT}/{name}_pae.npz")["pae"].astype(float)
    off = {"H": 0, "L": 120, "A": 227}
    blk = lambda sl: np.mean([Pm[a:b, a:b].mean() for a, b in sl])
    margins.append(blk([(0, 119), (119, 239), (239, 346)]) - blk([(0, 120), (120, 227), (227, 346)]))
    tk = lambda c, ks, o=0: [off[c] + (k - o) - 1 for k in sorted(ks)]
    fv_all, her_all = tk("H", H["n"]) + tk("L", L["n"]), tk("A", A["n"], OFF)
    fvc, hc = tk("H", fh) + tk("L", fl), tk("A", P, OFF)
    w, e = tk("H", WIN), tk("A", EP, OFF)
    sym = lambda x, y: float((Pm[np.ix_(x, y)].mean() + Pm[np.ix_(y, x)].mean()) / 2) if x and y else np.nan
    pi = r["pair_chains_iptm"]          # sorted chain ids: 0=A (HER2), 1=H, 2=L
    out.append(dict(name=name, status="ok", iptm=r["iptm"], ptm=r["ptm"], plddt=r["plddt"], ranking_score=r["ranking_score"],
                    has_clash=r["has_clash"], iptm_AH=pi["0"]["1"], iptm_AL=pi["0"]["2"], runtime_s=r["runtime_s"],
                    n_her2_contact=len(P), n_loop581_590=len(P & LOOP), ep_recall=len(P & EP) / len(EP),
                    ep_precision=len(P & EP) / max(1, len(P & resolved)), n_offep=len((P & resolved) - EP),
                    cdrh3_n_res_contact=len(wc), cdrh3_n_her2_res=len(pc), cdrh3_ep_recall=len(pc & EPc) / len(EPc),
                    fit_her2_rmsd=fit_x, lrmsd_xtal=l_x, drift_wt=l_wt,
                    pae_fv_her2=sym(fv_all, her_all), ipae=sym(fvc, hc), pae_cdrh3_ep=sym(w, e)))
d = pd.DataFrame(out)
print(f"\nPAE token-order check: within-chain PAE of the wrong segmentation minus the right one, min over folds = {min(margins):.2f} A "
      f"(positive = file order H,L,A is right); mean PAE {Pm.mean():.1f}")

man = pd.read_csv(MAN); var = pd.read_csv(VAR)[["rank", "n_mut", "mutations"]]
d = man[["name", "role", "window", "rank", "p_binder", "pred_log10_kd", "note"]].merge(d, on="name", how="right").merge(var, on="rank", how="left")
d.to_csv(OUT, index=False); print(f"wrote {OUT} ({len(d)} rows)")

# ---------------------------------------------------------------- report
cols = ["name", "iptm", "iptm_AH", "n_her2_contact", "ep_recall", "ep_precision", "n_offep", "cdrh3_ep_recall", "lrmsd_xtal", "drift_wt",
        "ipae", "pae_cdrh3_ep"]
ctl = d[d.role == "control"]
print("\n=== controls (SPR-measured; ctl_wt = trastuzumab) ===")
print(ctl[cols].round(3).to_string(index=False))
des = d[d.role == "design"].copy()
q = lambda s: pd.Series({"min": s.min(), "p10": s.quantile(.1), "median": s.median(), "p90": s.quantile(.9), "max": s.max()})
print(f"\n=== {len(des)} designs ===")
print(des[["iptm", "iptm_AH", "n_her2_contact", "ep_recall", "ep_precision", "n_offep", "cdrh3_n_res_contact", "cdrh3_ep_recall", "lrmsd_xtal",
           "drift_wt", "pae_fv_her2", "ipae", "pae_cdrh3_ep", "plddt"]].apply(q).T.round(3).to_string())
w0 = ctl[ctl.name == "ctl_wt"].iloc[0]
print(f"\nvs the WT fold: ep_recall >= WT {int((des.ep_recall >= w0.ep_recall).sum())}/{len(des)}; ipae <= WT {int((des.ipae <= w0.ipae).sum())}/{len(des)}; "
      f"pae_cdrh3_ep <= WT {int((des.pae_cdrh3_ep <= w0.pae_cdrh3_ep).sum())}/{len(des)}; drift_wt < 3 A {int((des.drift_wt < 3).sum())}/{len(des)}, "
      f"< 5 A {int((des.drift_wt < 5).sum())}/{len(des)}; lrmsd_xtal < 5 A {int((des.lrmsd_xtal < 5).sum())}/{len(des)}")
print("\n=== Spearman among designs (classifier / mutation count vs fold metrics) ===")
sp = {}
for a in ("p_binder", "pred_log10_kd", "n_mut"):
    sp[a] = {b: spearmanr(des[a], des[b]).statistic for b in ("iptm", "iptm_AH", "ep_recall", "n_her2_contact", "drift_wt", "ipae", "pae_cdrh3_ep")}
print(pd.DataFrame(sp).T.round(2).to_string())
print("\n=== top 8 designs by iptm ===")
print(des.sort_values("iptm", ascending=False)[["name", "mutations"] + cols[1:]].head(8).round(3).to_string(index=False))
print("\n=== 8 designs with the largest pose drift from the WT fold ===")
print(des.sort_values("drift_wt", ascending=False)[["name", "mutations"] + cols[1:]].head(8).round(3).to_string(index=False))
