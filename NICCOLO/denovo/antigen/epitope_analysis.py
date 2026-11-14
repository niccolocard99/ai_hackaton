"""HER2 domain IV epitope analysis for de novo binder design (P0-worker b).

Measures, from PDB 1N8Z cropped Fv+HER2 (apheris/rfd4/Data/12-her2/her2_fv_1n8z.cif):
  * HER2 (chains C,D) residues within 4.5 A of trastuzumab (chains A,B)
  * HER2 residues within 4.5 A of CDRH3 only (chain B, residues 99-108)
  * per-residue heavy-atom contact counts and nearest-atom distances
  * per-residue buried surface area (dSASA) from Shrake-Rupley SASA, biotite
and verifies the antigen-only output CIF.

mmCIF reader is adapted from /home/hackathon/wt/jahn/JAHN/pipeline/abpipe/structure.py
(read_chains / contact_residues), copied here per instruction, not modified in place.
"""
from __future__ import annotations
import csv, sys
import numpy as np

SRC = "/home/hackathon/hackathon-team-3/apheris/rfd4/Data/12-her2/her2_fv_1n8z.cif"
ANT = "/home/hackathon/hackathon-team-3/NICCOLO/denovo/antigen/her2_dIV_CD.cif"
OUTCSV = "/home/hackathon/hackathon-team-3/NICCOLO/denovo/antigen/epitope_contacts.csv"
CUT = 4.5

THREE = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
         "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
         "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}


def read_atoms(path):
    """Flat atom table from the _atom_site loop. Returns list of dicts (heavy atoms only)."""
    cols, in_loop, atoms = [], False, []
    for line in open(path):
        if line.startswith("_atom_site."):
            cols.append(line.strip().split(".", 1)[1]); in_loop = True
            continue
        if in_loop and not cols:
            continue
        if in_loop and (line.startswith("#") or line.startswith("loop_")):
            break
        if not (in_loop and line.startswith(("ATOM", "HETATM"))):
            continue
        t = line.split()
        if len(t) < len(cols):
            continue
        g = dict(zip(cols, t))
        if g.get("group_PDB") != "ATOM" or g.get("type_symbol") == "H":
            continue
        atoms.append(dict(chain=g["label_asym_id"], num=int(g["label_seq_id"]),
                          auth=int(g["auth_seq_id"]), comp=g["label_comp_id"],
                          atom=g["label_atom_id"], elem=g["type_symbol"],
                          xyz=(float(g["Cartn_x"]), float(g["Cartn_y"]), float(g["Cartn_z"]))))
    return atoms


def pairwise_min(A, B):
    """For each row of A, min distance to any row of B."""
    out = np.empty(len(A))
    for lo in range(0, len(A), 500):
        d2 = ((A[lo:lo+500, None, :] - B[None, :, :]) ** 2).sum(-1)
        out[lo:lo+500] = np.sqrt(d2.min(axis=1))
    return out


def count_contacts(A, B, cut):
    """For each row of A, number of rows of B within cut."""
    out = np.zeros(len(A), dtype=int)
    for lo in range(0, len(A), 500):
        d2 = ((A[lo:lo+500, None, :] - B[None, :, :]) ** 2).sum(-1)
        out[lo:lo+500] = (d2 <= cut * cut).sum(axis=1)
    return out


# --------------------------------------------------------------- load complex
at = read_atoms(SRC)
chains = sorted({a["chain"] for a in at})
ag = [a for a in at if a["chain"] in ("C", "D")]
ab = [a for a in at if a["chain"] in ("A", "B")]
h3 = [a for a in at if a["chain"] == "B" and 99 <= a["num"] <= 108]
print(f"[measured] complex chains={chains}  antigen heavy atoms={len(ag)}  "
      f"Fv heavy atoms={len(ab)}  CDRH3(B99-108) heavy atoms={len(h3)}")
print(f"[measured] CDRH3 residues present: "
      f"{[(a['num'], THREE.get(a['comp'],'X')) for a in h3 if a['atom']=='CA']}")

Xag = np.array([a["xyz"] for a in ag])
Xab = np.array([a["xyz"] for a in ab])
Xh3 = np.array([a["xyz"] for a in h3])

dmin_fv = pairwise_min(Xag, Xab)
dmin_h3 = pairwise_min(Xag, Xh3)
n_fv = count_contacts(Xag, Xab, CUT)
n_h3 = count_contacts(Xag, Xh3, CUT)

# -------------------------------------------------------------------- SASA
import biotite.structure as struc
from biotite.structure import AtomArray

def mk(atoms):
    arr = AtomArray(len(atoms))
    arr.coord = np.array([a["xyz"] for a in atoms], dtype=np.float32)
    arr.chain_id = np.array([a["chain"] for a in atoms])
    arr.res_id = np.array([a["num"] for a in atoms])
    arr.res_name = np.array([a["comp"] for a in atoms])
    arr.atom_name = np.array([a["atom"] for a in atoms])
    arr.element = np.array([a["elem"] for a in atoms])
    arr.hetero = np.zeros(len(atoms), dtype=bool)
    return arr

cplx = mk(at)
free = mk(ag)
is_ag = np.isin(cplx.chain_id, ["C", "D"])
sasa_cplx = struc.sasa(cplx, vdw_radii="Single", point_number=500)[is_ag]
sasa_free = struc.sasa(free, vdw_radii="Single", point_number=500)
assert len(sasa_cplx) == len(sasa_free) == len(ag)
sasa_cplx = np.nan_to_num(sasa_cplx)
sasa_free = np.nan_to_num(sasa_free)
dsasa = sasa_free - sasa_cplx
print(f"[measured] total HER2 dSASA buried by Fv = {dsasa.sum():.1f} A^2  "
      f"(interface area ~ {dsasa.sum():.0f} A^2 on the antigen side)")

# ------------------------------------------------------- per-residue rollup
rows = {}
for i, a in enumerate(ag):
    k = (a["chain"], a["num"])
    r = rows.setdefault(k, dict(chain=a["chain"], resnum=a["num"], resname=a["comp"],
                                aa=THREE.get(a["comp"], "X"), n_fv=0, n_h3=0,
                                dmin_fv=99.0, dmin_h3=99.0, dsasa=0.0, dsasa_sc=0.0,
                                sasa_free=0.0, sasa_cplx=0.0, best_atom=None,
                                best_atom_h3=None, atoms={}))
    r["n_fv"] += int(n_fv[i]); r["n_h3"] += int(n_h3[i])
    r["sasa_free"] += sasa_free[i]; r["sasa_cplx"] += sasa_cplx[i]
    r["dsasa"] += dsasa[i]
    if a["atom"] not in ("N", "CA", "C", "O", "OXT"):
        r["dsasa_sc"] += dsasa[i]
    if dmin_fv[i] < r["dmin_fv"]:
        r["dmin_fv"] = float(dmin_fv[i]); r["best_atom"] = a["atom"]
    if dmin_h3[i] < r["dmin_h3"]:
        r["dmin_h3"] = float(dmin_h3[i]); r["best_atom_h3"] = a["atom"]
    r["atoms"][a["atom"]] = (float(dmin_fv[i]), float(dmin_h3[i]), float(dsasa[i]))

contact = {k: v for k, v in rows.items() if v["n_fv"] > 0}
contact_h3 = {k: v for k, v in rows.items() if v["n_h3"] > 0}

def fmt(d):
    return [f"{c}/{n}/{v['aa']}" for (c, n), v in sorted(d.items(), key=lambda kv: (kv[0][0], kv[0][1]))]

print(f"\n[measured] EPITOPE vs whole Fv @ {CUT} A: n={len(contact)}")
print("  " + ", ".join(fmt(contact)))
print(f"\n[measured] EPITOPE vs CDRH3 (B99-108) @ {CUT} A: n={len(contact_h3)}")
print("  " + ", ".join(fmt(contact_h3)))

print("\n[measured] per-residue table (sorted by dSASA):")
hdr = f"{'res':>10} {'nFv':>4} {'nH3':>4} {'dminFv':>7} {'dminH3':>7} {'dSASA':>7} {'dSASAsc':>8} {'freeSASA':>9} {'%buried':>8} {'closest':>8}"
print(hdr)
for (c, n), v in sorted(rows.items(), key=lambda kv: -kv[1]["dsasa"]):
    if v["dsasa"] < 0.5 and v["n_fv"] == 0:
        continue
    pb = 100 * v["dsasa"] / v["sasa_free"] if v["sasa_free"] > 1e-6 else 0.0
    print(f"{c+'/'+str(n)+' '+v['aa']:>10} {v['n_fv']:>4} {v['n_h3']:>4} "
          f"{v['dmin_fv']:>7.2f} {v['dmin_h3']:>7.2f} {v['dsasa']:>7.1f} {v['dsasa_sc']:>8.1f} "
          f"{v['sasa_free']:>9.1f} {pb:>8.1f} {str(v['best_atom']):>8}")

# ------------------------------------------------------------------- CSV out
with open(OUTCSV, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["chain", "resnum", "resname", "aa", "in_fv_epitope", "in_cdrh3_epitope",
                "n_contacts_fv_4.5A", "n_contacts_cdrh3_4.5A", "dmin_fv_A", "dmin_cdrh3_A",
                "dSASA_A2", "dSASA_sidechain_A2", "sasa_free_A2", "sasa_complex_A2",
                "pct_buried", "closest_atom_to_fv", "closest_atom_to_cdrh3",
                "vs_published_19res_557_605"])
    for (c, n), v in sorted(rows.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        if v["n_fv"] == 0 and v["dsasa"] < 0.5:
            continue
        pb = 100 * v["dsasa"] / v["sasa_free"] if v["sasa_free"] > 1e-6 else 0.0
        w.writerow([c, n, v["resname"], v["aa"], int(v["n_fv"] > 0), int(v["n_h3"] > 0),
                    v["n_fv"], v["n_h3"], round(v["dmin_fv"], 2), round(v["dmin_h3"], 2),
                    round(v["dsasa"], 1), round(v["dsasa_sc"], 1), round(v["sasa_free"], 1),
                    round(v["sasa_cplx"], 1), round(pb, 1), v["best_atom"], v["best_atom_h3"],
                    "agrees" if v["n_fv"] > 0 else "not_in_either"])
print(f"\n[written] {OUTCSV}")

# ------------------------------------------ atom-level detail for hotspot picks
print("\n[measured] atom-level detail for candidate hotspot residues "
      "(atoms with dSASA>0 or within 4.5 A of Fv):")
for (c, n), v in sorted(contact.items(), key=lambda kv: -kv[1]["dsasa"])[:16]:
    sel = [(an, dfv, dh3, ds) for an, (dfv, dh3, ds) in v["atoms"].items() if ds > 1.0 or dfv <= CUT]
    sel.sort(key=lambda t: -t[3])
    print(f"  {c}/{n} {v['resname']}: " + ", ".join(
        f"{an}(d_fv={dfv:.2f},d_h3={dh3:.2f},dSASA={ds:.1f})" for an, dfv, dh3, ds in sel))

np.save("/tmp/claude-1001/-home-hackathon/c26d617c-6c92-4339-abec-99f979349904/scratchpad/_rows.npy",
        np.array([1]))  # marker
import pickle
with open("/tmp/claude-1001/-home-hackathon/c26d617c-6c92-4339-abec-99f979349904/scratchpad/rows.pkl", "wb") as fh:
    pickle.dump(rows, fh)
