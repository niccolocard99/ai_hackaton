#!/usr/bin/env python3
"""Per-design epitope check for generated VHH-HER2 complexes. Needs gemmi (pip install gemmi).

The design model is only *asked* to contact the hotspots. This measures what it did, for each generated
complex, against the 1N8Z structure of trastuzumab bound to HER2:

  frame_rmsd       HER2 CA RMSD to the 1N8Z crop before superposition. The model does NOT keep the input
                   frame (about 150 A on the smoke designs), so everything below superposes HER2 onto
                   1N8Z first; frame_rmsd_after_superposition (about 0.05 A) shows HER2 stayed rigid
  hs_<atom>        distance (A) from each hotspot atom to the closest VHH heavy atom
  n_hotspots_4p5   how many of the four hotspot atoms have a VHH atom within 4.5 A
  n_contacts       HER2 residues within 4.5 A of the VHH
  epitope_cov      share of trastuzumab's 19 epitope residues the VHH touches
  contacts_in_ep   share of the VHH's HER2 contacts that lie inside trastuzumab's epitope
  fv_min_dist      closest VHH heavy atom to trastuzumab's Fv (chains A, B of 1N8Z), after superposing HER2
  fv_pairs_3p0     VHH/Fv heavy-atom pairs closer than 3.0 A. Steric overlap is a proxy for competition
                   with trastuzumab, not a measurement of it: Adaptyv's assay decides

    ./epitope_check.py RUN_OR_CIF_OR_DIR [...] [-o epitope_check.csv]

A RUN is a folder under runs/ (as run_design.sh names them); a directory or a *_generated.cif path also
works. Output chains are C (HER2 489-580, numbered 1-92), D (591-607, numbered 1-17) and H (the VHH).
"""
import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import gemmi

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
HER2_CIF = REPO / "apheris" / "rfd4" / "Data" / "12-her2" / "her2_fv_1n8z.cif"
TARGET_CIF = HERE / "spec" / "target_her2_d4.cif"
CUT = 4.5
OFFSET = {"C": 488, "D": 590}  # output residue number + offset = 1N8Z (mature HER2) numbering

# hotspots as (chain, output residue number, residue name, atom); 1N8Z numbering in the label
HOTSPOTS = (("D", 3, "LYS", "NZ", "D593"), ("C", 85, "PHE", "CD1", "C573"),
            ("D", 1, "ILE", "CG2", "D591"), ("C", 84, "PRO", "CD", "C572"))


def load(path):
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    return st[0]


def heavy(chain, num_shift=0):
    return [(res.seqid.num + num_shift, res.name, a.name, a.pos)
            for res in chain if res.name != "HOH" for a in res if a.element.name != "H"]


def grid(atoms, cell):
    g = defaultdict(list)
    for i, (_, _, _, p) in enumerate(atoms):
        g[(int(p.x // cell), int(p.y // cell), int(p.z // cell))].append(i)
    return g


def near(a_atoms, b_atoms, cut=CUT):
    """Pairs (a_idx, b_idx, dist) with dist <= cut."""
    g = grid(b_atoms, cut)
    out = []
    for i, (_, _, _, p) in enumerate(a_atoms):
        cx, cy, cz = int(p.x // cut), int(p.y // cut), int(p.z // cut)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for j in g.get((cx + dx, cy + dy, cz + dz), ()):
                        d = p.dist(b_atoms[j][3])
                        if d <= cut:
                            out.append((i, j, d))
    return out


def vec(v):
    return gemmi.Position(v.x, v.y, v.z)


def ca_positions(model, chains):
    out = []
    for c in chains:
        for res in model[c]:
            a = res.find_atom("CA", "*")
            out.append(a.pos if a else None)
    return out


def trastuzumab_epitope(m1):
    """HER2 residues (1N8Z numbering) within 4.5 A of the Fv, same definition as epitope_patch.py."""
    fv = heavy(m1["A"]) + heavy(m1["B"])
    her2 = heavy(m1["C"]) + heavy(m1["D"])
    return {her2[j][0] for _, j, _ in near(fv, her2)}


def collect(args):
    items = []
    for a in args:
        p = Path(a)
        if p.is_file():
            items.append((p.parent.name, p))
        elif p.is_dir():
            items += [(p.name, f) for f in sorted(p.rglob("*_generated.cif"))]
        else:
            d = HERE / "runs" / a / "design" / "generation" / "inference_samples"
            if not d.is_dir():
                sys.exit(f"no such run, file or directory: {a}")
            items += [(a, f) for f in sorted(d.glob("*/*_generated.cif"))]
    return items


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("-o", "--out", default="epitope_check.csv")
    ap.add_argument("--her2-cif", default=str(HER2_CIF), help="1N8Z crop with trastuzumab (chains A, B, C, D)")
    a = ap.parse_args()

    m1 = load(a.her2_cif)
    epitope = trastuzumab_epitope(m1)
    fv = heavy(m1["A"]) + heavy(m1["B"])
    ref_ca = [p for p in ca_positions(m1, ("C", "D")) if p is not None]
    tgt = load(TARGET_CIF)
    tgt_ca = [p for p in ca_positions(tgt, ("C", "D")) if p is not None]
    if len(ref_ca) != len(tgt_ca):
        sys.exit("input target and 1N8Z crop disagree on the number of HER2 residues")

    rows = []
    for run, path in collect(a.inputs):
        m = load(path)
        sid = f"{run}__{path.stem.replace('_generated', '')}"
        gen_ca = [p for p in ca_positions(m, ("C", "D")) if p is not None]
        if len(gen_ca) != len(ref_ca):
            print(f"skip {sid}: HER2 has {len(gen_ca)} residues, expected {len(ref_ca)}", file=sys.stderr)
            continue
        raw = math.sqrt(sum(p.dist(q) ** 2 for p, q in zip(gen_ca, ref_ca)) / len(ref_ca))
        sup = gemmi.superpose_positions(ref_ca, gen_ca)  # moves the generated HER2 onto 1N8Z
        moved = [vec(sup.transform.apply(p)) for p in gen_ca]
        check_rmsd = math.sqrt(sum(p.dist(q) ** 2 for p, q in zip(moved, ref_ca)) / len(ref_ca))
        if abs(check_rmsd - sup.rmsd) > 0.5:  # guards the direction of the transform: every number below depends on it
            sys.exit(f"{sid}: superposition applied the wrong way round ({check_rmsd:.2f} A vs {sup.rmsd:.2f} A)")
        # HER2 residue keys in 1N8Z numbering: output numbering restarts at 1 per chain
        her2 = heavy(m["C"], OFFSET["C"]) + heavy(m["D"], OFFSET["D"])
        vhh = heavy(m["H"])

        row = {"sequence_id": sid, "frame_rmsd": round(raw, 2), "frame_rmsd_after_superposition": round(sup.rmsd, 2)}
        n_hot = 0
        for ch, num, rn, an, label in HOTSPOTS:
            res = next((r for r in m[ch] if r.seqid.num == num), None)
            atom = res.find_atom(an, "*") if res is not None else None
            if res is None or atom is None or res.name != rn:
                row[f"hs_{label}"] = ""
                continue
            d = min(atom.pos.dist(p) for (_, _, _, p) in vhh)
            row[f"hs_{label}"] = round(d, 2)
            n_hot += d <= CUT
        row["n_hotspots_4p5"] = n_hot

        contacts = {her2[i][0] for i, _, _ in near(her2, vhh)}
        row["n_contacts"] = len(contacts)
        row["epitope_cov"] = round(len(contacts & epitope) / len(epitope), 2)
        row["contacts_in_ep"] = round(len(contacts & epitope) / len(contacts), 2) if contacts else 0.0

        vhh_ref = [(n, rn, an, vec(sup.transform.apply(p))) for (n, rn, an, p) in vhh]
        pairs = near(vhh_ref, fv, cut=CUT)
        row["fv_min_dist"] = round(min((d for _, _, d in pairs), default=float("nan")), 2) if pairs else ">4.5"
        row["fv_pairs_3p0"] = sum(1 for _, _, d in pairs if d < 3.0)
        rows.append(row)

    if not rows:
        sys.exit("no designs found")
    with open(a.out, "w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} designs -> {a.out}; trastuzumab epitope = {len(epitope)} HER2 residues")
    keys = ["frame_rmsd", "n_hotspots_4p5", "n_contacts", "epitope_cov", "contacts_in_ep", "fv_min_dist", "fv_pairs_3p0"]
    print("  ".join(f"{k}" for k in ["sequence_id"] + keys))
    for r in rows:
        print("  ".join(str(r[k]) for k in ["sequence_id"] + keys))


if __name__ == "__main__":
    main()
