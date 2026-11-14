"""Is the trastuzumab epitope reachable by one VHH paratope, and where is its core?

Reads the cropped 1N8Z complex, takes the measured epitope (HER2 residues within
4.5 A of trastuzumab), and reports:
  - how many trastuzumab residues each epitope residue touches (contact weight)
  - the pairwise spread of the epitope, to see whether one small paratope can cover it
  - the best-covered patch for a given paratope radius

A VHH buries roughly 600-800 A^2 and its paratope spans about 18-22 A, so a patch
whose residues all sit inside a ~20 A sphere is the realistic target.

Usage: python3 epitope_patch.py [path/to/her2_fv_1n8z.cif]   (default: the organiser's copy in this repo)
"""
import math
import sys
from collections import defaultdict
from pathlib import Path


def parse_atom_site(path):
    lines = open(path, encoding="utf8", errors="replace").read().split("\n")
    i = 0
    while i < len(lines):
        if lines[i].strip() == "loop_":
            j = i + 1
            names = []
            while j < len(lines) and lines[j].lstrip().startswith("_"):
                names.append(lines[j].strip())
                j += 1
            if names and names[0].startswith("_atom_site."):
                cols = [n.split(".", 1)[1] for n in names]
                rows = []
                while j < len(lines):
                    s = lines[j].strip()
                    if s in ("#", "") or s.startswith("_") or s == "loop_":
                        break
                    rows.append(s.split())
                    j += 1
                return cols, rows
            i = j
        else:
            i += 1
    raise SystemExit("no _atom_site loop found")


def main():
    default = Path(__file__).resolve().parents[3] / "apheris" / "rfd4" / "Data" / "12-her2" / "her2_fv_1n8z.cif"
    path = sys.argv[1] if len(sys.argv) > 1 else str(default)
    cols, rows = parse_atom_site(path)
    idx = {c: k for k, c in enumerate(cols)}
    atoms = []
    for r in rows:
        try:
            atoms.append(
                dict(
                    ch=r[idx["label_asym_id"]],
                    seq=int(r[idx["label_seq_id"]]),
                    comp=r[idx["label_comp_id"]],
                    name=r[idx["label_atom_id"]],
                    x=float(r[idx["Cartn_x"]]),
                    y=float(r[idx["Cartn_y"]]),
                    z=float(r[idx["Cartn_z"]]),
                )
            )
        except (ValueError, KeyError):
            continue

    ab = [a for a in atoms if a["ch"] in ("A", "B")]  # trastuzumab Fv
    ag = [a for a in atoms if a["ch"] in ("C", "D")]  # HER2 domain IV

    CUT = 4.5
    # contact weight = number of distinct antibody residues an epitope residue touches
    weight = defaultdict(set)
    close_atoms = defaultdict(set)
    for a in ag:
        for b in ab:
            dx, dy, dz = a["x"] - b["x"], a["y"] - b["y"], a["z"] - b["z"]
            if dx * dx + dy * dy + dz * dz <= CUT * CUT:
                key = (a["ch"], a["seq"], a["comp"])
                weight[key].add((b["ch"], b["seq"], b["comp"]))
                close_atoms[key].add(a["name"])

    epitope = sorted(weight, key=lambda k: -len(weight[k]))
    print(f"=== measured epitope: {len(epitope)} HER2 residues within {CUT} A of trastuzumab ===")
    print(f"{'residue':>12}  {'n_ab_res':>8}  {'heavy/light':>12}  contact atoms")
    for k in epitope:
        partners = weight[k]
        h = sum(1 for p in partners if p[0] == "B")
        l = sum(1 for p in partners if p[0] == "A")
        print(
            f"{k[2]}{k[1]} ({k[0]}){'':>2}  {len(partners):>8}  "
            f"{'H' + str(h) + '/L' + str(l):>12}  {','.join(sorted(close_atoms[k]))}"
        )

    # centroid of each epitope residue's contacting atoms
    cen = {}
    for k in epitope:
        sel = [a for a in ag if (a["ch"], a["seq"], a["comp"]) == k]
        cen[k] = (
            sum(a["x"] for a in sel) / len(sel),
            sum(a["y"] for a in sel) / len(sel),
            sum(a["z"] for a in sel) / len(sel),
        )

    def dist(p, q):
        return math.dist(cen[p], cen[q])

    span = max(dist(p, q) for p in epitope for q in epitope)
    print(f"\n=== epitope spread ===")
    print(f"max residue-centroid separation across the whole epitope: {span:.1f} A")
    print("(a VHH paratope spans roughly 18-22 A, a Fab footprint much more)")

    # best patch: for each epitope residue as centre, which others fall inside R
    print(f"\n=== best patch per paratope radius (sum of contact weight covered) ===")
    for R in (9.0, 11.0, 13.0, 15.0):
        best = None
        for c in epitope:
            inside = [k for k in epitope if dist(c, k) <= R]
            score = sum(len(weight[k]) for k in inside)
            if best is None or score > best[0]:
                best = (score, c, inside)
        score, c, inside = best
        tot = sum(len(weight[k]) for k in epitope)
        hw = sum(sum(1 for p in weight[k] if p[0] == 'B') for k in inside)
        lw = sum(sum(1 for p in weight[k] if p[0] == 'A') for k in inside)
        htot = sum(sum(1 for p in weight[k] if p[0] == 'B') for k in epitope)
        ltot = sum(sum(1 for p in weight[k] if p[0] == 'A') for k in epitope)
        names = ", ".join(f"{k[2]}{k[1]}" for k in sorted(inside, key=lambda t: t[1]))
        print(f"\n  R={R:>4.1f} A around {c[2]}{c[1]}: covers {len(inside)}/{len(epitope)} residues, "
              f"{score}/{tot} ({100*score/tot:.0f}%) of contact weight")
        print(f"        {names}")
        print(f"        heavy-chain contacts {hw}/{htot} ({100*hw/htot:.0f}%), light-chain {lw}/{ltot} ({100*lw/ltot:.0f}%)")


main()
