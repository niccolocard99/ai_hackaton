#!/usr/bin/env python3
"""GEOMETRIC chemistry classification of every interface contact pair.

    python contact_chemistry.py <file.cif> <antibody_chains> <antigen_chains> [cutoff]

THIS IS NOT AN ENERGY CALCULATION. No force field is evaluated, nothing is
minimised, no solvation term exists here. FoldX, Rosetta, PyRosetta and OpenMM
were all checked for and none is installed in this environment, and nothing may
be installed, so instead of fabricating a dG this script reports what can be
read straight off the coordinates: which atom pairs are close enough, and of the
right chemistry, to be the interaction type named.

Rules applied (all purely geometric, hydrogens absent from the models):
  salt bridge   formally charged group of opposite sign, heavy atoms < 4.0 A
                (Asp OD1/OD2, Glu OE1/OE2, C-term OXT  vs  Lys NZ,
                 Arg NE/NH1/NH2, His ND1/NE2, N-terminal N)
  H-bond        one donor heavy atom and one acceptor heavy atom, < 3.5 A.
                Donor/acceptor roles are taken from a residue/atom table; they
                CANNOT be confirmed without hydrogens, and no angle criterion is
                applied, so these are "compatible with an H-bond", not proven
                H-bonds.
  hydrophobic   carbon-carbon (or carbon-sulfur) pair < 4.5 A
  polar/vdW     anything else inside the cutoff

Each residue pair is reported with its closest approach, the strongest
interaction class it supports, and the atom pair that justifies it.
"""
import math
import sys
from collections import defaultdict

from bsa_1n8z import load_atoms

ANION = {
    "ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"},
}
CATION = {
    "LYS": {"NZ"}, "ARG": {"NE", "NH1", "NH2"}, "HIS": {"ND1", "NE2"},
}

# Heavy atoms that carry a hydrogen and can donate.
DONOR = {
    "ARG": {"N", "NE", "NH1", "NH2"}, "ASN": {"N", "ND2"},
    "GLN": {"N", "NE2"}, "HIS": {"N", "ND1", "NE2"}, "LYS": {"N", "NZ"},
    "SER": {"N", "OG"}, "THR": {"N", "OG1"}, "TRP": {"N", "NE1"},
    "TYR": {"N", "OH"}, "CYS": {"N", "SG"},
}
DONOR_DEFAULT = {"N"}          # backbone amide, every residue except proline
# Heavy atoms with a lone pair that can accept.
ACCEPTOR = {
    "ASP": {"O", "OD1", "OD2"}, "GLU": {"O", "OE1", "OE2"},
    "ASN": {"O", "OD1"}, "GLN": {"O", "OE1"}, "HIS": {"O", "ND1", "NE2"},
    "SER": {"O", "OG"}, "THR": {"O", "OG1"}, "TYR": {"O", "OH"},
    "MET": {"O", "SD"}, "CYS": {"O", "SG"},
}
ACCEPTOR_DEFAULT = {"O"}       # backbone carbonyl
TERMINAL_ACID = {"OXT"}

HB_CUT = 3.5
SB_CUT = 4.0
HPHOB_CUT = 4.5


def is_donor(resn, atom):
    if atom in TERMINAL_ACID:
        return False
    if atom == "N" and resn == "PRO":
        return False
    return atom in DONOR.get(resn, set()) | DONOR_DEFAULT


def is_acceptor(resn, atom):
    return atom in ACCEPTOR.get(resn, set()) | ACCEPTOR_DEFAULT | TERMINAL_ACID


def charge_of(resn, atom):
    if atom in TERMINAL_ACID or atom in ANION.get(resn, set()):
        return -1
    if atom in CATION.get(resn, set()):
        return +1
    return 0


def main():
    path, gA, gB = sys.argv[1], sys.argv[2], sys.argv[3]
    cutoff = float(sys.argv[4]) if len(sys.argv) > 4 else 4.5
    setA = set(gA.split(",")) if "," in gA else set(gA)
    setB = set(gB.split(",")) if "," in gB else set(gB)

    print("=== CONTACT CHEMISTRY - GEOMETRIC CLASSIFICATION, NOT AN ENERGY ===")
    print(f"  file {path}, cutoff {cutoff} A, heavy atoms, first model")
    print("  FoldX / Rosetta / PyRosetta / OpenMM: none installed, none used.")
    print(f"  thresholds: salt bridge <{SB_CUT} A between opposite formal "
          f"charges; H-bond-compatible <{HB_CUT} A donor-acceptor;")
    print(f"              hydrophobic <{HPHOB_CUT} A C/S-C/S; no angles, no "
          f"hydrogens, no desolvation.\n")

    atoms, _seqs = load_atoms(path)
    A = [a for a in atoms if a[0] in setA]
    B = [a for a in atoms if a[0] in setB]

    cell = cutoff
    grid = defaultdict(list)
    for j, b in enumerate(B):
        x, y, z = b[6]
        grid[(int(x // cell), int(y // cell), int(z // cell))].append(j)

    pairs = defaultdict(list)
    c2 = cutoff ** 2
    offs = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1)
            for dz in (-1, 0, 1)]
    for a in A:
        x, y, z = a[6]
        kx, ky, kz = int(x // cell), int(y // cell), int(z // cell)
        for dx, dy, dz in offs:
            for j in grid.get((kx + dx, ky + dy, kz + dz), ()):
                b = B[j]
                bx, by, bz = b[6]
                d2 = (x - bx) ** 2 + (y - by) ** 2 + (z - bz) ** 2
                if d2 > c2:
                    continue
                pairs[(a[0], a[1], a[3], b[0], b[1], b[3])].append(
                    (a[4], a[5], b[4], b[5], math.sqrt(d2)))

    rank = {"SALT BRIDGE": 0, "H-BOND(compatible)": 1, "HYDROPHOBIC": 2,
            "polar/vdW": 3}
    rows = []
    for key, al in pairs.items():
        ach, anum, ares, bch, bnum, bres = key
        best = {}
        for aa, ael, ba, bel, d in al:
            cls = None
            qa, qb = charge_of(ares, aa), charge_of(bres, ba)
            if qa and qb and qa * qb < 0 and d < SB_CUT:
                cls = "SALT BRIDGE"
            elif (ael in "NO" and bel in "NO" and d < HB_CUT
                  and ((is_donor(ares, aa) and is_acceptor(bres, ba))
                       or (is_acceptor(ares, aa) and is_donor(bres, ba)))):
                cls = "H-BOND(compatible)"
            elif ael in ("C", "S") and bel in ("C", "S") and d < HPHOB_CUT:
                cls = "HYDROPHOBIC"
            else:
                cls = "polar/vdW"
            if cls not in best or d < best[cls][2]:
                best[cls] = (aa, ba, d)
        top = min(best, key=lambda c: (rank[c], best[c][2]))
        mind = min(t[4] for t in al)
        ncontact = len(al)
        rows.append((ach, anum, ares, bch, bnum, bres, mind, ncontact, top,
                     best))

    rows.sort(key=lambda r: (rank[r[8]], r[6]))
    print(f"  {len(rows)} residue pairs within {cutoff} A\n")
    print(f"  {'antibody':<16}{'antigen':<14}{'min d':>7}{'atoms':>6}  "
          f"{'class':<20} justifying atom pair")
    counts = defaultdict(int)
    for ach, anum, ares, bch, bnum, bres, mind, nc, top, best in rows:
        aa, ba, d = best[top]
        counts[top] += 1
        print(f"  {ach+' '+ares+str(anum):<16}{bch+' '+bres+str(bnum):<14}"
              f"{mind:>7.2f}{nc:>6}  {top:<20} {aa}...{ba} {d:.2f} A")
    print("\n  class totals: " + ", ".join(
        f"{k} {v}" for k, v in sorted(counts.items(), key=lambda t: rank[t[0]])))

    # roll-up per antibody residue
    print(f"\n=== PER-PARATOPE-RESIDUE CHEMISTRY ROLL-UP ===")
    per = defaultdict(lambda: defaultdict(int))
    for ach, anum, ares, bch, bnum, bres, mind, nc, top, best in rows:
        per[(ach, anum, ares)][top] += 1
    print(f"  {'residue':<16}{'salt':>6}{'hbond':>7}{'phobic':>8}{'vdW':>6}"
          f"   partners")
    for k in sorted(per, key=lambda t: (t[0], t[1])):
        c = per[k]
        partners = sorted({(r[4], r[5]) for r in rows
                           if (r[0], r[1], r[2]) == k})
        print(f"  {k[0]+' '+k[2]+str(k[1]):<16}"
              f"{c['SALT BRIDGE']:>6}{c['H-BOND(compatible)']:>7}"
              f"{c['HYDROPHOBIC']:>8}{c['polar/vdW']:>6}   "
              + ", ".join(f"{r}{n}" for n, r in partners))


if __name__ == "__main__":
    main()
