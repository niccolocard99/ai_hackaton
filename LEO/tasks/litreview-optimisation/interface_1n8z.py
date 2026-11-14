#!/usr/bin/env python3
"""Compute the trastuzumab(Herceptin)-HER2 interface directly from PDB 1N8Z.

Usage: python interface_1n8z.py <path-to-1n8z.cif> [cutoff_angstrom]

Why this exists: the literature abstracts in this task's cache never list
trastuzumab's paratope residues or its CDR sequences. The deposited structure
does, so we measure it rather than quote it. No external dependencies.

Outputs: chain inventory, one-letter sequences, and every inter-chain residue
contact pair under the cutoff (default 4.5 A, heavy atoms only).
"""
import math
import sys
from collections import defaultdict

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}


def parse_atom_site(path):
    """Yield dicts for each _atom_site row of an mmCIF file."""
    cols, rows, in_loop, in_atom = [], [], False, False
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s == "loop_":
                in_loop, in_atom, cols = True, False, []
                continue
            if in_loop and s.startswith("_atom_site."):
                cols.append(s.split(".", 1)[1])
                in_atom = True
                continue
            if in_atom and cols:
                if s.startswith("_") or s == "#" or not s:
                    if rows:
                        break
                    in_loop = in_atom = False
                    continue
                parts = s.split()
                if len(parts) != len(cols):
                    continue
                rows.append(dict(zip(cols, parts)))
    return rows


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    path = sys.argv[1]
    cutoff = float(sys.argv[2]) if len(sys.argv) > 2 else 4.5

    rows = parse_atom_site(path)
    if not rows:
        sys.exit("no _atom_site rows parsed - is this an mmCIF file?")

    # Heavy atoms of standard residues, first model only.
    atoms = []
    seqs = defaultdict(dict)
    for r in rows:
        if r.get("group_PDB") != "ATOM":
            continue
        if r.get("pdbx_PDB_model_num", "1") != "1":
            continue
        el = r.get("type_symbol", "")
        if el == "H" or el == "D":
            continue
        comp = r.get("label_comp_id", "")
        if comp not in THREE_TO_ONE:
            continue
        chain = r.get("auth_asym_id") or r.get("label_asym_id")
        seqid = r.get("auth_seq_id") or r.get("label_seq_id")
        try:
            xyz = (float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"]))
            num = int(seqid)
        except (KeyError, ValueError):
            continue
        atoms.append((chain, num, comp, xyz))
        seqs[chain][num] = THREE_TO_ONE[comp]

    print(f"parsed {len(atoms)} heavy atoms, cutoff {cutoff} A\n")
    print("=== CHAIN INVENTORY ===")
    for ch in sorted(seqs):
        nums = sorted(seqs[ch])
        print(f"chain {ch}: {len(nums)} residues, numbering {nums[0]}..{nums[-1]}")
    print()

    print("=== SEQUENCES (auth numbering, gaps marked '-') ===")
    for ch in sorted(seqs):
        nums = sorted(seqs[ch])
        seq = "".join(seqs[ch].get(i, "-") for i in range(nums[0], nums[-1] + 1))
        print(f">chain_{ch} start={nums[0]} end={nums[-1]} len={len(nums)}")
        for i in range(0, len(seq), 60):
            print(f"  {nums[0]+i:>4} {seq[i:i+60]}")
    print()

    # Spatial hash so we do not do an O(n^2) all-pairs comparison.
    cell = cutoff
    grid = defaultdict(list)
    for idx, (ch, num, comp, xyz) in enumerate(atoms):
        key = (int(xyz[0] // cell), int(xyz[1] // cell), int(xyz[2] // cell))
        grid[key].append(idx)

    best = {}  # (chainA,resA,chainB,resB) -> min distance
    c2 = cutoff * cutoff
    for idx, (ch, num, comp, xyz) in enumerate(atoms):
        kx, ky, kz = int(xyz[0] // cell), int(xyz[1] // cell), int(xyz[2] // cell)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for jdx in grid.get((kx + dx, ky + dy, kz + dz), ()):
                        if jdx <= idx:
                            continue
                        ch2, num2, comp2, xyz2 = atoms[jdx]
                        if ch2 == ch:
                            continue
                        d2 = ((xyz[0] - xyz2[0]) ** 2 + (xyz[1] - xyz2[1]) ** 2
                              + (xyz[2] - xyz2[2]) ** 2)
                        if d2 > c2:
                            continue
                        a, b = (ch, num, comp), (ch2, num2, comp2)
                        if a > b:
                            a, b = b, a
                        key = (a, b)
                        d = math.sqrt(d2)
                        if key not in best or d < best[key]:
                            best[key] = d

    pairs_by_chains = defaultdict(list)
    for (a, b), d in best.items():
        pairs_by_chains[(a[0], b[0])].append((a, b, d))

    print("=== INTER-CHAIN RESIDUE CONTACTS ===")
    for chains in sorted(pairs_by_chains):
        lst = sorted(pairs_by_chains[chains], key=lambda t: (t[0][1], t[1][1]))
        print(f"\n--- chains {chains[0]} <-> {chains[1]} : "
              f"{len(lst)} residue pairs ---")
        resA = sorted({(p[0][1], p[0][2]) for p in lst})
        resB = sorted({(p[1][1], p[1][2]) for p in lst})
        print(f"  {chains[0]} interface residues ({len(resA)}): "
              + ", ".join(f"{c}{n}" for n, c in resA))
        print(f"  {chains[1]} interface residues ({len(resB)}): "
              + ", ".join(f"{c}{n}" for n, c in resB))
        print("  pairs (min heavy-atom distance, A):")
        for a, b, d in lst:
            print(f"    {chains[0]} {a[2]}{a[1]:<5} -- {chains[1]} {b[2]}{b[1]:<5} {d:.2f}")


if __name__ == "__main__":
    main()
