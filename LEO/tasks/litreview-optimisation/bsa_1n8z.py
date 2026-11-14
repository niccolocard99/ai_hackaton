#!/usr/bin/env python3
"""Buried solvent-accessible surface area (BSA) across an antibody-antigen interface.

Written for the trastuzumab-HER2 interface of PDB 1N8Z, but generic: the chain
groups are command-line arguments, so the same code runs on the cryo-EM
cross-check structures.

    python bsa_1n8z.py <file.cif> <groupA_chains> <groupB_chains> [--cdr-chains H=B,L=A]

e.g.  python bsa_1n8z.py 1n8z.cif AB C --cdr-chains H=B,L=A

WHY THIS FILE EXISTS
--------------------
Neither biopython nor freesasa is installed in this environment and nothing may
be installed, so SASA is computed here from scratch with the **Shrake-Rupley
numerical method** (J Mol Biol 1973;79:351): each atom is given a sphere of
radius r_vdw + r_probe, that sphere is sampled with a fixed quasi-uniform point
set, and the accessible fraction is the share of sample points that fall outside
every neighbouring atom's expanded sphere.

    SASA_i = 4*pi*(r_i + r_probe)^2 * (accessible points / total points)

Point set: golden-spiral (Fibonacci) sphere sampling, 960 points/atom.
Probe radius: 1.40 A (water).
Radii: heavy-atom (united-atom) van der Waals radii, hydrogens absent from the
model and therefore folded into the radii, as is standard for X-ray coordinates.

BSA per residue is the standard difference definition:

    BSA(residue) = SASA(residue in its own chain group, alone)
                 - SASA(residue in the assembled complex)

Summing BSA over both sides of the interface gives the total area buried on
assembly; the conventional single number quoted as "interface area" is half
that (the area of one face). Both are printed, labelled.

Everything printed is MEASURED FROM THE STRUCTURE. No energies are computed
here and none are implied: BSA is a geometric quantity.
"""
import math
import sys
from collections import defaultdict

import numpy as np

# ---------------------------------------------------------------- constants

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}

# United-atom van der Waals radii, A. Shrake-Rupley / Chothia values.
VDW = {
    "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "SE": 1.90,
    "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98, "H": 1.20, "D": 1.20,
}
VDW_DEFAULT = 1.70

PROBE = 1.40
N_POINTS = 960

# Canonical trastuzumab CDRs (Kabat content). Located by substring search on
# each structure's own deposited sequence, so no numbering scheme is assumed
# and the same definitions transfer between structures with different numbering.
CDR_MOTIFS = {
    "H1": ("H", "DTYIH"),
    "H2": ("H", "RIYPTNGYTRYADSVK"),
    "H3": ("H", "WGGDGFYAMDY"),
    "L1": ("L", "RASQDVNTAVA"),
    "L2": ("L", "SASFLYS"),
    "L3": ("L", "QQHYTTPPT"),
}

# ---------------------------------------------------------------- mmCIF input


def parse_atom_site(path):
    """Yield one dict per _atom_site row of an mmCIF file.

    Tolerant of quoted values and of the loop appearing anywhere in the file.
    """
    cols, rows, in_loop, in_atom = [], [], False, False
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s == "loop_":
                if rows:
                    break
                in_loop, in_atom, cols = True, False, []
                continue
            if in_loop and s.startswith("_atom_site."):
                cols.append(s.split(".", 1)[1])
                in_atom = True
                continue
            if in_atom and cols:
                if not s or s == "#" or s.startswith("_") or s.startswith("loop_"):
                    if rows:
                        break
                    in_loop = in_atom = False
                    continue
                parts = _split_cif(s)
                if len(parts) != len(cols):
                    continue
                rows.append(dict(zip(cols, parts)))
    return rows


def _split_cif(line):
    """Whitespace split that respects single/double quoted CIF values."""
    out, buf, quote = [], [], None
    for ch in line:
        if quote:
            if ch == quote:
                quote = None
                out.append("".join(buf))
                buf = []
            else:
                buf.append(ch)
        elif ch in "'\"":
            quote = ch
        elif ch.isspace():
            if buf:
                out.append("".join(buf))
                buf = []
        else:
            buf.append(ch)
    if buf:
        out.append("".join(buf))
    return out


def load_atoms(path, altloc_first=True):
    """Return (atoms, seqs).

    atoms: list of (chain, resnum, insertion_code, resname, atomname, element, xyz)
    seqs:  {chain: {resnum: one-letter}} for the standard residues present.
    """
    rows = parse_atom_site(path)
    if not rows:
        sys.exit(f"no _atom_site rows parsed from {path}")
    atoms, seqs, seen = [], defaultdict(dict), set()
    for r in rows:
        if r.get("group_PDB") != "ATOM":
            continue
        if r.get("pdbx_PDB_model_num", "1") not in ("1", "."):
            continue
        el = (r.get("type_symbol") or "").upper()
        if el in ("H", "D"):
            continue
        comp = r.get("label_comp_id", "")
        if comp not in THREE_TO_ONE:
            continue
        alt = r.get("label_alt_id", ".")
        chain = r.get("auth_asym_id") or r.get("label_asym_id")
        seqid = r.get("auth_seq_id") or r.get("label_seq_id")
        icode = r.get("pdbx_PDB_ins_code", ".")
        icode = "" if icode in (".", "?") else icode
        name = r.get("label_atom_id", "")
        try:
            xyz = (float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"]))
            num = int(seqid)
        except (KeyError, ValueError):
            continue
        key = (chain, num, icode, name)
        if altloc_first:
            if key in seen:
                continue  # keep only the first altloc of each atom
            seen.add(key)
        elif alt not in (".", "?", "A"):
            continue
        atoms.append((chain, num, icode, comp, name, el, xyz))
        seqs[chain][num] = THREE_TO_ONE[comp]
    return atoms, seqs

# ---------------------------------------------------------------- Shrake-Rupley


def sphere_points(n):
    """n quasi-uniform points on the unit sphere (golden-spiral / Fibonacci)."""
    i = np.arange(n, dtype=float) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)
    theta = math.pi * (1.0 + 5.0 ** 0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1)


def sasa(coords, radii, probe=PROBE, n_points=N_POINTS):
    """Per-atom SASA in A^2 by the Shrake-Rupley numerical method.

    coords (m,3), radii (m,) van der Waals radii. Neighbour search uses a
    uniform cell list so the cost is linear in the number of atoms.
    """
    coords = np.asarray(coords, dtype=float)
    R = np.asarray(radii, dtype=float) + probe
    m = len(coords)
    unit = sphere_points(n_points)
    const = 4.0 * math.pi / n_points

    cell = 2.0 * R.max()
    keys = np.floor(coords / cell).astype(np.int64)
    grid = defaultdict(list)
    for idx in range(m):
        grid[(keys[idx, 0], keys[idx, 1], keys[idx, 2])].append(idx)

    out = np.empty(m)
    offs = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1)
            for dz in (-1, 0, 1)]
    for i in range(m):
        kx, ky, kz = keys[i]
        cand = []
        for dx, dy, dz in offs:
            cand.extend(grid.get((kx + dx, ky + dy, kz + dz), ()))
        cand = np.fromiter((j for j in cand if j != i), dtype=np.int64)
        if cand.size:
            d = np.linalg.norm(coords[cand] - coords[i], axis=1)
            cand = cand[d < (R[i] + R[cand])]
        if cand.size == 0:
            out[i] = const * R[i] * R[i] * n_points
            continue
        test = coords[i] + unit * R[i]                       # (n,3)
        diff = test[:, None, :] - coords[cand][None, :, :]    # (n,k,3)
        d2 = np.einsum("nkd,nkd->nk", diff, diff)
        buried = (d2 < (R[cand] ** 2)[None, :]).any(axis=1)
        out[i] = const * R[i] * R[i] * int(n_points - buried.sum())
    return out


def selftest():
    """The method must reproduce the analytic SASA of an isolated atom."""
    r = VDW["C"]
    got = sasa(np.array([[0.0, 0.0, 0.0]]), np.array([r]))[0]
    want = 4.0 * math.pi * (r + PROBE) ** 2
    err = abs(got - want) / want
    print(f"  self-test, isolated carbon: computed {got:.2f} A^2, "
          f"analytic 4*pi*(1.70+1.40)^2 = {want:.2f} A^2, error {err:.2e}")
    assert err < 1e-9, "sphere sampling is wrong"
    # Two bonded carbons: compare against the exact spherical-cap result.
    # Each expanded sphere loses the cap of height h = R - d/2 (equal radii),
    # area 2*pi*R*h, so SASA_atom = 4*pi*R^2 - 2*pi*R*h.
    d = 1.52
    R = r + PROBE
    h = R - d / 2.0
    exact = 2.0 * (4.0 * math.pi * R * R - 2.0 * math.pi * R * h)
    got2 = sasa(np.array([[0.0, 0, 0], [d, 0, 0]]), np.array([r, r])).sum()
    err2 = abs(got2 - exact) / exact
    print(f"  self-test, two carbons 1.52 A apart: computed {got2:.2f} A^2, "
          f"exact spherical-cap {exact:.2f} A^2, error {err2:.2%}")
    assert err2 < 0.01, "sphere sampling is too coarse"

# ---------------------------------------------------------------- residue roll-up


def residue_sasa(atoms, per_atom):
    """{(chain,resnum,icode,resname): SASA} from per-atom values."""
    agg = defaultdict(float)
    for (ch, num, ic, comp, _name, _el, _xyz), a in zip(atoms, per_atom):
        agg[(ch, num, ic, comp)] += a
    return agg


def locate_cdrs(seqs, cdr_chain_map):
    """{cdr: (chain, lo, hi, motif)} by substring search on each chain sequence.

    cdr_chain_map maps the logical chain ('H'/'L') to this structure's author
    chain id. Residue numbers returned are this structure's author numbering.
    """
    ranges = {}
    for name, (logical, motif) in CDR_MOTIFS.items():
        chain = cdr_chain_map.get(logical)
        if chain is None or chain not in seqs:
            continue
        nums = sorted(seqs[chain])
        # Sequence as deposited, in author-numbering order, with an index back
        # to the author number of every position (gaps simply absent).
        seq = "".join(seqs[chain][n] for n in nums)
        idx = seq.find(motif)
        if idx < 0:
            continue
        if seq.find(motif, idx + 1) > 0:
            print(f"  ! CDR {name} motif is not unique in chain {chain}; skipped")
            continue
        auth = nums[idx:idx + len(motif)]
        ranges[name] = (chain, auth[0], auth[-1], motif, auth)
    return ranges

# ---------------------------------------------------------------- main


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    path, gA, gB = sys.argv[1], sys.argv[2], sys.argv[3]
    # Default: 1N8Z, where chain B is heavy and chain A is light.
    cdr_map = {"H": "B", "L": "A"}
    rest = sys.argv[4:]
    if "--cdr-chains" in rest:
        spec = rest[rest.index("--cdr-chains") + 1]
        cdr_map = {}
        for part in spec.split(","):
            k, v = part.split("=", 1)
            cdr_map[k.strip()] = v.strip()

    # Chain groups may be given as "AB" (single-char ids) or "A,B".
    setA = set(gA.split(",")) if "," in gA else set(gA)
    setB = set(gB.split(",")) if "," in gB else set(gB)

    print(f"=== BSA by Shrake-Rupley, implemented here (no external SASA library) ===")
    print(f"  file: {path}")
    print(f"  group A (antibody): {sorted(setA)}   group B (antigen): {sorted(setB)}")
    print(f"  probe radius {PROBE} A, {N_POINTS} sample points/atom, "
          f"heavy atoms only (united-atom radii)")
    selftest()

    atoms, seqs = load_atoms(path)
    keep = [a for a in atoms if a[0] in setA | setB]
    if not keep:
        sys.exit("no atoms in the requested chains")
    idxA = [i for i, a in enumerate(keep) if a[0] in setA]
    idxB = [i for i, a in enumerate(keep) if a[0] in setB]
    if not idxA or not idxB:
        sys.exit(f"one group is empty: |A|={len(idxA)} |B|={len(idxB)}")
    print(f"\n  atoms used: {len(keep)} total "
          f"({len(idxA)} in A, {len(idxB)} in B)")
    for ch in sorted(setA | setB):
        nums = sorted(seqs.get(ch, {}))
        if nums:
            print(f"    chain {ch}: {len(nums)} residues, "
                  f"author numbering {nums[0]}..{nums[-1]}")

    coords = np.array([a[6] for a in keep])
    radii = np.array([VDW.get(a[5], VDW_DEFAULT) for a in keep])

    print("\n  running SASA on the complex ...", flush=True)
    s_cplx = sasa(coords, radii)
    print("  running SASA on group A alone ...", flush=True)
    sA = sasa(coords[idxA], radii[idxA])
    print("  running SASA on group B alone ...", flush=True)
    sB = sasa(coords[idxB], radii[idxB])

    free = np.zeros(len(keep))
    free[idxA] = sA
    free[idxB] = sB
    dsasa = free - s_cplx          # >= 0; area buried on assembly

    res_free = residue_sasa(keep, free)
    res_cplx = residue_sasa(keep, s_cplx)
    res_bsa = {k: res_free[k] - res_cplx[k] for k in res_free}

    totA = sum(v for k, v in res_bsa.items() if k[0] in setA)
    totB = sum(v for k, v in res_bsa.items() if k[0] in setB)
    print(f"\n=== TOTAL INTERFACE BSA (measured from structure) ===")
    print(f"  buried on the antibody side (group A): {totA:8.1f} A^2")
    print(f"  buried on the antigen  side (group B): {totB:8.1f} A^2")
    print(f"  total area buried on assembly        : {totA+totB:8.1f} A^2")
    print(f"  interface area (one face, = total/2) : {(totA+totB)/2:8.1f} A^2")

    print(f"\n=== BSA PER CHAIN ===")
    print(f"  {'chain':<7}{'BSA A^2':>10}{'share of own side':>20}  residues buried >1 A^2")
    for ch in sorted(setA | setB):
        tot = sum(v for k, v in res_bsa.items() if k[0] == ch)
        side = totA if ch in setA else totB
        nres = sum(1 for k, v in res_bsa.items() if k[0] == ch and v > 1.0)
        pct = 100 * tot / side if side else 0.0
        print(f"  {ch:<7}{tot:>10.1f}{pct:>19.1f}%  {nres}")

    ranges = locate_cdrs(seqs, cdr_map)
    if ranges:
        print(f"\n=== CDR LOCATIONS (this structure's author numbering) ===")
        for name in ("H1", "H2", "H3", "L1", "L2", "L3"):
            if name in ranges:
                ch, lo, hi, motif, _auth = ranges[name]
                print(f"  CDR-{name}: chain {ch} {lo}-{hi}  {motif}")

        def cdr_of(ch, num):
            for name, (c, lo, hi, _m, auth) in ranges.items():
                if c == ch and num in auth:
                    return name
            return None

        print(f"\n=== BSA PER CDR (antibody side) ===")
        by = defaultdict(float)
        nres = defaultdict(int)
        members = defaultdict(list)
        fw = 0.0
        for (ch, num, ic, comp), v in res_bsa.items():
            if ch not in setA:
                continue
            name = cdr_of(ch, num)
            if name is None:
                fw += v
                continue
            by[name] += v
            if v > 1.0:
                nres[name] += 1
                members[name].append((num, comp, v))
        print(f"  {'CDR':<6}{'BSA A^2':>10}{'% of paratope BSA':>20}{'res>1A^2':>10}  "
              f"contributing residues (BSA A^2)")
        for name in sorted(by, key=lambda k: -by[k]):
            lst = sorted(members[name], key=lambda t: -t[2])
            desc = ", ".join(f"{c}{n} {v:.0f}" for n, c, v in lst)
            print(f"  {name:<6}{by[name]:>10.1f}{100*by[name]/totA:>19.1f}%"
                  f"{nres[name]:>10}  {desc}")
        print(f"  {'FRWK':<6}{fw:>10.1f}{100*fw/totA:>19.1f}%"
              f"{'':>10}  (framework / constant region)")
        hv = sum(v for k, v in by.items() if k.startswith("H"))
        lt = sum(v for k, v in by.items() if k.startswith("L"))
        cdr_tot = hv + lt
        print(f"\n  heavy-chain CDRs {hv:.1f} A^2 "
              f"({100*hv/cdr_tot:.0f}% of CDR BSA) | "
              f"light-chain CDRs {lt:.1f} A^2 ({100*lt/cdr_tot:.0f}%)")
        print(f"  all six CDRs {cdr_tot:.1f} A^2 = "
              f"{100*cdr_tot/totA:.0f}% of the antibody-side BSA")

        if "H3" in ranges:
            ch, lo, hi, motif, auth = ranges["H3"]
            print(f"\n=== BSA PER CDRH3 POSITION ({motif}, chain {ch} {lo}-{hi}) ===")
            print(f"  {'pos':>5} {'aa':>3} {'BSA A^2':>9} {'free SASA':>10} "
                  f"{'% of its SASA buried':>22}")
            for aa, num in zip(motif, auth):
                k = [kk for kk in res_bsa if kk[0] == ch and kk[1] == num]
                if not k:
                    continue
                k = k[0]
                b, f = res_bsa[k], res_free[k]
                pct = 100 * b / f if f > 0.5 else 0.0
                flag = "  <== buried" if b > 10 else ("  (trace)" if b > 1 else "")
                print(f"  {num:>5} {aa:>3} {b:>9.1f} {f:>10.1f} {pct:>21.1f}%{flag}")

    print(f"\n=== ANTIGEN-SIDE RESIDUES, BSA > 1 A^2 ===")
    lst = sorted(((k, v) for k, v in res_bsa.items() if k[0] in setB and v > 1.0),
                 key=lambda t: -t[1])
    print(f"  {len(lst)} residues; top by BSA:")
    for k, v in lst:
        print(f"    {k[0]} {k[3]}{k[1]}{k[2]:<2} {v:8.1f} A^2")


if __name__ == "__main__":
    main()
