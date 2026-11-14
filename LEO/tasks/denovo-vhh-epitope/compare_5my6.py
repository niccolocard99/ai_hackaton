"""Where does the 5MY6 nanobody (2Rs15d) bind HER2, and does it touch trastuzumab's epitope (1N8Z)?

1. Contacts between the nanobody (chain B) and HER2 (chain A, plus chain C glycan) at 4.5 A.
2. Residue-number offset between 5MY6 (UniProt numbering) and the organisers' 1N8Z crop (mature).
3. Superpose HER2 domain IV of 5MY6 onto 1N8Z, then compare footprints and test the nanobody for
   steric overlap with trastuzumab's Fv (chains A, B of 1N8Z).
Also re-derives trastuzumab's epitope with the same code as a cross-check (expect 19 residues).

Why this exists: 5MY6 was proposed as the VHH-HER2 structure to take hotspots from. This shows its nanobody
(2Rs15d) binds HER2 domain I, not trastuzumab's site in domain IV, so it cannot define hotspots for a
VHH that must compete with trastuzumab. Needs `gemmi` (pip install gemmi).

Fetch 5MY6 into its own new directory first (untrusted download, read as data only):
    mkdir dl_5my6 && curl -fsSL -o dl_5my6/5MY6.cif https://files.rcsb.org/download/5MY6.cif
Usage: python -I -X utf8 compare_5my6.py dl_5my6/5MY6.cif <repo>/apheris/rfd4/Data/12-her2/her2_fv_1n8z.cif
"""
import math
import re
import sys
from collections import defaultdict

import gemmi

CUT = 4.5


def load(path):
    st = gemmi.read_structure(path)
    st.setup_entities()
    return st[0]


def atoms_of(chain, seq_index=False):
    """(resnum, resname, atom, pos) for heavy atoms. seq_index=True numbers residues 1..n in chain order
    (author numbering can carry insertion codes or gaps)."""
    out = []
    for i, res in enumerate(chain, 1):
        if res.name == "HOH":
            continue
        num = i if seq_index else res.seqid.num
        for a in res:
            if a.element.name == "H":
                continue
            out.append((num, res.name, a.name, a.pos))
    return out


def grid(atoms, cell):
    g = defaultdict(list)
    for i, (_, _, _, p) in enumerate(atoms):
        g[(int(p.x // cell), int(p.y // cell), int(p.z // cell))].append(i)
    return g


def contacts(binder, target, cut=CUT):
    """{(target_resnum, target_resname): {binder_resnum: min_dist}}, plus atom names touching."""
    g = grid(target, cut)
    pairs = {}
    names = defaultdict(set)
    for bn, brn, ban, bp in binder:
        cx, cy, cz = int(bp.x // cut), int(bp.y // cut), int(bp.z // cut)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for ti in g.get((cx + dx, cy + dy, cz + dz), ()):
                        tn, trn, tan, tp = target[ti]
                        d = bp.dist(tp)
                        if d <= cut:
                            k = ((tn, trn), (bn, brn))
                            if k not in pairs or d < pairs[k]:
                                pairs[k] = d
                            names[(tn, trn)].add(tan)
    return pairs, names


def by_target(pairs):
    out = defaultdict(dict)
    for (t, b), d in pairs.items():
        out[t][b] = d
    return out


def by_binder(pairs):
    out = defaultdict(dict)
    for (t, b), d in pairs.items():
        out[b][t] = d
    return out


ONE = {r: gemmi.find_tabulated_residue(r).one_letter_code.upper()
       for r in ("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL").split()}


def cdr_labels(seq):
    """1-based index -> label for a VHH, from sequence landmarks (heuristic, fine for reporting)."""
    lab = {i + 1: "FR" for i in range(len(seq))}
    c1 = seq.find("C")
    w = seq.find("W", c1)
    m3 = re.search(r"[YF][YF]C", seq)
    j = seq.find("RFT", w)
    e3 = seq.find("WGQG")
    spans = {"CDR1": (c1 + 4, w - 1)}
    if j > 0:
        spans["CDR2"] = (w + 14, j - 1)
    if m3 and e3 > 0:
        spans["CDR3"] = (m3.end() + 2, e3 - 1)
    for name, (a, b) in spans.items():
        for i in range(a, b + 1):
            lab[i + 1] = name
    return lab, spans


def centroid(points):
    n = len(points)
    return gemmi.Position(sum(p.x for p in points) / n, sum(p.y for p in points) / n, sum(p.z for p in points) / n)


m5 = load(sys.argv[1])
m1 = load(sys.argv[2])

# ---------------------------------------------------------------- 1. 5MY6 interface
A5, B5, C5 = m5["A"], m5["B"], m5["C"]
her2_5 = atoms_of(A5)
vhh_5 = atoms_of(B5, seq_index=True)
glyc_5 = atoms_of(C5)
vhh_seq = "".join(ONE.get(r.name, "X") for r in B5 if r.het_flag == "A")
auth_b = {i: f"{r.seqid.num}{r.seqid.icode.strip()}" for i, r in enumerate(B5, 1)}
dup = len(set(r.seqid.num for r in A5)) != len(A5)
print(f"nanobody chain: {len(B5)} residues, author labels {auth_b[1]}..{auth_b[len(B5)]}; "
      f"HER2 chain has duplicate numbers: {dup}")
labels, spans = cdr_labels(vhh_seq)

pairs5, names5 = contacts(vhh_5, her2_5)
t5 = by_target(pairs5)
b5 = by_binder(pairs5)

# offset between 5MY6 author numbering and the 1N8Z (mature) numbering
d5 = {r.seqid.num: r.name for r in A5 if r.het_flag == "A"}  # polymer residues only (no waters/glycans)
d1 = {}
for ch in ("C", "D"):
    for r in m1[ch]:
        d1[r.seqid.num] = r.name
best = max(range(-60, 61), key=lambda off: sum(1 for a, rn in d5.items() if d1.get(a + off) == rn))
shared = [a for a in d5 if d1.get(a + best) == d5[a]]
tot_over = sum(1 for a in d5 if (a + best) in d1)
OFF = best
print(f"=== numbering: mature = 5MY6 author {OFF:+d}; {len(shared)}/{tot_over} overlapping residues identical ===")
print(f"5MY6 HER2 modelled: author {min(d5)}..{max(d5)} (mature {min(d5)+OFF}..{max(d5)+OFF}); "
      f"1N8Z crop: {min(d1)}..{max(d1)}")

print(f"\n=== nanobody 2Rs15d: {len(vhh_seq)} aa, CDR spans (1-based): "
      + ", ".join(f"{k} {a+1}-{b+1} {vhh_seq[a:b+1]}" for k, (a, b) in spans.items()))

print(f"\n=== 5MY6 epitope: HER2 residues within {CUT} A of the nanobody ===")
print(f"{'author':>7} {'mature':>7} {'aa':>3} {'n_vhh':>6}  nanobody residues (region)")
for (tn, trn), partners in sorted(t5.items()):
    desc = ", ".join(f"{ONE.get(brn, brn)}{auth_b[bn]}({labels[bn]})" for (bn, brn) in sorted(partners))
    print(f"{tn:>7} {tn+OFF:>7} {ONE.get(trn, trn):>3} {len(partners):>6}  {desc}")
print(f"-> {len(t5)} HER2 residues; mature numbering span {min(k[0] for k in t5)+OFF}..{max(k[0] for k in t5)+OFF}")

print(f"\n=== nanobody paratope (residues touching HER2) ===")
for (bn, brn), partners in sorted(b5.items()):
    print(f"  {ONE.get(brn, brn)}{auth_b[bn]:<5} {labels[bn]:<5} touches {len(partners)} HER2 residues")
reg = defaultdict(int)
for (bn, brn), partners in b5.items():
    reg[labels[bn]] += len(partners)
tot = sum(reg.values())
print("contact share by region: " + ", ".join(f"{k} {v}/{tot} ({100*v/tot:.0f}%)" for k, v in sorted(reg.items())))

gp, _ = contacts(vhh_5, glyc_5)
print(f"\nnanobody-glycan (chain C NAG) contacts within {CUT} A: {len(gp)} residue pairs")

# ---------------------------------------------------------------- 2. trastuzumab epitope, same code
fv_1 = atoms_of(m1["A"]) + atoms_of(m1["B"])
her2_1 = atoms_of(m1["C"]) + atoms_of(m1["D"])
pairs1, _ = contacts(fv_1, her2_1)
t1 = by_target(pairs1)
print(f"\n=== cross-check: trastuzumab epitope by this code = {len(t1)} HER2 residues (earlier measurement: 19) ===")

# ---------------------------------------------------------------- 3. superpose 5MY6 onto 1N8Z
ref, mov = [], []
for r in m1["C"]:
    a = r.seqid.num - OFF
    if a in d5 and d5[a] == r.name:
        ca1 = r.find_atom("CA", "*")
        res5 = A5[str(a)][0] if False else None
        for r5 in A5:
            if r5.seqid.num == a:
                ca5 = r5.find_atom("CA", "*")
                if ca1 and ca5:
                    ref.append(ca1.pos)
                    mov.append(ca5.pos)
                break
print(f"\n=== superposition on {len(ref)} shared domain IV CA atoms ===")
sup = gemmi.superpose_positions(ref, mov)


def rmsd(a, b):
    return math.sqrt(sum(p.dist(q) ** 2 for p, q in zip(a, b)) / len(a))


def mv(tf, p):
    v = tf.apply(p)  # gemmi returns Vec3
    return gemmi.Position(v.x, v.y, v.z)


moved = [mv(sup.transform, p) for p in mov]
r_fwd = rmsd(ref, moved)
inv = sup.transform.inverse()
moved_inv = [mv(inv, p) for p in mov]
r_inv = rmsd(ref, moved_inv)
tf = sup.transform if r_fwd <= r_inv else inv
print(f"gemmi reports rmsd {sup.rmsd:.2f}; direct check: forward {r_fwd:.2f}, inverse {r_inv:.2f}; using "
      f"{'forward' if tf is sup.transform else 'inverse'}")

vhh_in_1n8z = [(n, rn, an, mv(tf, p)) for (n, rn, an, p) in vhh_5]

# footprint comparison (HER2 residues in mature numbering)
ep5 = {k[0] + OFF for k in t5}
ep1 = {k[0] for k in t1}
present = {n for n in ep1 if (n - OFF) in d5}
inter = sorted(ep5 & ep1)
print(f"\n=== footprint overlap (mature numbering) ===")
print(f"2Rs15d epitope : {len(ep5)} residues")
print(f"trastuzumab    : {len(ep1)} residues, of which {len(present)} are modelled in 5MY6")
print(f"shared residues: {len(inter)} -> {inter}")

# steric test: nanobody (in 1N8Z frame) against trastuzumab Fv
gfv = grid(fv_1, 4.5)
n2 = n3 = n4 = 0
touch_res = set()
dmin = 1e9
for bn, brn, ban, bp in vhh_in_1n8z:
    cx, cy, cz = int(bp.x // 4.5), int(bp.y // 4.5), int(bp.z // 4.5)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                for ti in gfv.get((cx + dx, cy + dy, cz + dz), ()):
                    d = bp.dist(fv_1[ti][3])
                    dmin = min(dmin, d)
                    if d <= 2.0:
                        n2 += 1
                    if d <= 3.0:
                        n3 += 1
                        touch_res.add(bn)
                    if d <= 4.0:
                        n4 += 1
print(f"\n=== steric test: nanobody placed in the 1N8Z frame vs trastuzumab Fv ===")
print(f"grid search (4.5 A cells): atom pairs <=2.0 A: {n2}; <=3.0 A: {n3}; <=4.0 A: {n4}; "
      f"nanobody residues with an atom <=3.0 A of the Fv: {len(touch_res)}")

true_min = min(bp.dist(fp) for (_, _, _, bp) in vhh_in_1n8z for (_, _, _, fp) in fv_1)
print(f"true closest approach, any nanobody atom to any trastuzumab Fv atom (all pairs): {true_min:.1f} A")
c5 = centroid([a[3] for a in her2_5 if (a[0] in {k[0] for k in t5}) and a[2] == "CA"])
c1 = centroid([a[3] for a in her2_1 if (a[0] in ep1) and a[2] == "CA"])
c5t = mv(tf, c5)
print(f"distance between the two epitope centroids (CA, 1N8Z frame): {c5t.dist(c1):.1f} A")
