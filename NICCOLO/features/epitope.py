import math, collections, sys

CIF = "/Users/niccoloc/Desktop/Other_projects/Hackaton/ai_hackaton/apheris/rfd4/Data/12-her2/her2_fv_1n8z.cif"

three = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
         'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
         'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

atoms = []
for l in open(CIF):
    if l.startswith('ATOM'):
        f = l.split()
        try:
            x, y, z = float(f[19]), float(f[20]), float(f[21])
        except (ValueError, IndexError):
            continue
        # chain, auth_seq, comp, atom, igregion, coords
        atoms.append((f[11], int(f[9]), f[4], f[2], f[17].strip("'"), (x, y, z)))

ab = [a for a in atoms if a[0] in ('A', 'B')]      # VL, VH
ag = [a for a in atoms if a[0] in ('C', 'D')]      # HER2 domain IV

CUT = 4.5
# epitope: HER2 residue -> {antibody chain -> {igregion -> pair count}}
epi = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
# paratope: antibody residue -> pair count
para = collections.Counter()

for a in ab:
    for g in ag:
        if math.dist(a[5], g[5]) < CUT:
            key = (g[0], g[1], three.get(g[2], 'X'))
            epi[key][a[0]][a[4] or 'framework'] += 1
            para[(a[0], a[1], three.get(a[2], 'X'), a[4] or 'framework')] += 1

print(f"=== trastuzumab-HER2 epitope, 1N8Z, heavy-atom cutoff {CUT} A ===\n")
print(f"{'HER2':>10}  {'pairs':>5}  {'VH':>4} {'VL':>4}   {'reachable by VHH?':<18} regions")
vh_only = vl_any = 0
for key in sorted(epi, key=lambda k: (k[0], k[1])):
    d = epi[key]
    nvh = sum(d.get('B', {}).values())
    nvl = sum(d.get('A', {}).values())
    regs = ','.join(sorted({r for c in d for r in d[c]}))
    if nvl == 0:
        verdict = "yes (VH only)"
        vh_only += 1
    elif nvh == 0:
        verdict = "NO (VL only)"
        vl_any += 1
    else:
        verdict = "partial (VH+VL)"
        vl_any += 1
    print(f"{key[0]}{key[1]:>4}{key[2]:>4}  {nvh+nvl:>5}  {nvh:>4} {nvl:>4}   {verdict:<18} {regs}")

print(f"\nepitope residues: {len(epi)}   VH-only: {vh_only}   involving VL: {vl_any}")
tot = sum(sum(sum(r.values()) for r in d.values()) for d in epi.values())
tvh = sum(sum(d.get('B', {}).values()) for d in epi.values())
print(f"contact pairs: {tot}   from VH: {tvh} ({100*tvh/tot:.0f}%)   from VL: {tot-tvh} ({100*(tot-tvh)/tot:.0f}%)")

print("\n=== paratope by region ===")
byreg = collections.Counter()
for (ch, num, aa, reg), n in para.items():
    byreg[(ch, reg)] += n
for k in sorted(byreg, key=lambda k: -byreg[k]):
    print(f"  chain {k[0]} {k[1]:<10} {byreg[k]:>4} pairs")

print("\n=== damage from the unresolved HER2 581-590 loop ===")
# how close is each epitope residue to the chain break (C ends 580, D starts 591)?
endpoints = [a for a in ag if (a[0] == 'C' and a[1] == 580) or (a[0] == 'D' and a[1] == 591)]
for key in sorted(epi, key=lambda k: (k[0], k[1])):
    res = [a for a in ag if a[0] == key[0] and a[1] == key[1]]
    dmin = min(math.dist(r[5], e[5]) for r in res for e in endpoints)
    seqgap = min(abs(key[1] - 580) if key[0] == 'C' else 999,
                 abs(key[1] - 591) if key[0] == 'D' else 999)
    flag = "  <-- neighbour counts biased low" if dmin < 12 else ""
    print(f"{key[0]}{key[1]:>4}{key[2]:>4}  {dmin:5.1f} A from the break, {seqgap:>3} res away in sequence{flag}")
