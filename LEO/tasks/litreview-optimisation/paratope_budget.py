#!/usr/bin/env python3
"""Per-CDR paratope budget for trastuzumab from the measured 1N8Z interface.

Usage: python paratope_budget.py interface_1n8z.txt

Reads the contact listing produced by interface_1n8z.py, assigns each
antibody residue to a CDR (CDRs located by substring search on the chain
sequences, so no numbering scheme is assumed), and reports how much of the
HER2 interface each CDR actually carries.

The question this answers: if CDRH3 must be rewritten to <70% identity, how
much of the binding interface is at risk, and which CDRH3 positions are free?
"""
import re
import sys
from collections import defaultdict

# Chain sequences as deposited in 1N8Z (auth numbering starts at 1).
LIGHT = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPS"
         "RFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
HEAVY = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRY"
         "ADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")

# Canonical trastuzumab CDR sequences (Kabat definition content).
CDRS = {
    "L1": ("A", "RASQDVNTAVA"),
    "L2": ("A", "SASFLYS"),
    "L3": ("A", "QQHYTTPPT"),
    "H1": ("B", "DTYIH"),
    "H2": ("B", "RIYPTNGYTRYADSVK"),
    "H3": ("B", "WGGDGFYAMDY"),
}
SEQ = {"A": LIGHT, "B": HEAVY}

# Locate each CDR in its chain -> inclusive auth residue range.
ranges = {}
for name, (chain, motif) in CDRS.items():
    idx = SEQ[chain].find(motif)
    if idx < 0:
        sys.exit(f"CDR {name} motif {motif} not found in chain {chain}")
    if SEQ[chain].find(motif, idx + 1) > 0:
        sys.exit(f"CDR {name} motif {motif} is not unique in chain {chain}")
    ranges[name] = (chain, idx + 1, idx + len(motif), motif)

print("=== CDR LOCATIONS (1N8Z auth numbering, sequential) ===")
for name in ("H1", "H2", "H3", "L1", "L2", "L3"):
    ch, lo, hi, motif = ranges[name]
    print(f"  CDR-{name}: chain {ch} residues {lo}-{hi}  len {len(motif):>2}  {motif}")


def cdr_of(chain, num):
    for name, (ch, lo, hi, _m) in ranges.items():
        if ch == chain and lo <= num <= hi:
            return name
    return "framework"


# Parse the contact pairs from the interface listing.
text = open(sys.argv[1] if len(sys.argv) > 1 else "interface_1n8z.txt",
            encoding="utf-8", errors="replace").read()

pair_re = re.compile(
    r"^\s+([A-Z])\s+([A-Z]{3})(\d+)\s+--\s+([A-Z])\s+([A-Z]{3})(\d+)\s+([\d.]+)\s*$",
    re.M)

ab_her2 = []  # (ab_chain, ab_num, ab_res, her2_num, her2_res, dist)
for m in pair_re.finditer(text):
    c1, r1, n1, c2, r2, n2, d = m.groups()
    chains = {c1, c2}
    if chains == {"A", "C"} or chains == {"B", "C"}:
        if c1 == "C":
            c1, r1, n1, c2, r2, n2 = c2, r2, n2, c1, r1, n1
        ab_her2.append((c1, int(n1), r1, int(n2), r2, float(d)))

if not ab_her2:
    sys.exit("no antibody-HER2 contact pairs parsed")

by_cdr = defaultdict(list)
for ch, num, res, hnum, hres, d in ab_her2:
    by_cdr[cdr_of(ch, num)].append((ch, num, res, hnum, hres, d))

total = len(ab_her2)
print(f"\n=== PER-CDR SHARE OF THE HER2 INTERFACE ({total} contact pairs) ===")
print(f"{'CDR':<10}{'pairs':>6}{'share':>8}{'residues':>10}  contacting residues")
order = sorted(by_cdr, key=lambda k: -len(by_cdr[k]))
for name in order:
    lst = by_cdr[name]
    res = sorted({(n, r) for _c, n, r, _hn, _hr, _d in lst})
    print(f"{name:<10}{len(lst):>6}{100*len(lst)/total:>7.0f}%{len(res):>10}  "
          + ", ".join(f"{r}{n}" for n, r in res))

heavy = sum(len(v) for k, v in by_cdr.items() if k.startswith("H"))
light = sum(len(v) for k, v in by_cdr.items() if k.startswith("L"))
print(f"\n  heavy-chain CDRs: {heavy} pairs ({100*heavy/total:.0f}%)")
print(f"  light-chain CDRs: {light} pairs ({100*light/total:.0f}%)")

# The design question: which CDRH3 positions are free of HER2 contact?
ch, lo, hi, motif = ranges["H3"]
contacting = {n for _c, n, _r, _hn, _hr, _d in by_cdr.get("H3", [])}
print(f"\n=== CDRH3 POSITION-BY-POSITION ({motif}, {lo}-{hi}) ===")
print(f"{'pos':>5} {'aa':>3}  {'HER2 contact?':<14} partners (min dist A)")
free = []
for i, aa in enumerate(motif):
    pos = lo + i
    hits = sorted({(hr, hn, d) for _c, n, _r, hn, hr, d in by_cdr.get("H3", [])
                   if n == pos}, key=lambda t: t[2])
    if hits:
        desc = ", ".join(f"{hr}{hn} {d:.2f}" for hr, hn, d in hits)
        print(f"{pos:>5} {aa:>3}  {'CONTACT':<14} {desc}")
    else:
        free.append((pos, aa))
        print(f"{pos:>5} {aa:>3}  {'free':<14} -")

n = len(motif)
need = n - int(n * 0.70)          # substitutions needed to fall below 70%
need = next(k for k in range(1, n + 1) if (n - k) / n < 0.70)
print(f"\n  CDRH3 length {n}; identity after k substitutions = (={n}-k)/{n}")
print(f"  substitutions needed for <70% identity: {need} "
      f"(gives {100*(n-need)/n:.1f}%)")
print(f"  HER2-contacting positions: {len(contacting)} -> "
      + ", ".join(f"{motif[p-lo]}{p}" for p in sorted(contacting)))
print(f"  contact-free positions:    {len(free)} -> "
      + ", ".join(f"{aa}{pos}" for pos, aa in free))
if len(free) >= need:
    print(f"\n  => <70% CDRH3 identity is reachable using contact-free "
          f"positions ALONE ({need} of {len(free)} available).")
else:
    print(f"\n  => contact-free positions are insufficient ({len(free)} < "
          f"{need}); at least {need-len(free)} contacting position(s) must change.")

# HER2 epitope span.
hres = sorted({(hn, hr) for _c, _n, _r, hn, hr, _d in ab_her2})
print(f"\n=== HER2 EPITOPE ({len(hres)} residues) ===")
print("  " + ", ".join(f"{r}{n}" for n, r in hres))
print(f"  span: {hres[0][0]}-{hres[-1][0]} "
      f"({hres[-1][0]-hres[0][0]+1} residues end to end, "
      f"{len(hres)} in actual contact -> discontinuous)")
