#!/usr/bin/env python3
"""Consensus trastuzumab paratope / HER2 epitope across every available structure.

    python epitope_consensus.py <1n8z.cif> <dir-with-other-cifs> [cutoff]

Two things make a naive cross-structure comparison wrong, and this script fixes
both:

1. **HER2 numbering differs between depositions.** Some use mature-protein
   numbering (1N8Z: 1-607), others UniProt P04626 precursor numbering (+22).
   The offset is determined here by sequence alignment of each HER2 chain
   against 1N8Z chain C, and is reported with its identity so the mapping is
   auditable, not assumed. Everything is then reported in 1N8Z mature numbering.

2. **"Not a contact" and "not modelled" are different.** A disordered loop is
   absent from the coordinates, so it can never be scored as a contact. Each
   cell of the output therefore distinguishes:
       X  = contact at the cutoff
       .  = modelled, but no contact
       ?  = NOT MODELLED in this structure, so uninformative

Only cells marked X or . are evidence. Disagreements between X and . are real;
a ? is not a disagreement.
"""
import os
import sys
from collections import defaultdict

from bsa_1n8z import load_atoms
from crosscheck_structures import (classify, contacts, locate_cdrs_in,
                                   tras_start_index, HER2_TAGS)

H3_MOTIF = "WGGDGFYAMDY"


def her2_offset(S, R):
    """Integer k such that this chain's residue n corresponds to 1N8Z's n-k."""
    best = None
    for k in range(-80, 81):
        common = [n for n in S if (n - k) in R]
        if len(common) < 50:
            continue
        m = sum(1 for n in common if S[n] == R[n - k])
        if best is None or m > best[1]:
            best = (k, m, len(common))
    return best


def analyse(path, refC, cutoff):
    """-> (name, offset, identity, epitope_set_1n8z_numbering, modelled_set,
            h3_contacts, h3_partners, pairs, her2_chain)"""
    name = os.path.basename(path).split(".")[0].upper()
    atoms, seqs = load_atoms(path)
    roles = classify(seqs)
    if not roles.get("tras_H"):
        return None

    tras_res, cdr_auth = set(), {}
    for logical, role in (("H", "tras_H"), ("L", "tras_L")):
        for ch, nums, seq in roles.get(role, []):
            start = tras_start_index(seq, logical)
            tras_res |= {(ch, n) for n in nums[start:]}
            cdr_auth.update({k: (ch,) + v for k, v in
                             locate_cdrs_in(seq, nums).items()
                             if k.startswith(logical)})

    # Pick the HER2 chain that trastuzumab actually touches most.
    outs = []
    for ch, nums, _seq in roles.get("HER2", []):
        off = her2_offset(seqs[ch], refC)
        if off is None:
            continue
        k, m, tot = off
        cmap = contacts(atoms, tras_res, {(ch, n) for n in nums}, cutoff)
        outs.append((len(cmap), ch, nums, k, m, tot, cmap))
    if not outs:
        return None
    outs.sort(reverse=True)
    npairs, ch, nums, k, m, tot, cmap = outs[0]
    extra = [(o[1], o[0]) for o in outs[1:] if o[0] > 0]

    epi = {n - k for (_ac, _an, _ai, _ar, _hc, n, _ic, _r) in cmap}
    modelled = {n - k for n in nums}
    h3 = {}
    if "H3" in cdr_auth:
        hch, lo, hi, auth, motif = cdr_auth["H3"]
        for aa, num in zip(motif, auth):
            hits = sorted(((kk[5] - k, kk[7], d) for kk, (d, _a) in cmap.items()
                           if kk[0] == hch and kk[1] == num), key=lambda t: t[2])
            h3[num] = hits
        h3_present = set(auth)
    else:
        h3_present = set()

    per_cdr = defaultdict(int)

    def cdr_of(c, n):
        for cn, (cc, _l, _h, au, _mm) in cdr_auth.items():
            if cc == c and n in au:
                return cn
        return "framework"
    for kk in cmap:
        per_cdr[cdr_of(kk[0], kk[1])] += 1

    return dict(name=name, chain=ch, offset=k, ident=(m, tot), epi=epi,
                modelled=modelled, h3=h3, h3_present=h3_present,
                pairs=npairs, per_cdr=dict(per_cdr), extra=extra,
                cdr_auth=cdr_auth)


def main():
    ref_path, other_dir = sys.argv[1], sys.argv[2]
    cutoff = float(sys.argv[3]) if len(sys.argv) > 3 else 4.5

    _a, refseqs = load_atoms(ref_path)
    refC = refseqs["C"]

    results = [analyse(ref_path, refC, cutoff)]
    for f in sorted(os.listdir(other_dir)):
        if not f.lower().endswith(".cif"):
            continue
        try:
            r = analyse(os.path.join(other_dir, f), refC, cutoff)
        except Exception as exc:
            print(f"!! {f} FAILED: {type(exc).__name__}: {exc}")
            continue
        if r is None:
            print(f"-- {f}: no trastuzumab/HER2 pair, skipped")
            continue
        results.append(r)

    print(f"\n=== STRUCTURES USED (cutoff {cutoff} A, heavy atoms) ===")
    print(f"  {'PDB':<6}{'HER2 ch':>8}{'numbering':>12}{'seq identity to 1N8Z':>23}"
          f"{'pairs':>7}  HER2 modelled range (1N8Z numbering)")
    for r in results:
        mn, mx = min(r["modelled"]), max(r["modelled"])
        nm = sorted(set(range(mn, mx + 1)) - r["modelled"])
        gaps = []
        if nm:
            s = nm[0]
            for a, b in zip(nm, nm[1:] + [None]):
                if b != (a + 1 if a is not None else None):
                    gaps.append(f"{s}-{a}" if s != a else f"{a}")
                    s = b
        idtxt = "%d/%d" % r["ident"]
        print(f"  {r['name']:<6}{r['chain']:>8}{'1N8Z+' + str(r['offset']):>12}"
              f"{idtxt:>23}"
              f"{r['pairs']:>7}  {mn}-{mx}"
              + (f"   GAPS: {', '.join(gaps)}" if gaps else "   (no gaps)"))
        if r["extra"]:
            print(f"         (also contacts a second HER2 chain: "
                  + ", ".join(f"{c}:{n} pairs" for c, n in r["extra"]) + ")")

    # ---------------- epitope consensus ----------------
    allres = sorted(set().union(*[r["epi"] for r in results]))
    print(f"\n=== HER2 EPITOPE CONSENSUS ({len(allres)} residues contacted in "
          f">=1 structure), 1N8Z mature numbering ===")
    print("  X = contact   . = modelled, no contact   ? = NOT MODELLED "
          "(uninformative)\n")
    hdr = "  res    aa  " + " ".join(f"{r['name']:>5}" for r in results) + "   n_X/n_informative"
    print(hdr)
    agree, disagree, partial = [], [], []
    for n in allres:
        aa = refC.get(n, "?")
        cells, nx, ninf = [], 0, 0
        for r in results:
            if n not in r["modelled"]:
                cells.append("?")
            elif n in r["epi"]:
                cells.append("X")
                nx += 1
                ninf += 1
            else:
                cells.append(".")
                ninf += 1
        print(f"  {n:<6} {aa}   " + " ".join(f"{c:>5}" for c in cells)
              + f"     {nx}/{ninf}")
        if ninf and nx == ninf:
            agree.append(n)
        elif nx == 1 and ninf > 2:
            disagree.append(n)
        else:
            partial.append(n)
    print(f"\n  CORE EPITOPE (contact in EVERY structure that models it, "
          f"{len(agree)} residues):")
    print("    " + ", ".join(f"{refC.get(n,'?')}{n}" for n in agree))
    print(f"  PERIPHERAL (contact in some, not others, {len(partial)}):")
    print("    " + ", ".join(f"{refC.get(n,'?')}{n}" for n in partial))
    print(f"  SINGLE-STRUCTURE ONLY (contact in exactly one of >2 informative, "
          f"{len(disagree)}):")
    print("    " + ", ".join(f"{refC.get(n,'?')}{n}" for n in disagree))

    # ---------------- CDRH3 consensus ----------------
    print(f"\n=== CDRH3 {H3_MOTIF} POSITION CONSENSUS (trastuzumab numbering "
          f"99-109, identical in every structure) ===")
    print("  X = contacts HER2   . = modelled, no HER2 contact   "
          "? = position not modelled\n")
    print("  pos aa  " + " ".join(f"{r['name']:>5}" for r in results)
          + "   n_X/n_informative   HER2 partners (union, 1N8Z numbering)")
    for i, aa in enumerate(H3_MOTIF):
        pos = 99 + i
        cells, nx, ninf, partners = [], 0, 0, {}
        for r in results:
            if pos not in r["h3_present"]:
                cells.append("?")
                continue
            ninf += 1
            hits = r["h3"].get(pos, [])
            if hits:
                cells.append("X")
                nx += 1
                for hn, hr, d in hits:
                    if hn not in partners or d < partners[hn][1]:
                        partners[hn] = (hr, d)
            else:
                cells.append(".")
        pl = ", ".join(f"{partners[hn][0]}{hn}" for hn in sorted(partners))
        print(f"  {pos:<4}{aa}  " + " ".join(f"{c:>5}" for c in cells)
              + f"     {nx}/{ninf}          {pl}")

    print(f"\n  CDRH3 positions contacting HER2 in EVERY structure that models them:")
    always = [99 + i for i, _ in enumerate(H3_MOTIF)
              if all((99 + i) not in r["h3_present"] or r["h3"].get(99 + i)
                     for r in results)
              and any((99 + i) in r["h3_present"] for r in results)]
    print("    " + ", ".join(f"{H3_MOTIF[p-99]}{p}" for p in always))
    never = [99 + i for i, _ in enumerate(H3_MOTIF)
             if all((99 + i) not in r["h3_present"] or not r["h3"].get(99 + i)
                    for r in results)]
    print(f"  CDRH3 positions contacting HER2 in NO structure "
          f"(the only safe free set):")
    print("    " + ", ".join(f"{H3_MOTIF[p-99]}{p}" for p in never))
    some = [99 + i for i, _ in enumerate(H3_MOTIF)
            if (99 + i) not in always and (99 + i) not in never]
    print(f"  CDRH3 positions contacting HER2 in SOME structures only "
          f"(1N8Z alone would call these free or occupied - the disagreements):")
    print("    " + ", ".join(f"{H3_MOTIF[p-99]}{p}" for p in some))

    print(f"\n=== PER-CDR CONTACT SHARE BY STRUCTURE ===")
    names = ["H1", "H2", "H3", "L1", "L2", "L3", "framework"]
    print(f"  {'PDB':<7}" + "".join(f"{n:>11}" for n in names))
    for r in results:
        tot = r["pairs"] or 1
        row = ""
        for n in names:
            v = r["per_cdr"].get(n, 0)
            row += f"{v:>4} {100*v/tot:>4.0f}%" if v else f"{'-':>10} "
        notfound = [n for n in ("H1", "H2", "H3", "L1", "L2", "L3")
                    if n not in r["cdr_auth"]]
        print(f"  {r['name']:<7}" + row
              + ("   [CDR motif mutated/unmodelled, not located: "
                 + ",".join(notfound) + "]" if notfound else ""))


if __name__ == "__main__":
    main()
