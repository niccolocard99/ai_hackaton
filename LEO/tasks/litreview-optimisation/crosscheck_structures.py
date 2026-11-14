#!/usr/bin/env python3
"""Cross-check the trastuzumab paratope / HER2 epitope across several structures.

    python crosscheck_structures.py <dir-with-cifs> [cutoff]

Identifies, in each structure and without being told the chain ids:
  * trastuzumab heavy chain  - the chain whose sequence contains CDRH3 WGGDGFYAMDY
  * trastuzumab light chain  - the chain containing CDRL1 RASQDVNTAVA
  * HER2                     - the chain containing the 1N8Z epitope-region peptide
  * anything else            - reported as "other", by sequence, so a fused
                               domain (e.g. the 7PKL anti-idiotypic VHH, which is
                               deposited as part of the light-chain polymer) is
                               not silently mistaken for trastuzumab.

A trastuzumab light chain carrying a fused extra domain is split: the residue
range matching the trastuzumab light-chain sequence is treated as trastuzumab
and the remainder as a separate binding partner.

Then reports, per structure and in that structure's own author numbering:
  * per-CDR contact counts and contacting residues
  * the CDRH3 position-by-position contact map
  * the partner-side (HER2, or anti-idiotypic domain) contact residue set

Everything printed is MEASURED FROM THE DEPOSITED COORDINATES. Heavy atoms only,
first model, altloc A. No energies.
"""
import math
import os
import sys
from collections import defaultdict

from bsa_1n8z import load_atoms, CDR_MOTIFS

# Deposited trastuzumab variable-domain sequences (1N8Z chains A and B).
TRAS_L = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPS"
          "RFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
TRAS_H = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRY"
          "ADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")

H3 = "WGGDGFYAMDY"
L1 = "RASQDVNTAVA"
H1 = "DTYIH"

# Peptides unique to the HER2 ECD, used only to recognise which chain is HER2.
# Checked against the HER3 chains present in 7MN8/8YRY: none of these match.
HER2_TAGS = (
    "TQVCTGTDMKLRLPASPETHLDML",   # domain I, N-terminus
    "CPYNYLSTDVGSCTLVCP",         # domain III
    "PESFDGDPASNTAPLQPEQLQVFET",  # domain IV, the trastuzumab epitope region
)

# Offset of each anchor motif inside the trastuzumab chain, from 1N8Z: CDR-H1
# starts at heavy-chain residue 31 and CDR-L1 at light-chain residue 24. Used to
# find where a trastuzumab chain begins inside a longer fusion polymer (7PKL
# deposits the anti-idiotypic VHH and the light chain as one chain).
ANCHOR = {"H": (H1, TRAS_H.find(H1)), "L": (L1, TRAS_L.find(L1))}


def chain_seq(seqs, ch):
    """(sequence string, [author numbers]) for a chain, in numbering order."""
    nums = sorted(seqs[ch])
    return "".join(seqs[ch][n] for n in nums), nums


def classify(seqs):
    """{role: [(chain, [author numbers])]} plus a 'split' map for fusions."""
    roles = defaultdict(list)
    for ch in sorted(seqs):
        seq, nums = chain_seq(seqs, ch)
        if H3 in seq:
            roles["tras_H"].append((ch, nums, seq))
        elif L1 in seq:
            roles["tras_L"].append((ch, nums, seq))
        elif any(t in seq for t in HER2_TAGS):
            roles["HER2"].append((ch, nums, seq))
        else:
            roles["other"].append((ch, nums, seq))
    return roles


def tras_start_index(seq, logical):
    """Sequence index at which the trastuzumab chain itself begins.

    Found from a CDR anchor motif and that motif's known offset inside the
    trastuzumab chain, so it is unaffected by point mutations elsewhere and by
    author-numbering jumps. Returns 0 if the chain is not a fusion.
    """
    motif, off = ANCHOR[logical]
    i = seq.find(motif)
    if i < 0:
        return 0
    return max(0, i - off)


def contacts(atoms, groupA, groupB, cutoff):
    """{(chainA,numA,icA,resA,chainB,numB,icB,resB): (min_dist, [atom pairs])}."""
    A = [a for a in atoms if (a[0], a[1]) in groupA]
    B = [a for a in atoms if (a[0], a[1]) in groupB]
    if not A or not B:
        return {}
    cell = cutoff
    grid = defaultdict(list)
    for j, a in enumerate(B):
        x, y, z = a[6]
        grid[(int(x // cell), int(y // cell), int(z // cell))].append(j)
    out = {}
    c2 = cutoff * cutoff
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
                d = math.sqrt(d2)
                key = (a[0], a[1], a[2], a[3], b[0], b[1], b[2], b[3])
                rec = out.setdefault(key, [float("inf"), []])
                rec[1].append((a[4], a[5], b[4], b[5], d))
                if d < rec[0]:
                    rec[0] = d
    return out


def locate_cdrs_in(seq, nums):
    """{cdr: (lo, hi, [author numbers])} for whichever CDR motifs are present."""
    found = {}
    for name, (_logical, motif) in CDR_MOTIFS.items():
        i = seq.find(motif)
        if i < 0 or seq.find(motif, i + 1) > 0:
            continue
        auth = nums[i:i + len(motif)]
        found[name] = (auth[0], auth[-1], auth, motif)
    return found


def analyse(path, cutoff):
    name = os.path.basename(path).split(".")[0].upper()
    print("\n" + "=" * 78)
    print(f"=== {name}   cutoff {cutoff} A ===")
    print("=" * 78)
    atoms, seqs = load_atoms(path)
    roles = classify(seqs)

    for role in ("tras_H", "tras_L", "HER2", "other"):
        for ch, nums, seq in roles.get(role, []):
            print(f"  chain {ch:<3} {role:<7} {len(nums):>4} res  "
                  f"auth {nums[0]}..{nums[-1]}  {seq[:44]}...")

    if not roles.get("tras_H"):
        print("  !! no trastuzumab heavy chain found (CDRH3 motif absent) - skipped")
        return None

    # Trastuzumab residues, with any fused extra domain peeled off the N-terminus.
    tras_res, fused = set(), defaultdict(set)
    cdr_auth = {}
    for logical, role in (("H", "tras_H"), ("L", "tras_L")):
        for ch, nums, seq in roles.get(role, []):
            start = tras_start_index(seq, logical)
            tras_res |= {(ch, n) for n in nums[start:]}
            if start > 0:
                extra = nums[:start]
                label = f"{ch}:{extra[0]}-{extra[-1]} (FUSED non-trastuzumab domain)"
                fused[label] |= {(ch, n) for n in extra}
                print(f"  chain {ch} is a FUSION: author residues "
                      f"{extra[0]}-{extra[-1]} ({len(extra)} res) are NOT "
                      f"trastuzumab; trastuzumab {logical} chain starts at "
                      f"author {nums[start]}")
            cdr_auth.update({k: (ch,) + v for k, v in
                             locate_cdrs_in(seq, nums).items()
                             if k.startswith(logical)})

    # Candidate partners: every HER2 chain, every 'other' chain, plus any fused
    # domain carved out of a trastuzumab chain.
    partners = {}
    for ch, nums, _seq in roles.get("HER2", []):
        partners[f"{ch} (HER2)"] = {(ch, n) for n in nums}
    for ch, nums, _seq in roles.get("other", []):
        partners[f"{ch} (other)"] = {(ch, n) for n in nums}
    partners.update(fused)

    results = {}
    for label, pset in sorted(partners.items()):
        cmap = contacts(atoms, tras_res, pset, cutoff)
        if not cmap:
            continue
        results[label] = cmap
        pairs = len(cmap)
        ab = sorted({(k[0], k[1], k[3]) for k in cmap})
        pt = sorted({(k[5], k[7]) for k in cmap})
        print(f"\n  --- trastuzumab <-> {label}: {pairs} residue pairs, "
              f"{len(ab)} paratope res, {len(pt)} partner res ---")
        print(f"      partner contact residues ({len(pt)}): "
              + ", ".join(f"{r}{n}" for n, r in pt))
        if pt:
            print(f"      partner span: {pt[0][0]}-{pt[-1][0]}")

        # per-CDR breakdown
        def cdr_of(ch, num):
            for cname, (c, _lo, _hi, auth, _m) in cdr_auth.items():
                if c == ch and num in auth:
                    return cname
            return "framework"
        by = defaultdict(list)
        for k, (d, _al) in cmap.items():
            by[cdr_of(k[0], k[1])].append((k, d))
        print(f"      {'CDR':<11}{'pairs':>6}{'share':>8}  contacting residues")
        for cname in sorted(by, key=lambda x: -len(by[x])):
            res = sorted({(k[1], k[3]) for k, _d in by[cname]})
            print(f"      {cname:<11}{len(by[cname]):>6}"
                  f"{100*len(by[cname])/pairs:>7.0f}%  "
                  + ", ".join(f"{r}{n}" for n, r in res))

        # CDRH3 position map
        if "H3" in cdr_auth:
            ch3, lo, hi, auth, motif = cdr_auth["H3"]
            print(f"      CDRH3 {motif} (chain {ch3} {lo}-{hi}) vs {label}:")
            for aa, num in zip(motif, auth):
                hits = sorted({(k[7], k[5], d) for k, (d, _a) in cmap.items()
                               if k[0] == ch3 and k[1] == num}, key=lambda t: t[2])
                if hits:
                    print(f"        {num:>5} {aa}  CONTACT   "
                          + ", ".join(f"{r}{n} {d:.2f}" for r, n, d in hits))
                else:
                    print(f"        {num:>5} {aa}  free")
    if not results:
        print("  (trastuzumab makes no inter-chain contact in this structure)")
    return results


def main():
    d = sys.argv[1]
    cutoff = float(sys.argv[2]) if len(sys.argv) > 2 else 4.5
    files = sorted(f for f in os.listdir(d) if f.lower().endswith(".cif"))
    if not files:
        sys.exit(f"no .cif files in {d}")
    for f in files:
        try:
            analyse(os.path.join(d, f), cutoff)
        except Exception as exc:          # keep going; report the failure
            print(f"\n!! {f} FAILED: {type(exc).__name__}: {exc}")


if __name__ == "__main__":
    main()
