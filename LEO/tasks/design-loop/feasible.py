#!/usr/bin/env python3
"""How big is the region a legal optimisation design can live in?

Three constraints act on the designed 10-mer at once, and the ledger has only
ever costed the first:

  identity      >= 4 substitutions from trastuzumab  (13-mer CDRH3 < 70%)
  novelty       >= 3 mutations from all 525,063 reference CDRH3s (abguard)
  ranker scope  <= 4 mutations from a measured SPR binder, or pred_log10_kd
                is uninformative (AGENT_GUIDE.md section 3: Spearman -0.057
                at Hamming >= 5)

The second and third pull in opposite directions: novelty pushes designs away
from everything measured, the ranker's validity pulls them back. This script
measures whether a region satisfying all three exists, and checks what the
constraints do to the 21 variants EXP-008 called "legal".

Writes nothing; prints. Run from the task folder.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from her2lib import TZ_LOOP  # noqa: E402
from novelty import NoveltyIndex, encode, load_reference, min_dist_to  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA, REF = HERE / "data", HERE / "ref"


def main() -> int:
    print("building the reference index ...")
    ref_codes, counts = load_reference(DATA, REF)
    for k, v in counts.items():
        print(f"  {k:<16} {v:>9,}")
    idx = NoveltyIndex(ref_codes)
    print(f"  masks            {len(idx.masks):>9}  (1 + 10 + 45)")

    tz = encode(TZ_LOOP)
    print(f"\ntrastuzumab's own loop {TZ_LOOP}: min distance to reference = "
          f"{idx.min_distance_capped(tz)[0]} (expected 0, it is in the tables)")

    # ---- the measured binders, and the band EXP-008 called legal -------------
    with (DATA / "trastuzumab_her2_spr_all_768.csv").open(encoding="utf-8") as fh:
        spr = [r for r in csv.DictReader(fh)
               if len(r["cdrh3"]) == 10 and set(r["cdrh3"]) <= set("ACDEFGHIKLMNPQRSTVWY")]
    binders = [r for r in spr if r["binding_result"] == "quantified_binding"]
    b_loops = [r["cdrh3"] for r in binders]
    b_kd = np.array([float(r["kd_nM"]) for r in binders])
    b_codes = encode(b_loops)
    d_tz = (b_codes != tz[0]).sum(1)
    print(f"\nmeasured binders with a KD: {len(binders)}")

    legal = d_tz >= 4
    print(f"  of those, >= 4 substitutions from trastuzumab (EXP-008's legal band): {legal.sum()}")
    nd = idx.min_distance_capped(b_codes[legal])
    print(f"  ... and how many of those pass abguard's novelty gate (>= 3)?  {(nd == 3).sum()}")
    print("  -> every measured variant is itself in the reference tables, so the")
    print("     'legal band' EXP-008 ranked is entirely unsubmittable. It is a")
    print("     source of design rules, never of candidates.")

    tight = b_kd < 0.152
    print(f"\n  measured binders tighter than trastuzumab (0.152 nM): {tight.sum()}")
    print(f"  of those, >= 4 substitutions from trastuzumab:        {(tight & legal).sum()}")

    # ---- does a feasible region exist at all? -------------------------------
    # Walk outward from the measured binders: all loops at exactly Hamming r
    # from a seed, restricted to residues seen at that position among binders.
    print("\n--- is there anywhere legal to design? ---")
    obs = [sorted({b_loops[i][p] for i in range(len(b_loops))}) for p in range(10)]
    print("residues observed per position among measured binders:")
    for p in range(10):
        print(f"  pos {99 + p:>3} {TZ_LOOP[p]} -> {''.join(obs[p]):<21} ({len(obs[p])})")

    rng = np.random.default_rng(0)
    print(f"\n{'from':<6} {'r':>2} {'sampled':>9} {'>=4 from tz':>12} {'novel':>8} "
          f"{'both':>8} {'yield':>7}")
    for r in (3, 4, 5):
        cand = set()
        target = 60000
        while len(cand) < target:
            s = rng.integers(0, len(b_loops), 4000)
            pos = np.array([rng.permutation(10)[:r] for _ in range(4000)])
            for si, ps in zip(s, pos):
                lp = list(b_loops[si])
                for p in ps:
                    alt = [a for a in obs[p] if a != lp[p]]
                    if alt:
                        lp[p] = alt[rng.integers(0, len(alt))]
                cand.add("".join(lp))
            if len(cand) > target * 3:
                break
        cl = sorted(cand)
        cc = encode(cl)
        d_t = (cc != tz[0]).sum(1)
        ok_id = d_t >= 4
        nov = idx.min_distance_capped(cc) == 3
        both = ok_id & nov
        print(f"{'binder':<6} {r:>2} {len(cl):>9,} {ok_id.sum():>12,} {nov.sum():>8,} "
              f"{both.sum():>8,} {both.mean():>6.1%}")

        if both.any():
            d_seed = min_dist_to(cc[both], b_codes)
            print(f"       -> of the {both.sum():,} legal ones, min distance to a measured "
                  f"binder: {np.bincount(d_seed, minlength=6)[:6].tolist()} (index = distance)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
