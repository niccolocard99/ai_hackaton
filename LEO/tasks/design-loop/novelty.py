"""Exact novelty index over the 525,063 reference CDRH3s.

abguard's binding rule for the optimisation track is `cdrh3_min_distance >= 3`
against its whole reference set, not just against trastuzumab. The probe in
probe_abguard.py showed that is the constraint that actually bites: a design 6
substitutions from trastuzumab still sat 2 mutations from a FACS library member
and was rejected. So we need to answer "is this loop >= 3 from everything
known?" for millions of candidates, exactly, and fast.

Method. For fixed length L, two loops are within Hamming <= k iff there is a set
of k positions outside which they agree. So for k = 2 we build, for each of the
C(10,2) = 45 position pairs, the set of reference loops with those two positions
masked out. A candidate is within 2 of some reference iff one of its 45 masked
keys is present in the corresponding set. The same trick at k = 1 (10 masks) and
k = 0 (1 mask) gives the exact minimum distance whenever it is below 3, which is
all the gate needs.

Keys are base-20 integers over the unmasked positions, so a key fits in an
int64 and membership is a sorted-array `searchsorted`, fully vectorised.
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"
L = 10
_CODE = np.full(256, -1, np.int8)
for _i, _c in enumerate(AA):
    _CODE[ord(_c)] = _i


def encode(seqs) -> np.ndarray:
    """(N, L) int8 residue codes. Raises on a non-standard residue."""
    if isinstance(seqs, str):
        seqs = [seqs]
    buf = np.frombuffer("".join(seqs).encode(), np.uint8)
    if buf.size % L:
        raise ValueError(f"not all sequences are {L} aa")
    out = _CODE[buf].reshape(-1, L)
    if out.min() < 0:
        raise ValueError("non-standard residue in input")
    return out


def _key(codes: np.ndarray, keep: tuple[int, ...]) -> np.ndarray:
    """Base-20 integer over the kept positions."""
    k = np.zeros(len(codes), np.int64)
    for p in keep:
        k = k * 20 + codes[:, p]
    return k


class NoveltyIndex:
    """Exact 'is this loop within Hamming <= 2 of any reference loop' test."""

    def __init__(self, ref_codes: np.ndarray):
        self.n_ref = len(ref_codes)
        self.masks: list[tuple[int, ...]] = [()]
        self.masks += [(i,) for i in range(L)]
        self.masks += list(combinations(range(L), 2))
        self._tables = []
        for m in self.masks:
            keep = tuple(p for p in range(L) if p not in m)
            self._tables.append(np.unique(_key(ref_codes, keep)))
        self._keeps = [tuple(p for p in range(L) if p not in m) for m in self.masks]

    def _hit(self, codes: np.ndarray, lo: int, hi: int) -> np.ndarray:
        """Does any mask in masks[lo:hi] match a reference? -> (N,) bool"""
        out = np.zeros(len(codes), bool)
        for t in range(lo, hi):
            tab = self._tables[t]
            k = _key(codes, self._keeps[t])
            pos = np.searchsorted(tab, k)
            np.clip(pos, 0, len(tab) - 1, out=pos)
            out |= tab[pos] == k
        return out

    def min_distance_capped(self, codes: np.ndarray) -> np.ndarray:
        """Exact min Hamming distance to the reference set, capped at 3.

        A returned 3 means ">= 3", which is exactly what abguard's gate asks.
        """
        n = len(codes)
        out = np.full(n, 3, np.int8)
        d0 = self._hit(codes, 0, 1)                       # exact match
        out[d0] = 0
        rest = ~d0
        if rest.any():
            d1 = self._hit(codes[rest], 1, 1 + L)         # within 1
            idx = np.flatnonzero(rest)
            out[idx[d1]] = 1
            rest2 = idx[~d1]
            if rest2.size:
                d2 = self._hit(codes[rest2], 1 + L, len(self.masks))  # within 2
                out[rest2[d2]] = 2
        return out

    def is_novel(self, codes: np.ndarray) -> np.ndarray:
        """abguard's gate: min distance to every reference CDRH3 >= 3."""
        return self.min_distance_capped(codes) == 3


def load_reference(data_dir: Path, ref_dir: Path) -> tuple[np.ndarray, dict]:
    """Every 10-mer CDRH3 abguard compares against, as residue codes.

    Two sources, the same two abguard names in its rejection messages:
    `trastuzumab_her2_spr_all_768` and `tz_her2_facs_3classes_524346`.
    """
    import csv

    loops: list[str] = []
    counts = {}
    spr = data_dir / "trastuzumab_her2_spr_all_768.csv"
    with spr.open(encoding="utf-8") as fh:
        got = [r["cdrh3"] for r in csv.DictReader(fh)]
    keep = [s for s in got if len(s) == L and set(s) <= set(AA)]
    counts["spr_768"] = len(keep)
    loops += keep

    facs = ref_dir / "tz_her2_facs_3classes_524346.csv"
    with facs.open(encoding="utf-8") as fh:
        rd = csv.reader(fh)
        next(rd)
        got = [row[0] for row in rd]
    keep = [s for s in got if len(s) == L and set(s) <= set(AA)]
    counts["facs_524346"] = len(keep)
    loops += keep

    codes = encode(loops)
    uniq = np.unique(_key(codes, tuple(range(L))))
    counts["total_rows"] = len(loops)
    counts["unique_loops"] = len(uniq)
    return codes, counts


def min_dist_to(codes: np.ndarray, seed_codes: np.ndarray, batch: int = 20000) -> np.ndarray:
    """Exact min Hamming distance from each candidate to a small seed set."""
    out = np.empty(len(codes), np.int8)
    for i in range(0, len(codes), batch):
        blk = codes[i : i + batch]
        d = (blk[:, None, :] != seed_codes[None, :, :]).sum(2)
        out[i : i + batch] = d.min(1)
    return out
