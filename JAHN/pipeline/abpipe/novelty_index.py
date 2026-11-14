"""Exact CDRH3 novelty distances without abguard's 610k-wide cross join.

abguard grades novelty in :func:`abguard.novelty._nearest` with
``queries.join(known, how="cross")``: every design against all 610,109 reference
loops. Measured on the team VM (4 polars threads, niced): **37 designs/s**, i.e.
26 ms/design of pure novelty work, so a 1M-candidate screen would take ~7.3 hours.
Numbering, by contrast, runs at ~1,500 designs/s -- the cross join is the whole cost.

The gate never needs the true nearest neighbour. It asks one question: *is the
closest reference closer than a small threshold?* That is answerable in time
independent of the reference size.

How it works
------------
1. **Deletion-neighbourhood index (SymSpell).** If ``levenshtein(q, r) <= d`` then the
   sets of subsequences obtainable by deleting up to ``d`` characters from ``q`` and
   from ``r`` must intersect. So an index over the references' deletion
   neighbourhoods turns the search into a handful of hash lookups.
2. **Verification.** The converse does not hold -- an intersection only proves
   ``distance <= 2d`` -- so every candidate is verified with a real Levenshtein
   computation. That is what makes the result *exact*, and it is also why hashing
   the keys is safe: a hash collision adds a candidate, and verification drops it.
3. **Length buckets, built lazily.** ``levenshtein(q, r) >= |len(q) - len(r)|``, so a
   query of length ``L`` can only be matched by references of length ``L-d .. L+d``.
   Buckets are indexed on first use, which keeps the reference's long tail (loops up
   to 60 residues, where the number of deletion combinations explodes) out of the
   index entirely when the queries are 13-mers.

The result is exact agreement with brute-force Levenshtein -- pinned by
``tests/test_novelty_index.py``, which compares against rapidfuzz over the real
reference set at thresholds 1, 2 and 3.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np

__all__ = ["NoveltyIndex"]

_FNV_OFFSET = np.uint64(14695981039346656037)
_FNV_PRIME = np.uint64(1099511628211)

#: Candidate (query, reference) pairs verified per batch. Bounds peak memory; the
#: verification itself is O(d * len) per pair in rapidfuzz's bit-parallel Levenshtein.
_VERIFY_CHUNK = 2_000_000


def _char_matrix(sequences: list[str], length: int) -> np.ndarray:
    """Equal-length sequences as an ``(n, length)`` uint8 matrix."""
    if not sequences:
        return np.empty((0, length), dtype=np.uint8)
    joined = "".join(sequences).encode("ascii")
    return np.frombuffer(joined, dtype=np.uint8).reshape(len(sequences), length)


def _hash_rows(matrix: np.ndarray, key_length: int) -> np.ndarray:
    """FNV-1a over the rows, with the key length folded in.

    Collisions are harmless: they only add candidates, which verification rejects.
    """
    out = np.full(matrix.shape[0], _FNV_OFFSET, dtype=np.uint64)
    with np.errstate(over="ignore"):
        for column in range(matrix.shape[1]):
            out = (out ^ matrix[:, column].astype(np.uint64)) * _FNV_PRIME
        out = (out ^ np.uint64(key_length)) * _FNV_PRIME
    return out


def _deletion_masks(length: int, max_deletions: int) -> list[tuple[int, ...]]:
    """Column selections left after deleting up to `max_deletions` positions."""
    masks: list[tuple[int, ...]] = []
    for deletions in range(min(max_deletions, length) + 1):
        for dropped in combinations(range(length), deletions):
            dropped_set = set(dropped)
            masks.append(tuple(p for p in range(length) if p not in dropped_set))
    return masks


def _neighbourhood_hashes(matrix: np.ndarray, max_deletions: int) -> tuple[np.ndarray, np.ndarray]:
    """Hashes of every up-to-`max_deletions` deletion subsequence of every row.

    Returns ``(hashes, row_index)``, both flat and aligned.
    """
    rows, length = matrix.shape
    if rows == 0:
        return np.empty(0, dtype=np.uint64), np.empty(0, dtype=np.int64)
    hashes: list[np.ndarray] = []
    for keep in _deletion_masks(length, max_deletions):
        sub = matrix[:, keep] if len(keep) != length else matrix
        hashes.append(_hash_rows(sub, len(keep)))
    stacked = np.concatenate(hashes)
    row_index = np.tile(np.arange(rows, dtype=np.int64), len(hashes))
    return stacked, row_index


@dataclass
class _Bucket:
    """References of one length, with their deletion index once it is built."""

    reference_ids: np.ndarray
    matrix: np.ndarray
    hashes: np.ndarray | None = None
    owners: np.ndarray | None = None

    def build(self, max_deletions: int) -> None:
        if self.hashes is not None:
            return
        hashes, rows = _neighbourhood_hashes(self.matrix, max_deletions)
        order = np.argsort(hashes, kind="stable")
        self.hashes = hashes[order]
        self.owners = self.reference_ids[rows[order]]


class NoveltyIndex:
    """Exact minimum Levenshtein distance to a reference set, capped at `max_distance`.

    Build once, query in batches::

        index = NoveltyIndex.build(reference_cdrh3s, max_distance=1)
        novel = index.min_distance(candidate_cdrh3s) >= 2   # the official rule

    `max_distance` is the largest distance the index can resolve. Distances above it
    saturate at ``max_distance + 1``, which is all a threshold gate needs and is why
    the index is cheap: cost grows with ``max_distance``, not with the reference size.
    """

    def __init__(self, references: list[str], max_distance: int) -> None:
        if max_distance < 0:
            raise ValueError("max_distance must be >= 0")
        self.references = references
        self.max_distance = max_distance
        self._buckets: dict[int, _Bucket] = {}
        by_length: dict[int, list[int]] = {}
        for i, sequence in enumerate(references):
            by_length.setdefault(len(sequence), []).append(i)
        for length, ids in by_length.items():
            identifiers = np.asarray(ids, dtype=np.int64)
            self._buckets[length] = _Bucket(
                reference_ids=identifiers,
                matrix=_char_matrix([references[i] for i in ids], length),
            )

    @classmethod
    def build(cls, references, max_distance: int) -> "NoveltyIndex":
        return cls(list(references), max_distance)

    # -- querying ------------------------------------------------------------------

    def min_distance(self, queries) -> np.ndarray:
        """Exact distance to the nearest reference, saturating at ``max_distance + 1``."""
        return self.scan(queries)[0]

    def nearest(self, queries) -> tuple[np.ndarray, list[str | None]]:
        """The nearest reference and its distance; ``None`` where the distance saturates.

        A design that *fails* the gate is by construction within ``max_distance`` of
        something, so its nearest reference is always resolved -- which is what a
        rejection message needs. Passing designs report ``None``: the index stopped
        looking once it knew they were far enough, which is the point.
        """
        distance, reference_index = self.scan(queries)
        nearest = [
            self.references[int(i)] if i >= 0 else None for i in reference_index
        ]
        return distance, nearest

    def scan(self, queries) -> tuple[np.ndarray, np.ndarray]:
        """``(distance, reference_index)``; reference_index is -1 where saturated.

        The index-level primitive: :meth:`min_distance` and :meth:`nearest` are views
        of it, and ``abpipe.abguard_fast`` uses it to rebuild abguard's nearest-neighbour
        columns without a second pass.
        """
        queries = list(queries)
        saturate = self.max_distance + 1
        # (distance, reference) packed into one integer so the reduction over candidate
        # pairs stays a single vectorised ufunc.at, with ties broken by reference order.
        stride = len(self.references) + 1
        packed = np.full(len(queries), saturate * stride, dtype=np.int64)
        if not queries:
            return packed // stride, np.full(len(queries), -1, dtype=np.int64)

        by_length: dict[int, list[int]] = {}
        for i, sequence in enumerate(queries):
            by_length.setdefault(len(sequence), []).append(i)

        for length, ids in by_length.items():
            positions = np.asarray(ids, dtype=np.int64)
            matrix = _char_matrix([queries[i] for i in ids], length)
            key_hashes, key_rows = _neighbourhood_hashes(matrix, self.max_distance)

            query_pairs: list[np.ndarray] = []
            reference_pairs: list[np.ndarray] = []
            for candidate_length in range(max(length - self.max_distance, 0), length + self.max_distance + 1):
                bucket = self._buckets.get(candidate_length)
                if bucket is None:
                    continue
                bucket.build(self.max_distance)
                lo = np.searchsorted(bucket.hashes, key_hashes, side="left")
                hi = np.searchsorted(bucket.hashes, key_hashes, side="right")
                counts = hi - lo
                total = int(counts.sum())
                if total == 0:
                    continue
                offsets = np.concatenate(([0], np.cumsum(counts)[:-1]))
                flat = np.arange(total) - np.repeat(offsets, counts) + np.repeat(lo, counts)
                query_pairs.append(positions[np.repeat(key_rows, counts)])
                reference_pairs.append(bucket.owners[flat])

            if not query_pairs:
                continue
            self._verify(np.concatenate(query_pairs), np.concatenate(reference_pairs), queries, packed, stride)

        distance = packed // stride
        reference_index = np.where(distance > self.max_distance, -1, packed % stride)
        return distance, reference_index

    def _verify(self, query_index: np.ndarray, reference_index: np.ndarray, queries,
                packed: np.ndarray, stride: int) -> None:
        """Replace candidate pairs with real Levenshtein distances and keep the minimum."""
        from rapidfuzz import process
        from rapidfuzz.distance import Levenshtein

        # One (query, reference) pair may be proposed by several deletion keys.
        combined = query_index.astype(np.uint64) * np.uint64(len(self.references)) + reference_index.astype(np.uint64)
        keep = np.unique(combined, return_index=True)[1]
        query_index = query_index[keep]
        reference_index = reference_index[keep]

        for start in range(0, query_index.size, _VERIFY_CHUNK):
            stop = start + _VERIFY_CHUNK
            left = [queries[i] for i in query_index[start:stop]]
            right = [self.references[i] for i in reference_index[start:stop]]
            distances = process.cpdist(
                left, right, scorer=Levenshtein.distance,
                score_cutoff=self.max_distance, workers=-1,
            ).astype(np.int64)
            candidate = distances * stride + reference_index[start:stop]
            # Above the cutoff rapidfuzz reports max_distance + 1; such a pair must not
            # claim a nearest reference, so it packs with the "no reference" slot.
            candidate = np.where(
                distances > self.max_distance, (self.max_distance + 1) * stride, candidate
            )
            np.minimum.at(packed, query_index[start:stop], candidate)

    def any_within(self, queries, distance: int) -> np.ndarray:
        """Whether some reference lies within `distance` edits of each query."""
        if distance > self.max_distance:
            raise ValueError(
                f"index resolves distances up to {self.max_distance}; rebuild with max_distance >= {distance}"
            )
        return self.min_distance(queries) <= distance
