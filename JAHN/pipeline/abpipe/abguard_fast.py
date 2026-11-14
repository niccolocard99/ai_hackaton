"""Make abguard's novelty gate fast enough to screen a million candidates.

abguard grades CDRH3 novelty with a cross join in :func:`abguard.novelty._nearest`:
every distinct design CDRH3 against all 610,109 reference loops. Measured on the team
VM (4 polars threads, niced):

    numbering only            ~1,500 designs/s
    numbering + novelty          37 designs/s   <- 26 ms/design of cross join

So novelty is ~97% of the cost, and a 1M-candidate screen would need ~7.5 hours.

:func:`install` swaps that one function for :class:`abpipe.novelty_index.NoveltyIndex`,
which answers the same question with a deletion-neighbourhood index plus exact
Levenshtein verification: measured **233,578 designs/s**, a ~6,300x speedup on the
novelty stage, taking 1M candidates from 7.5 hours to under 10 seconds.

Nothing else about abguard changes. The numbering, the hard checks, the liability scan
and the grading logic are untouched, and ``tests/test_abguard_fast.py`` runs the real
abguard both ways to pin that every verdict, class and rejection reason is identical.

What the patch does and does not report
---------------------------------------
The index resolves distances up to a cap -- by default ``max(MUTATIONS_REQUIRED)``, so
3 -- and reports ``cap + 1`` for anything beyond it. Every gate comparison stays exact,
because the gate only ever asks ``cdrh3_min_distance >= required`` for a ``required``
of 2 or 3, and so does every rejection message, because a design that fails is within
the cap by construction. What you lose is the *true* distance of designs that comfortably
pass (reported as 4 rather than, say, 7) and their nearest-neighbour fields. That is
information the gate never uses. Set ``ABGUARD_FAST_MAX_DISTANCE`` higher if you want it,
at a cost that grows with the cap.

Usage::

    from abpipe import abguard_fast
    abguard_fast.install()      # idempotent; uninstall() restores the original

The reference index is cached per reference content, so repeated ``annotate`` calls in
one process build it once (~3.4 s for the bundled 610k loops at cap 1).
"""
from __future__ import annotations

import os

import polars as pl

from abpipe.novelty_index import NoveltyIndex

__all__ = ["install", "uninstall", "is_installed", "max_distance"]

_original = None
_cache: dict[tuple[int, int, int], NoveltyIndex] = {}


def max_distance() -> int:
    """The largest distance the patched gate resolves exactly.

    Defaults to the strictest rule abguard applies (``MUTATIONS_REQUIRED``), which is
    what keeps every verdict exact; ``ABGUARD_FAST_MAX_DISTANCE`` overrides it.
    """
    override = os.environ.get("ABGUARD_FAST_MAX_DISTANCE")
    if override:
        return int(override)
    from abguard.novelty import MUTATIONS_REQUIRED

    return max(MUTATIONS_REQUIRED.values())


def _index_for(reference: pl.DataFrame, cap: int) -> NoveltyIndex:
    """Build (or reuse) the index for this reference set."""
    loops = reference.get_column("cdrh3")
    digest = int(loops.hash().sum() or 0) if loops.len() else 0
    key = (loops.len(), digest, cap)
    index = _cache.get(key)
    if index is None:
        index = NoveltyIndex.build(loops.to_list(), max_distance=cap)
        _cache[key] = index
    return index


def _fast_nearest(reference, designs):
    """Drop-in for ``abguard.novelty._nearest``: same schema, same values, no cross join."""
    cap = max_distance()
    known = reference.select("cdrh3", "source", "id").collect().drop_nulls("cdrh3")
    queries = designs.select("cdrh3").drop_nulls().unique().collect()

    empty = pl.DataFrame(
        schema={
            "cdrh3": pl.String, "cdrh3_min_distance": pl.UInt32, "cdrh3_nearest": pl.String,
            "cdrh3_nearest_source": pl.String, "cdrh3_nearest_id": pl.String,
        }
    )
    if queries.height == 0 or known.height == 0:
        return empty.lazy()

    index = _index_for(known, cap)
    distance, reference_index = index.scan(queries.get_column("cdrh3").to_list())

    # Gather the matched reference row for the resolved queries; the rest stay null.
    sources = known.get_column("source").to_list()
    identifiers = known.get_column("id").to_list()
    loops = known.get_column("cdrh3").to_list()
    nearest, nearest_source, nearest_id = [], [], []
    for position in reference_index.tolist():
        if position < 0:
            nearest.append(None)
            nearest_source.append(None)
            nearest_id.append(None)
        else:
            nearest.append(loops[position])
            nearest_source.append(sources[position])
            nearest_id.append(identifiers[position])

    return pl.DataFrame(
        {
            "cdrh3": queries.get_column("cdrh3"),
            "cdrh3_min_distance": pl.Series(distance, dtype=pl.UInt32),
            "cdrh3_nearest": pl.Series(nearest, dtype=pl.String),
            "cdrh3_nearest_source": pl.Series(nearest_source, dtype=pl.String),
            "cdrh3_nearest_id": pl.Series(nearest_id, dtype=pl.String),
        }
    ).lazy()


def install() -> None:
    """Replace abguard's novelty cross join with the index. Idempotent."""
    global _original
    from abguard import novelty

    if _original is None:
        _original = novelty._nearest
    novelty._nearest = _fast_nearest


def uninstall() -> None:
    """Restore abguard's own implementation."""
    global _original
    if _original is None:
        return
    from abguard import novelty

    novelty._nearest = _original
    _original = None


def is_installed() -> bool:
    return _original is not None
