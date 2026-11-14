"""The fast novelty index must agree with brute-force Levenshtein, exactly.

abguard grades CDRH3 novelty with a cross join: every design against all 610,109
reference loops (measured 26 ms/design, so ~7.3 h for 1M candidates). The gate only
ever asks whether the nearest reference is closer than a small threshold, so a
deletion-neighbourhood index answers it in O(1) per query instead of O(610k).

That is only worth anything if the verdicts are identical. These tests pin that:
the hand-written cases fix the semantics, and `test_matches_bruteforce_*` compares
against rapidfuzz over the real reference set, including deliberately
boundary-adjacent queries.
"""
from __future__ import annotations

import os
import pathlib
import random

import numpy as np
import pytest

from abpipe.novelty_index import NoveltyIndex

REF_PARQUET_ENV = "ABGUARD_CDRH3_PARQUET"


def _reference_path() -> pathlib.Path | None:
    candidates = []
    if os.environ.get(REF_PARQUET_ENV):
        candidates.append(pathlib.Path(os.environ[REF_PARQUET_ENV]))
    job = os.environ.get("CLAUDE_JOB_DIR")
    if job:
        candidates.append(pathlib.Path(job) / "tmp/refdata/abguard_cdrh3.parquet")
    for pattern in ("~/.venvs/abguard/lib/python*/site-packages/abguard/reference_data/abguard_cdrh3.parquet",):
        candidates.extend(sorted(pathlib.Path("/").glob(os.path.expanduser(pattern).lstrip("/"))))
    for path in candidates:
        if path.exists():
            return path
    return None


@pytest.fixture(scope="module")
def reference() -> list[str]:
    path = _reference_path()
    if path is None:
        pytest.skip(f"no abguard CDRH3 reference parquet (set {REF_PARQUET_ENV})")
    import polars as pl

    return pl.read_parquet(path).get_column("cdrh3").to_list()


# --- semantics -------------------------------------------------------------------


def test_identical_sequence_is_within_distance_zero():
    index = NoveltyIndex.build(["SRWGGDGFYAMDY"], max_distance=1)
    assert index.any_within(["SRWGGDGFYAMDY"], 0).tolist() == [True]


def test_single_substitution_is_within_one_but_not_zero():
    index = NoveltyIndex.build(["SRWGGDGFYAMDY"], max_distance=1)
    query = ["SRWGGDGFYAMDA"]  # final Y -> A
    assert index.any_within(query, 1).tolist() == [True]
    assert index.any_within(query, 0).tolist() == [False]


def test_insertion_is_within_one():
    index = NoveltyIndex.build(["SRWGGDGFYAMDY"], max_distance=1)
    assert index.any_within(["SRWGGDGFYAMDYK"], 1).tolist() == [True]


def test_deletion_is_within_one():
    index = NoveltyIndex.build(["SRWGGDGFYAMDY"], max_distance=1)
    assert index.any_within(["SRWGGDGFYAMD"], 1).tolist() == [True]


def test_unrelated_sequence_is_not_within_threshold():
    index = NoveltyIndex.build(["SRWGGDGFYAMDY"], max_distance=2)
    assert index.any_within(["AAAAAAAAAAAAA"], 2).tolist() == [False]


def test_min_distance_is_exact_below_cap_and_saturates_above():
    index = NoveltyIndex.build(["SRWGGDGFYAMDY"], max_distance=2)
    got = index.min_distance(["SRWGGDGFYAMDY", "SRWGGDGFYAMDA", "SRWGGDGFYAMAA", "AAAAAAAAAAAAA"])
    assert got.tolist() == [0, 1, 2, 3]  # 3 is the saturation sentinel (max_distance + 1)


def test_nearest_reference_differing_only_outside_the_window_is_found():
    index = NoveltyIndex.build(["AAAAAAAAAAAAA", "SRWGGDGFYAMDY"], max_distance=1)
    assert index.any_within(["SRWGGDGFYAMDY"], 0).tolist() == [True]


# --- the property that matters ---------------------------------------------------


def _bruteforce_min_distance(queries: list[str], references: list[str]) -> np.ndarray:
    from rapidfuzz import process
    from rapidfuzz.distance import Levenshtein

    scores = process.cdist(queries, references, scorer=Levenshtein.distance, workers=-1)
    return scores.min(axis=1)


def _boundary_queries(references: list[str], n: int, seed: int) -> list[str]:
    """Queries that straddle the threshold: real loops mutated 0-5 edits, plus noise."""
    rng = random.Random(seed)
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    picks = [references[rng.randrange(len(references))] for _ in range(n)]
    out: list[str] = []
    for i, seq in enumerate(picks):
        edits = i % 6
        chars = list(seq)
        for _ in range(edits):
            kind = rng.choice(("sub", "ins", "del"))
            if kind == "sub" and chars:
                chars[rng.randrange(len(chars))] = rng.choice(alphabet)
            elif kind == "ins":
                chars.insert(rng.randrange(len(chars) + 1), rng.choice(alphabet))
            elif len(chars) > 1:
                del chars[rng.randrange(len(chars))]
        out.append("".join(chars))
    out.extend("".join(rng.choice(alphabet) for _ in range(13)) for _ in range(n // 4))
    return out


@pytest.mark.parametrize("max_distance", [1, 2, 3])
def test_matches_bruteforce_on_real_reference(reference, max_distance):
    queries = _boundary_queries(reference, 160, seed=20261114 + max_distance)
    index = NoveltyIndex.build(reference, max_distance=max_distance)

    expected_min = _bruteforce_min_distance(queries, reference)
    got_min = index.min_distance(queries)

    saturated = expected_min > max_distance
    assert np.array_equal(got_min[~saturated], expected_min[~saturated]), "exact distances disagree"
    assert np.all(got_min[saturated] == max_distance + 1), "saturated rows must report the sentinel"

    for threshold in range(max_distance + 1):
        assert np.array_equal(
            index.any_within(queries, threshold), expected_min <= threshold
        ), f"verdict disagrees at threshold {threshold}"


def test_nearest_reports_which_reference_matched():
    index = NoveltyIndex.build(["AAAAAAAAAAAAA", "SRWGGDGFYAMDY"], max_distance=2)
    distance, nearest = index.nearest(["SRWGGDGFYAMDA"])
    assert distance.tolist() == [1]
    assert nearest == ["SRWGGDGFYAMDY"]


def test_nearest_is_none_beyond_the_cap():
    index = NoveltyIndex.build(["SRWGGDGFYAMDY"], max_distance=1)
    distance, nearest = index.nearest(["AAAAAAAAAAAAA"])
    assert distance.tolist() == [2]  # the saturation sentinel
    assert nearest == [None]


def test_nearest_reference_really_is_at_the_reported_distance(reference):
    from rapidfuzz.distance import Levenshtein

    queries = _boundary_queries(reference, 120, seed=99)
    index = NoveltyIndex.build(reference, max_distance=3)
    distance, nearest = index.nearest(queries)

    for query, dist, match in zip(queries, distance.tolist(), nearest):
        if match is None:
            assert dist == 4, "a null nearest must mean saturation"
        else:
            assert Levenshtein.distance(query, match) == dist


def test_official_gate_verdict_matches_bruteforce(reference):
    """The real rule: a design is novel iff no reference is within 1 edit (min >= 2)."""
    queries = _boundary_queries(reference, 200, seed=7)
    index = NoveltyIndex.build(reference, max_distance=1)

    expected_novel = _bruteforce_min_distance(queries, reference) >= 2
    got_novel = index.min_distance(queries) >= 2

    assert np.array_equal(got_novel, expected_novel)
    assert expected_novel.any() and not expected_novel.all(), "test set must contain both verdicts"
