"""The fast novelty patch must not change a single abguard verdict.

`abpipe.abguard_fast` replaces `abguard.novelty._nearest` -- the 610k-wide cross join
that costs 26 ms/design -- with the deletion-neighbourhood index. These tests run the
real abguard both ways on the same designs and compare the columns the gate and the
report are built from. They need abguard and the organiser data, so they skip
elsewhere (see tests/conftest.py).
"""
from __future__ import annotations

import random

import polars as pl
import pytest

pytestmark = pytest.mark.needs_data

TZ_HEAVY = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAY"
    "LQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
TZ_LIGHT = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQP"
    "EDFATYYCQQHYTTPPTFGQGTKVEIK"
)
WINDOW = slice(99, 109)  # the 10 mutable CDRH3 residues inside TZ_HEAVY


def _designs(n: int, seed: int) -> pl.DataFrame:
    """Heavy chains carrying 0-8 substitutions in the CDRH3 window, plus the parent."""
    rng = random.Random(seed)
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    rows = [{"sequence_id": "trastuzumab", "kind": "paired", "heavy_sequence": TZ_HEAVY,
             "light_sequence": TZ_LIGHT, "design_description": "parent"}]
    for i in range(n):
        chars = list(TZ_HEAVY)
        for position in rng.sample(range(WINDOW.start, WINDOW.stop), rng.randint(0, 8)):
            chars[position] = rng.choice(alphabet)
        rows.append({
            "sequence_id": f"design_{i:04d}", "kind": "paired",
            "heavy_sequence": "".join(chars), "light_sequence": TZ_LIGHT,
            "design_description": "test",
        })
    return pl.DataFrame(rows)


@pytest.fixture(scope="module")
def reference():
    from abguard import reference as refmod

    bundle = refmod.bundled_reference()
    if bundle is None:
        pytest.skip("no bundled abguard reference")
    return refmod.cdrh3_reference(bundle)


def _annotate(frame: pl.DataFrame, reference) -> pl.DataFrame:
    from abguard.annotate import annotate
    from abguard.config import PipelineConfig

    return annotate(frame, PipelineConfig(), include_liabilities=False,
                    reference=reference, gate_novelty=True)


COMPARED = [
    "sequence_id", "cdrh3", "cdrh3_min_distance", "is_novel_enough",
    "novelty_class", "required_cdrh3_mutations", "errors", "is_valid",
]


def test_patch_preserves_every_novelty_verdict(reference):
    from abpipe import abguard_fast

    frame = _designs(400, seed=14)

    abguard_fast.uninstall()
    baseline = _annotate(frame, reference).select(COMPARED)

    abguard_fast.install()
    try:
        patched = _annotate(frame, reference).select(COMPARED)
    finally:
        abguard_fast.uninstall()

    # The gate decision, the class and the hard errors must be byte-identical.
    for column in ("sequence_id", "cdrh3", "is_novel_enough", "novelty_class",
                   "required_cdrh3_mutations", "errors", "is_valid"):
        assert patched.get_column(column).to_list() == baseline.get_column(column).to_list(), column

    # Distances agree exactly wherever the index resolves them, which is everywhere
    # the gate could reject: a failing design is within the cap by construction.
    cap = abguard_fast.max_distance()
    both = baseline.select("cdrh3_min_distance").rename({"cdrh3_min_distance": "base"}).hstack(
        patched.select("cdrh3_min_distance").rename({"cdrh3_min_distance": "fast"})
    )
    resolved = both.filter(pl.col("base") <= cap)
    assert resolved.height > 0, "test set must contain designs near the reference"
    assert resolved.get_column("base").to_list() == resolved.get_column("fast").to_list()
    saturated = both.filter(pl.col("base") > cap)
    assert (saturated.get_column("fast") == cap + 1).all(), "far designs must report the sentinel"


def test_patch_reports_the_same_rejection_reason(reference):
    """Designs that fail novelty must name the same colliding reference loop."""
    from abpipe import abguard_fast

    frame = _designs(400, seed=3)

    abguard_fast.uninstall()
    baseline = _annotate(frame, reference)
    abguard_fast.install()
    try:
        patched = _annotate(frame, reference)
    finally:
        abguard_fast.uninstall()

    failed = baseline.get_column("is_novel_enough").fill_null(False).not_()
    assert failed.sum() > 0, "test set must contain novelty failures"
    for column in ("cdrh3_nearest", "cdrh3_nearest_source", "novelty_reason"):
        if column not in baseline.columns:
            continue
        assert (baseline.filter(failed).get_column(column).to_list()
                == patched.filter(failed).get_column(column).to_list()), column


def test_install_is_idempotent_and_reversible(reference):
    from abguard import novelty

    from abpipe import abguard_fast

    original = novelty._nearest
    abguard_fast.install()
    abguard_fast.install()
    assert novelty._nearest is not original
    abguard_fast.uninstall()
    assert novelty._nearest is original
