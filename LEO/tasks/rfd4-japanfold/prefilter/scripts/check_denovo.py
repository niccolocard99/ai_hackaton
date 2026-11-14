#!/usr/bin/env python3
"""Check designs against the DE NOVO track rules (ai-antibodies-hackathon-apheris.md, Figure 2).

Shared hard checks (abguard): numbering confidence >= 0.1, format VHH/scFv/paired,
length <= 250, no unpaired cysteines, unique ids.
De novo: identity to nearest known antibody < 70% (abguard measures it against
trastuzumab's variable domains, i.e. its `de_novo` class), and CDRH3 at least
2 mutations from trastuzumab and from every CDRH3 in the reference.

    ./check_denovo.py designs.csv [--reference-db bundle/] [-o report.csv] [--strict-exit]
"""
from _track_common import *  # noqa: F401,F403
from _track_common import annotated, build_parser, cdrh3_identity, finish, pl, trastuzumab_cdrh3

MIN_CDRH3_MUTATIONS = 2


def main() -> int:
    args = build_parser(__doc__).parse_args()
    df = annotated(args.submission, args.reference_db)
    tz = trastuzumab_cdrh3()
    failures = []
    for r in df.iter_rows(named=True):
        f = [f"hard check: {e}" for e in r["errors"]]
        if not f:
            if r["novelty_class"] != "de_novo":
                f.append(f"not de novo: {r['identity_to_trastuzumab']:.0%} identity to a known antibody (needs < 70%); "
                         "it would be an optimisation design and is rejected here, not reclassified")
            elif r["cdrh3_min_distance"] < MIN_CDRH3_MUTATIONS:
                f.append(f"CDRH3 only {r['cdrh3_min_distance']} mutation(s) from {r['cdrh3_nearest_source']} "
                         f"({r['cdrh3_nearest']}); needs >= {MIN_CDRH3_MUTATIONS}")
        failures.append(f)
    df = df.with_columns(
        track_assigned=pl.col("novelty_class"),
        cdrh3_identity_to_trastuzumab=pl.col("cdrh3").map_elements(
            lambda c: None if c is None else cdrh3_identity(c, tz), return_dtype=pl.Float64),
        failures=pl.Series(failures, dtype=pl.List(pl.String)),
    ).with_columns(passes=pl.col("failures").list.len() == 0)
    return finish(df, "de novo", args)


if __name__ == "__main__":
    raise SystemExit(main())
