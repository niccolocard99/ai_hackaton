#!/usr/bin/env python3
"""Check designs against the OPTIMISATION track rules (ai-antibodies-hackathon-apheris.md, Figure 2).

Shared hard checks (abguard) plus:
  * identity to nearest known antibody >= 70% (abguard's `novel` class; else it is de novo)
  * CDRH3 identity to trastuzumab strictly < 70%
  * strictly MORE than 2 point mutations (>= 3) from any provided data point:
    CDRH3 vs the abguard reference (all 525,063 HER2-table CDRH3s + prior art), and,
    if --provided is given, whole VH+VL vs every antibody in those tables.

    ./check_optimisation.py designs.csv [--provided trastuzumab_her2_spr_all_768.csv] \\
        [--reference-db bundle/] [-o report.csv] [--strict-exit]
"""
from pathlib import Path

from _track_common import *  # noqa: F401,F403
from _track_common import (IDENTITY_THRESHOLD, annotated, build_parser, cdrh3_identity, finish,
                           levenshtein, pl, read_table, trastuzumab_cdrh3)

MIN_POINT_MUTATIONS = 3  # "more than 2"; see AGENTS.md off-by-one note


def provided_pairs(path: Path) -> list[tuple[str, str]]:
    t = read_table(path)
    if "vh_sequence" not in t.columns:
        raise SystemExit(f"{path}: needs a vh_sequence column (and optionally vl_sequence)")
    vl = t["vl_sequence"] if "vl_sequence" in t.columns else pl.Series([""] * t.height)
    return list({(h or "", l or "") for h, l in zip(t["vh_sequence"], vl.fill_null(""))})


def main() -> int:
    p = build_parser(__doc__)
    p.add_argument("--provided", type=Path, help="table with vh_sequence/vl_sequence of provided antibodies (e.g. the SPR csv)")
    args = p.parse_args()
    df = annotated(args.submission, args.reference_db)
    tz = trastuzumab_cdrh3()
    provided = provided_pairs(args.provided) if args.provided else []
    ident, failures = [], []
    for r in df.iter_rows(named=True):
        f = [f"hard check: {e}" for e in r["errors"]]
        c = r["cdrh3"]
        ident.append(None if c is None else cdrh3_identity(c, tz))
        if not f:
            if r["novelty_class"] != "novel":
                f.append(f"not an optimisation design: {r['identity_to_trastuzumab']:.0%} identity to a known antibody "
                         "(< 70% is de novo); rejected here, not reclassified")
            else:
                if ident[-1] >= IDENTITY_THRESHOLD:
                    f.append(f"CDRH3 identity to trastuzumab {ident[-1]:.0%}; needs < 70%")
                if r["cdrh3_min_distance"] < MIN_POINT_MUTATIONS:
                    f.append(f"CDRH3 only {r['cdrh3_min_distance']} mutation(s) from {r['cdrh3_nearest_source']} "
                             f"({r['cdrh3_nearest']}); needs >= {MIN_POINT_MUTATIONS}")
                if provided:
                    h, l = r["heavy_sequence"], r["light_sequence"] or ""
                    best = min(levenshtein(h, ph) + levenshtein(l, pl_) for ph, pl_ in provided) if r["kind"] == "paired" \
                        else min(levenshtein(h, ph) for ph, _ in provided)
                    if best < MIN_POINT_MUTATIONS:
                        f.append(f"only {best} mutation(s) from a provided antibody; needs >= {MIN_POINT_MUTATIONS}")
        failures.append(f)
    df = df.with_columns(
        track_assigned=pl.col("novelty_class"),
        cdrh3_identity_to_trastuzumab=pl.Series(ident, dtype=pl.Float64),
        failures=pl.Series(failures, dtype=pl.List(pl.String)),
    ).with_columns(passes=pl.col("failures").list.len() == 0)
    return finish(df, "optimisation", args)


if __name__ == "__main__":
    raise SystemExit(main())
