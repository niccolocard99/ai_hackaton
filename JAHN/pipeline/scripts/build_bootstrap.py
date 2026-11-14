"""Build the warm-start tables for the orchestrated optimisation loop from existing team work.

Sources (all committed on main):
  - LEO/tasks/design-loop/out/FINAL.csv       the team's ranked final shortlist (Tm, epitope, fold-aware)
  - LEO/tasks/design-loop/out/pool.csv        every window LEO's screen already evaluated
  - LEO/tasks/rfd4-japanfold/prefilter/leo_candidates_scored.csv    geometry-driven RFD4 windows (new)

Outputs (under JAHN/pipeline/bootstrap/):
  parents.csv   ranked elites: FINAL by its rank, then the RFD4 candidates by predicted KD. The
                generator's CEM pulls its residue/position/substitution distributions towards these.
  history.csv   every window already evaluated (pool + RFD4 gated + FINAL): the generator never
                proposes any of them again.

The window (10-mer) is grafted back onto trastuzumab's heavy chain at the parent offset, which is
exact for every source: all of them are framework-fixed CDRH3 variants (LEO's RFD4 verified the
framework held residue by residue; rfd4-japanfold's RESULTS.md). Rerun this script after any source
changes; the config pins both files by sha256, so a stale bootstrap cannot be used silently.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

REPO = Path(__file__).resolve().parents[3]
LEO = REPO / "LEO/tasks/design-loop/out"
RFD4 = REPO / "LEO/tasks/rfd4-japanfold/prefilter"
OUT = Path(__file__).resolve().parents[1] / "bootstrap"

COLUMNS = ["sequence_id", "kind", "heavy_sequence", "light_sequence", "design_description"]


def graft(cdrh3: str, heavy: str, start: int, L: int = 10) -> str:
    return heavy[:start] + cdrh3 + heavy[start + L:]


def main() -> None:
    from abpipe import _niccolo
    from abpipe.data import load_parent

    parent = load_parent({"table": _niccolo.SPR_REL, "window_column": "cdrh3", "id": "trastuzumab"},
                         REPO / "apheris/hackathon-data")
    s = parent.window_start

    final = pl.read_csv(LEO / "FINAL.csv")
    rows, rank = [], 0
    for r in final.sort("rank").iter_rows(named=True):
        rank += 1
        rows.append({"sequence_id": r["id"], "kind": "paired", "heavy_sequence": r["heavy_sequence"],
                     "light_sequence": r["light_sequence"],
                     "design_description": f"bootstrap: LEO final slot {r['slot']} ({r['addons']}); "
                                          f"kd_mean {r['kd_mean']:.2f}, dtm {r['dtm_vs_wt']}",
                     "source": "leo-final", "rank": rank})
    rfd4 = pl.read_csv(RFD4 / "leo_candidates_scored.csv").sort("pred_kd_nM")
    for i, r in enumerate(rfd4.iter_rows(named=True), 1):
        rank += 1
        rows.append({"sequence_id": f"rfd4_{r['cdrh3']}", "kind": "paired",
                     "heavy_sequence": graft(r["cdrh3"], parent.heavy, s), "light_sequence": parent.light,
                     "design_description": f"bootstrap: RFD4-on-JapanFold geometry design ({r['who']}); "
                                          f"pred KD {r['pred_kd_nM']:.2f} nM, p_binder {r['p_binder']:.2f}",
                     "source": "rfd4-japanfold", "rank": rank})
    parents = pl.DataFrame(rows).select(*COLUMNS, "source", "rank").sort("rank")
    bad = parents.filter(pl.col("heavy_sequence").str.len_chars() != len(parent.heavy))
    assert bad.height == 0, f"graft produced wrong-length chains: {bad['sequence_id'].to_list()[:3]}"

    seen, hist = set(), []
    for path, label in ((LEO / "pool.csv", "leo-pool"), (RFD4 / "rfd4_designs_gated.csv", "rfd4-gated"),
                        (LEO / "FINAL.csv", "leo-final"), (LEO / "tm_arm.csv", "leo-tm-arm")):
        t = pl.read_csv(path)
        for r in t.iter_rows(named=True):
            w = r.get("cdrh3") or r.get("base_cdrh3")
            if not w or len(w) != 10 or w in seen or not set(w) <= set(_niccolo.AA):
                continue
            seen.add(w)
            hist.append({"sequence_id": f"hist_{label}_{w}", "kind": "paired",
                         "heavy_sequence": graft(w, parent.heavy, s), "light_sequence": parent.light,
                         "design_description": f"already evaluated by {label}"})
    history = pl.DataFrame(hist, schema={c: pl.String for c in COLUMNS})
    assert history.height == len(seen) and history["heavy_sequence"].n_unique() == len(seen)

    OUT.mkdir(exist_ok=True)
    parents.write_csv(OUT / "parents.csv")
    history.write_csv(OUT / "history.csv")
    import hashlib
    for f in ("parents.csv", "history.csv"):
        print(f, (OUT / f).stat().st_size, "bytes, sha256",
              hashlib.sha256((OUT / f).read_bytes()).hexdigest()[:16])
    print(f"parents: {parents.height} elites ({parents['source'].value_counts().sort('source').to_dicts()})")
    print(f"history: {history.height} distinct evaluated windows")


if __name__ == "__main__":
    main()
