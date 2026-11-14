"""Round driver: generate -> filter -> score -> rank, then feed the best into the next round's generator.

    abpipe iterate CONFIG --rounds 3 --select 60 --out runs/loop [--propose-next]

Round 0 starts from the priors. In round r >= 1 every `scored_evolution` generator in the config is given the previous
round's selected windows (`parents`) and every window already evaluated (`history`). The new survivors are ranked
TOGETHER with the carried-over elites on their absolute scores (elitism). Only the single best candidate is
guaranteed not to get worse: the selected set as a whole can shift, because the diversity constraints re-pick around
better neighbours (a child within 2 substitutions of an elite displaces it if it scores better). Selection looks only
at the best `--pool-frac` of the ranking and returns fewer than `--select` rather than pad the list with mid-pool
candidates. Other generators in the config run unchanged each round. Re-running with a larger --rounds resumes after
the last finished round.

What improves is the SURROGATE's score. The loop cannot see binding, epitope or thermostability; only the scorers
in the config can, and the pre-filter in particular was validated on library variants, not on novel designs.
"""
from __future__ import annotations

import copy
import dataclasses
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
import polars as pl

from abpipe import _niccolo, novelty, registry, runner
from abpipe import config as cfgmod
from abpipe.generators.scored_evolution import file_sha256, history_digest

_RANK = ("rank", "rank_score")
EVOLVER = "scored_evolution"


def _strip_rank(df: pl.DataFrame) -> pl.DataFrame:
    return df.drop([c for c in df.columns if c in _RANK or c.startswith("rank__")])


def select_diverse(ranked: pl.DataFrame, parent, n_select: int, min_pair_hamming: int = 3,
                   max_sub_share: float = 0.25, pool_frac: float = 0.10) -> pl.DataFrame:
    """Best-first walk over the top `pool_frac` of a ranked table, keeping a window only if it differs from every one
    already kept in >= min_pair_hamming positions and no single substitution (e.g. A106E) appears in more than
    max_sub_share of the RETURNED list (never fewer than 2 uses).

    Two failure modes this is built to avoid, both seen in real runs. A cap that grows with the list as it is built
    allows only two uses of any substitution among the first eight picks, so when strong candidates share a few favoured
    residues the walk is starved until it reaches the middle of the pool. A cap tied to the requested size instead lets
    15 of a 17-long list share one substitution. So the cap is tied to the size that is actually returned: if the walk
    comes back shorter than assumed, the assumed size shrinks and the walk is repeated until the list and its cap agree.
    The pool limit means fewer than n_select rows come back, rather than padding the list with mid-pool candidates."""
    wins = _niccolo.windows(ranked["heavy_sequence"].to_list(), parent)
    limit = max(n_select, math.ceil(pool_frac * ranked.height))
    rows = [i for i, w in enumerate(wins) if w is not None and i < limit]
    if not rows:
        return ranked.head(0)
    S = novelty.encode([wins[i] for i in rows])
    wt = novelty.encode([parent.window])[0]
    subs_of = [[(int(p), int(S[r, p])) for p in np.flatnonzero(S[r] != wt)] for r in range(len(S))]

    def walk(size: int) -> list:
        cap = max(2, math.ceil(max_sub_share * size))
        chosen, counts = [], Counter()
        for r in range(len(S)):
            if chosen and (S[chosen] != S[r]).sum(1).min() < min_pair_hamming:
                continue
            if any(counts[s] + 1 > cap for s in subs_of[r]):
                continue
            chosen.append(r)
            counts.update(subs_of[r])
            if len(chosen) == size:
                break
        return chosen

    size = n_select
    while True:
        chosen = walk(size)
        if len(chosen) >= size or size <= 1:
            break
        size = len(chosen)
    return ranked[[rows[i] for i in chosen]]


def _write_selected(sel: pl.DataFrame, rd: Path) -> None:
    sel.write_parquet(rd / "selected.parquet")
    scalar = [c for c, t in sel.schema.items() if not isinstance(t, (pl.List, pl.Struct, pl.Array))]
    sel.select(scalar).write_csv(rd / "selected.csv")


def _round_stats(man: dict, sel: pl.DataFrame, scored: pl.DataFrame, parent, pool_n: int, carried: int) -> dict:
    gen = man["stages"]["generate"]
    flt = man["stages"]["filter"]["modules"]
    wins = novelty.encode([w for w in _niccolo.windows(sel["heavy_sequence"].to_list(), parent) if w])
    wt = novelty.encode([parent.window])[0]
    nm = (wins != wt).sum(1)
    pair = [int((wins[i] != wins[j]).sum()) for i in range(len(wins)) for j in range(i + 1, len(wins))] if len(wins) > 1 else []
    st = {"generated": gen["rows_out"],
          "filters": [{"name": m["name"], "in": m["rows_in"], "out": m["rows_out"]} for m in flt],
          "survivors": man["stages"]["filter"]["rows_out"], "pool_ranked": pool_n, "selected": sel.height,
          "carried_over": carried, "new_in_selected": sel.height - carried,
          "n_mut_selected": {int(k): int(v) for k, v in zip(*np.unique(nm, return_counts=True))},
          "mean_pairwise_hamming": round(float(np.mean(pair)), 2) if pair else None,
          "origin_of_selected": dict(Counter(sel["origin"].to_list())) if "origin" in sel.columns else {}}
    for col in (c for c in sel.columns if c.startswith("score__")):
        v = sel[col].drop_nulls()
        sv = scored[col].drop_nulls()
        st[col] = {"selected_median": round(float(v.median()), 4) if v.len() else None,
                   "selected_best": round(float(v.min() if col.endswith("kd") else v.max()), 4) if v.len() else None,
                   "survivors_median": round(float(sv.median()), 4) if sv.len() else None}
    return st


def run_loop(config_path, rounds: int = 3, n_select: int = 60, out_dir="runs/loop", seed: int | None = None,
             n: int | None = None, min_pair_hamming: int = 3, max_sub_share: float = 0.25,
             pool_frac: float = 0.10, propose_next: bool = False, log=print) -> dict:
    base = cfgmod.load(config_path)
    loop_dir = Path(out_dir)
    loop_dir.mkdir(parents=True, exist_ok=True)
    if not any(g.name == EVOLVER for g in base.generators):
        log(f"note: no `{EVOLVER}` generator in the config, so rounds will not learn from each other")
    summary_path = loop_dir / "loop_summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {"rounds": []}
    summary.update({"config": str(config_path) if isinstance(config_path, (str, Path)) else "<inline dict>",
                    "select": n_select})
    cfg_sha = hashlib.sha256(json.dumps(base.raw, sort_keys=True, default=str).encode()).hexdigest()[:16]
    prev, history = None, []
    for r in range(rounds):
        name = f"round_{r:02d}"
        rd = loop_dir / name
        done = rd / "selected.parquet"
        if done.exists():
            meta = json.loads((rd / "loop_round.json").read_text())
            same = (meta["config_sha256"] == cfg_sha and meta["select"] == n_select and meta.get("pool_frac") == pool_frac
                    and meta.get("max_sub_share") == max_sub_share and meta.get("min_pair_hamming") == min_pair_hamming)
            if not same:
                raise SystemExit(f"{name} was produced by a different config or selection settings; use a new --out directory")
            log(f"{name}: already finished, resuming after it")
            prev, history = pl.read_parquet(done), history + [rd / "01_generated.parquet"]
            continue
        cfg = dataclasses.replace(base, name=name, out_dir=loop_dir, cache_dir=loop_dir / "_cache",
                                  seed=(seed if seed is not None else base.seed) + 1000 * r,
                                  generators=copy.deepcopy(base.generators))
        for spec in cfg.generators:
            if n is not None:
                spec.n = n
            if spec.name == EVOLVER:
                spec.params = {**spec.params, "n_elite": n_select, "history": [str(p) for p in history],
                               "history_sha256": history_digest(history) if history else None,
                               "parents": str(loop_dir / f"round_{r - 1:02d}" / "selected.parquet") if prev is not None else None,
                               "parents_sha256": file_sha256(loop_dir / f"round_{r - 1:02d}" / "selected.parquet") if prev is not None else None}
        log(f"=== {name}  (seed {cfg.seed}) ===")
        man = runner.run(cfg, log=lambda m: log(m))
        if man.get("stopped"):
            summary["stopped"] = f"{name}: {man['stopped']}"
            break
        scored = pl.read_parquet(rd / "03_scored.parquet")
        if cfg.validators:
            # `03_scored` is upstream of the validate stage, which runs after ranking so
            # the expensive gate only sees the top of the list. Selection must still never
            # promote a design that gate rejected -- it would seed the next round and reach
            # selected.csv with no symptom until submission -- so restrict to its survivors.
            validated = pl.read_parquet(rd / "05_validated.parquet")["sequence_id"].to_list()
            scored = scored.filter(pl.col("sequence_id").is_in(validated))
        directions = json.loads((rd / "03_directions.json").read_text())
        ctx = runner.build_context(cfg)
        agg = registry.get("aggregator", cfg.aggregator.name)(**cfg.aggregator.params)
        pool = scored if prev is None else pl.concat([_strip_rank(prev), scored], how="diagonal_relaxed")
        ranked = agg.combine(pool, directions, ctx)
        sel = select_diverse(ranked, ctx.parent, n_select, min_pair_hamming, max_sub_share, pool_frac)
        carried = int(sel["sequence_id"].is_in(prev["sequence_id"].to_list()).sum()) if prev is not None else 0
        _write_selected(sel, rd)
        (rd / "loop_round.json").write_text(json.dumps({"config_sha256": cfg_sha, "select": n_select, "pool_frac": pool_frac,
                                                        "max_sub_share": max_sub_share, "min_pair_hamming": min_pair_hamming}))
        st = {"round": r, **_round_stats(man, sel, scored, ctx.parent, pool.height, carried)}
        summary["rounds"] = [x for x in summary["rounds"] if x["round"] != r] + [st]
        summary["rounds"].sort(key=lambda x: x["round"])
        summary_path.write_text(json.dumps(summary, indent=2, default=str))
        log(f"{name}: {st['generated']} generated -> {st['survivors']} survive filters -> {sel.height} of {n_select} selected "
            f"({carried} carried over)" + ("" if sel.height == n_select else
                                           f"; the top {pool_frac:.0%} held only {sel.height} mutually diverse candidates"))
        prev, history = sel, history + [rd / "01_generated.parquet"]
    if propose_next and prev is not None and not summary.get("stopped"):
        r = len(history)
        name = f"round_{r:02d}"
        cfg = dataclasses.replace(base, name=name, out_dir=loop_dir, cache_dir=loop_dir / "_cache",
                                  seed=(seed if seed is not None else base.seed) + 1000 * r,
                                  generators=copy.deepcopy(base.generators))
        last = loop_dir / f"round_{r - 1:02d}" / "selected.parquet"
        for spec in cfg.generators:
            if n is not None:
                spec.n = n
            if spec.name == EVOLVER:
                spec.params = {**spec.params, "n_elite": n_select, "history": [str(p) for p in history],
                               "history_sha256": history_digest(history), "parents": str(last),
                               "parents_sha256": file_sha256(last)}
        log(f"=== {name}: proposals for the next iteration (generate + filter only, not scored) ===")
        runner.run(cfg, stages=("generate", "filter"), log=lambda m: log(m))
        kept = pl.read_parquet(loop_dir / name / "02_kept.parquet")
        scalar = [c for c, t in kept.schema.items() if not isinstance(t, (pl.List, pl.Struct, pl.Array))]
        kept.select(scalar).write_csv(loop_dir / name / "proposals.csv")
        summary["next_iteration_proposals"] = {"round": r, "n": kept.height, "file": f"{name}/proposals.csv"}
        summary_path.write_text(json.dumps(summary, indent=2, default=str))
    return summary
