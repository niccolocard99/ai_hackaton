"""The dynamic loop driver: rounds of generate -> filter -> score -> select, steered between rounds.

    abpipe orchestrate configs/orchestrated_loop.yaml --rounds 6 --select 60 --out runs/orchestrated

Each round runs the ordinary abpipe stages (runner.run), then the survivors are ranked together with
the previous round's selected set (elitism) and the top `select` are kept with a diversity constraint
— the same contract as `abpipe iterate`. What iterate does not do is CHANGE the config between
rounds: before every round this driver

  1. probes the machine (GPU, heavy containers, disk, load; abpipe/probe.py),
  2. asks a policy for a patch to the next round's config (abpipe/policy.py: deterministic rules,
     or an external agent through the same allowlist),
  3. checks the time budget against the measured cost of a round, and stops when one more round
     cannot fit.

Every decision is journaled with its reasoning: orchestration.json (state), decisions.jsonl (one
line per round) and DECISIONS.md (human-readable, for the pitch's "clear reasoning" criterion).

What improves is the SURROGATE's score; read NICCOLO/eda/AGENT_GUIDE.md sections 3-6 before
trusting a ranking. Thermostability, epitope preservation and expression are not scored unless a
scorer for them is added to the config.
"""
from __future__ import annotations

import copy
import dataclasses
import hashlib
import json
import statistics
import time
from collections import Counter
from datetime import datetime, timedelta
from pathlib import Path

import polars as pl

from abpipe import policy as policymod, probe as probemod, registry, report, runner
from abpipe import config as cfgmod
from abpipe.iterate import EVOLVER, _round_stats, _strip_rank, _write_selected, select_diverse
from abpipe.generators.scored_evolution import file_sha256, history_digest

VHH_GEN = "vhh_variants"


# --------------------------------------------------------------------------- helpers

def parse_deadline(x) -> datetime | None:
    """ISO datetime, 'HH:MM' (today, or tomorrow if already past) or None."""
    if x is None or x == "":
        return None
    if isinstance(x, datetime):
        return x
    s = str(x).strip()
    try:
        if ":" in s and len(s) <= 5:
            t = datetime.strptime(s, "%H:%M")
            now = datetime.now()
            t = t.replace(year=now.year, month=now.month, day=now.day)
            if t <= now:
                t += timedelta(days=1)
            return t
        return datetime.fromisoformat(s.replace("Z", "+00:00"))
    except ValueError as e:
        raise ValueError(f"--deadline {s!r}: use ISO '2026-11-14T17:30' or '17:30' ({e})") from None


def _primary(cfg, override=None):
    """The score column the plateau rule watches: the aggregator's first column, or an explicit override."""
    directions: dict[str, str] = {}
    first = None
    for spec in cfg.scorers:
        outs = dict(registry.get("scorer", spec.name)(**spec.params).outputs)
        directions.update(outs)
        first = first or next(iter(outs), None)
    col = override or (list(cfg.aggregator.params.get("columns") or [])[0]
                       if cfg.aggregator.params.get("columns") else first)
    if col is None:
        return None, None
    if col not in directions:
        raise ValueError(f"primary score {col!r} is not produced by any scorer; available {sorted(directions)}")
    return col, directions[col]


def _effective(raw: dict) -> dict:
    """name -> current params (generators also carry n), for the policy state."""
    out: dict[str, dict] = {}
    for section in ("generators", "filters", "scorers"):
        for spec in raw.get(section) or []:
            spec = {"name": spec} if isinstance(spec, str) else spec
            merged = dict(spec.get("params") or {})
            if section == "generators" and spec.get("n") is not None:
                merged["n"] = spec["n"]
            out.setdefault(spec["name"], {}).update(merged)
    return out


def select_generic_diverse(ranked: pl.DataFrame, n_select: int, min_dist: int = 3) -> pl.DataFrame:
    """Best-first selection with a pairwise Levenshtein floor, for configs without a parent window
    (e.g. the de novo VHH lane, where candidates do not share the parent's length)."""
    if ranked.height == 0:
        return ranked
    from abpipe.distance import nearest
    buckets: dict = {}
    keep = []
    for i, row in enumerate(ranked.iter_rows(named=True)):
        seq = row.get("heavy_sequence") or ""
        best = nearest(seq, buckets, "levenshtein")
        if best is not None and best[0] < min_dist:
            continue
        seqs, ids = buckets.setdefault(len(seq), ([], []))
        seqs.append(seq)
        ids.append(row["sequence_id"])
        keep.append(i)
        if len(keep) == n_select:
            break
    return ranked[keep] if keep else ranked.head(0)


def _lite_stats(man: dict, sel: pl.DataFrame, scored: pl.DataFrame, pool_n: int, carried: int) -> dict:
    """Round stats for parentless configs: everything except the CDRH3-window diagnostics."""
    st = {"generated": man["stages"]["generate"]["rows_out"],
          "filters": [{"name": m["name"], "in": m["rows_in"], "out": m["rows_out"]}
                      for m in man["stages"]["filter"]["modules"]],
          "survivors": man["stages"]["filter"]["rows_out"], "pool_ranked": pool_n, "selected": sel.height,
          "carried_over": carried, "new_in_selected": sel.height - carried,
          "mean_pairwise_hamming": None,
          "origin_of_selected": dict(Counter(sel["origin"].to_list())) if "origin" in sel.columns else {}}
    for col in (c for c in sel.columns if c.startswith("score__")):
        v, sv = sel[col].drop_nulls(), scored[col].drop_nulls()
        st[col] = {"selected_median": round(float(v.median()), 4) if v.len() else None,
                   "selected_best": round(float(v.min() if col.endswith("kd") else v.max()), 4) if v.len() else None,
                   "survivors_median": round(float(sv.median()), 4) if sv.len() else None}
    return st


def _round_entry(r, snap, decision, budget, st=None) -> dict:
    entry = {"round": r, "utc": snap["utc"], "source": decision.source, "reasoning": decision.reasoning,
             "patch": decision.patch, "stop": decision.stop, "stop_reason": decision.stop_reason,
             "probe": {k: snap[k] for k in ("gpu", "heavy_images", "disk_free_gb", "load1",
                                            "japanfold_balance_usd") if k in snap},
             "budget": budget}
    if st:
        entry.update({k: st[k] for k in ("generated", "survivors", "selected", "carried_over",
                                         "new_in_selected", "wall_seconds", "selected_medians",
                                         "mean_pairwise_hamming") if k in st})
    return entry


def _md(path: Path, text: str) -> None:
    with open(path, "a") as fh:
        fh.write(text)


def _fmt_patch(patch: dict) -> str:
    parts = [f"loop {json.dumps(patch.get('loop') or {}, sort_keys=True)}"] if patch.get("loop") else []
    parts += [f"{name} {json.dumps(params, sort_keys=True)}"
              for name, params in sorted((patch.get("modules") or {}).items())]
    return "; ".join(parts) or "no change"


def _render_entry(e: dict) -> str:
    lines = [f"\n## Round {e['round']:02d} — {e['utc']}", f"- source: {e['source']}"]
    gpu = (e.get("probe") or {}).get("gpu")
    if gpu:
        lines.append(f"- probe: GPU {gpu['util_pct']}% ({gpu['mem_used_mb']}/{gpu['mem_total_mb']} MB), "
                     f"heavy containers {e['probe'].get('heavy_images') or 'none'}, "
                     f"disk {e['probe'].get('disk_free_gb')} GB free, load {e['probe'].get('load1')}")
    else:
        lines.append(f"- probe: {json.dumps(e.get('probe') or {}, default=str)}")
    b = e.get("budget") or {}
    if b.get("deadline_remaining_s") is not None:
        lines.append(f"- budget: {b['deadline_remaining_s'] / 60:.0f} min left, "
                     f"round costs about {(b.get('est_round_seconds') or 0) / 60:.0f} min")
    lines.append(f"- patch: {_fmt_patch(e['patch'])}")
    for why in e["reasoning"]:
        lines.append(f"- why: {why}")
    if e.get("generated") is not None:
        lines.append(f"- funnel: {e['generated']} generated -> {e['survivors']} survive -> "
                     f"{e['selected']} selected ({e['carried_over']} carried over) in {e.get('wall_seconds')}s")
        if e.get("selected_medians"):
            lines.append(f"- selected medians: {json.dumps(e['selected_medians'], sort_keys=True)}")
    if e["stop"]:
        lines.append(f"- **STOP: {e['stop_reason']}**")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- driver

def run_orchestrated(config, rounds=None, select=None, out_dir="runs/orchestrated", seed=None,
                     policy_kind=None, agent_cmd=None, deadline=None, patience=None, primary=None,
                     n=None, propose_next=False, dry_run=False, log=print) -> dict:
    base = cfgmod.load(config)
    registry.load_plugins(base.plugins)
    o = dict(base.raw.get("orchestrate") or {})
    base_dir = Path(config).parent if isinstance(config, (str, Path)) else None
    rounds = int(rounds if rounds is not None else o.get("rounds", 6))
    select = int(select if select is not None else o.get("select", 60))
    patience = int(patience if patience is not None else o.get("patience", 2))
    deadline = parse_deadline(deadline if deadline is not None else o.get("deadline"))
    pol = policymod.make_policy(policy_kind or o.get("policy", "rules"), agent_cmd or o.get("agent_cmd"))
    seed0 = seed if seed is not None else base.seed
    pcol, pdir = _primary(base, primary or o.get("primary"))
    lp = {"select": select, "min_pair_hamming": 3, "max_sub_share": 0.25, "pool_frac": 0.10}
    loop_dir = Path(out_dir)
    loop_dir.mkdir(parents=True, exist_ok=True)
    state_path, jsonl_path, md_path = loop_dir / "orchestration.json", loop_dir / "decisions.jsonl", loop_dir / "DECISIONS.md"

    state_doc = (json.loads(state_path.read_text()) if state_path.exists() else
                 {"config": str(config) if isinstance(config, (str, Path)) else "<inline>", "rounds": []})
    state_doc.update({"policy": {"kind": policy_kind or o.get("policy", "rules"), "agent_cmd": agent_cmd or o.get("agent_cmd")},
                       "deadline": deadline.isoformat() if deadline else None, "primary": pcol,
                       "primary_direction": pdir, "patience": patience, "select": select})
    prev, history = None, []
    for rd in sorted(p for p in loop_dir.glob("round_*") if (p / "selected.parquet").exists()):
        prev, history = pl.read_parquet(rd / "selected.parquet"), history + [rd / "01_generated.parquet"]
    walls = [r["wall_seconds"] for r in state_doc["rounds"] if r.get("wall_seconds")]
    est = statistics.mean(walls) if walls else None

    if not md_path.exists():
        _md(md_path, f"# Orchestrator decisions — {state_doc['config']}\n\n"
                     f"Policy: {state_doc['policy']['kind']}. Primary score: {pcol} "
                     f"({'lower is better' if pdir == 'lower' else 'higher is better'}).\n")

    raw_current = copy.deepcopy(base.raw)
    boot_hist = [Path(str(p)) for g in base.generators if g.name == EVOLVER
                 for p in (g.params.get("history") or [])]
    history = boot_hist + history          # config-declared evaluated windows are never re-proposed
    for r in range(rounds):
        rd = loop_dir / f"round_{r:02d}"
        if (rd / "selected.parquet").exists():
            log(f"round_{r:02d}: already finished, resuming after it")
            continue
        t0 = time.time()
        snap = probemod.probe(loop_dir)
        budget = {"deadline_remaining_s": probemod.seconds_until(deadline), "est_round_seconds": est}
        state = {"round": r, "rounds_total": rounds, "loop": lp, "effective": _effective(raw_current),
                 "rounds": state_doc["rounds"], "primary": pcol, "primary_direction": pdir,
                 "patience": patience, "probe": snap, **budget}
        decision = pol(state)
        entry = {"round": r, "utc": snap["utc"], "source": decision.source, "reasoning": decision.reasoning,
                 "patch": decision.patch, "probe": snap, "budget": budget}
        if decision.stop:
            entry.update({"stop": True, "stop_reason": decision.stop_reason})
            state_doc["stopped"] = f"round_{r:02d}: {decision.stop_reason}"
            state_doc["rounds"].append(_round_entry(r, snap, decision, budget))
            with open(jsonl_path, "a") as fh:
                fh.write(json.dumps(entry, default=str) + "\n")
            _md(md_path, _render_entry(_round_entry(r, snap, decision, budget)))
            state_path.write_text(json.dumps(state_doc, indent=2, default=str))
            log(f"round_{r:02d}: STOP — {decision.stop_reason}")
            break
        if dry_run:
            log(_render_entry({**entry, "stop": False, "stop_reason": ""}))
            log("dry run: nothing executed")
            return state_doc

        raw_current = policymod.apply_module_patch(raw_current, decision.patch)
        lp = policymod.apply_loop_patch(lp, decision.patch)
        cfg = cfgmod.load(raw_current, base_dir=base_dir)
        cfg = dataclasses.replace(cfg, name=f"round_{r:02d}", out_dir=loop_dir, cache_dir=loop_dir / "_cache",
                                  seed=seed0 + 1000 * r, generators=copy.deepcopy(cfg.generators))
        prev_dir = loop_dir / f"round_{r - 1:02d}" if prev is not None else None
        for spec in cfg.generators:
            if n is not None:
                spec.n = n
            if spec.name == EVOLVER:
                spec.params = {**spec.params, "n_elite": lp["select"], "history": [str(p) for p in history],
                               "history_sha256": history_digest(history) if history else None,
                               "parents": str(prev_dir / "selected.parquet") if prev_dir else spec.params.get("parents"),
                               "parents_sha256": (file_sha256(prev_dir / "selected.parquet") if prev_dir
                                                  else spec.params.get("parents_sha256"))}
            if spec.name == VHH_GEN and prev_dir is not None:
                spec.params = {**spec.params, "parents": str(prev_dir / "selected.csv"),
                               "parents_sha256": file_sha256(prev_dir / "selected.csv")}
        log(f"=== round_{r:02d} (seed {cfg.seed}) ===")
        rd.mkdir(parents=True, exist_ok=True)
        slog = report.StepLog(rd, t0, log)
        man = runner.run(cfg, log=slog)
        if man.get("stopped"):
            entry.update({"stop": True, "stop_reason": man["stopped"]})
            state_doc["stopped"] = f"round_{r:02d}: {man['stopped']}"
            state_doc["rounds"].append(entry)
            with open(jsonl_path, "a") as fh:
                fh.write(json.dumps(entry, default=str) + "\n")
            report.iteration_report(rd, entry, man)
            _md(md_path, _render_entry(entry))
            state_path.write_text(json.dumps(state_doc, indent=2, default=str))
            log(f"round_{r:02d}: STOP — {man['stopped']}")
            break
        scored = pl.read_parquet(rd / "03_scored.parquet")
        directions = json.loads((rd / "03_directions.json").read_text())
        ctx = runner.build_context(cfg)
        agg = registry.get("aggregator", cfg.aggregator.name)(**cfg.aggregator.params)
        pool = scored if prev is None else pl.concat([_strip_rank(prev), scored], how="diagonal_relaxed")
        ranked = agg.combine(pool, directions, ctx)
        if ctx.parent is not None:
            sel = select_diverse(ranked, ctx.parent, lp["select"], lp["min_pair_hamming"],
                                 lp["max_sub_share"], lp["pool_frac"])
        else:
            sel = select_generic_diverse(ranked, lp["select"], lp["min_pair_hamming"])
        carried = int(sel["sequence_id"].is_in(prev["sequence_id"].to_list()).sum()) if prev is not None else 0
        _write_selected(sel, rd)
        cfg_sha = hashlib.sha256(json.dumps(raw_current, sort_keys=True, default=str).encode()).hexdigest()[:16]
        (rd / "loop_round.json").write_text(json.dumps({"round": r, "config_sha256": cfg_sha, "select": lp["select"],
                                                         "min_pair_hamming": lp["min_pair_hamming"],
                                                         "max_sub_share": lp["max_sub_share"],
                                                         "pool_frac": lp["pool_frac"],
                                                         "patch_applied": decision.patch}))
        if ctx.parent is not None:
            st = _round_stats(man, sel, scored, ctx.parent, pool.height, carried)
        else:
            st = _lite_stats(man, sel, scored, pool.height, carried)
        st["selected_medians"] = {c: st[c]["selected_median"] for c in st if isinstance(st.get(c), dict)
                                  and "selected_median" in st[c]}
        st["round"], st["wall_seconds"] = r, round(time.time() - t0, 1)
        entry = _round_entry(r, snap, decision, budget, st)
        report.iteration_report(rd, entry, man)
        state_doc["rounds"] = [x for x in state_doc["rounds"] if x["round"] != r] + [entry]
        state_doc["rounds"].sort(key=lambda x: x["round"])
        walls.append(st["wall_seconds"])
        est = statistics.mean(walls[-3:])
        with open(jsonl_path, "a") as fh:
            fh.write(json.dumps(entry, default=str) + "\n")
        _md(md_path, _render_entry(entry))
        state_path.write_text(json.dumps(state_doc, indent=2, default=str))
        log(f"round_{r:02d}: {st['generated']} generated -> {st['survivors']} survive -> {sel.height} selected "
            f"({carried} carried over) in {st['wall_seconds']}s")
        prev, history = sel, history + [rd / "01_generated.parquet"]

    if propose_next and prev is not None and not state_doc.get("stopped"):
        r = len(history)
        rd = loop_dir / f"round_{r:02d}"
        cfg = cfgmod.load(raw_current, base_dir=base_dir)
        cfg = dataclasses.replace(cfg, name=f"round_{r:02d}", out_dir=loop_dir, cache_dir=loop_dir / "_cache",
                                  seed=seed0 + 1000 * r, generators=copy.deepcopy(cfg.generators))
        last = loop_dir / f"round_{r - 1:02d}"
        for spec in cfg.generators:
            if n is not None:
                spec.n = n
            if spec.name == EVOLVER:
                spec.params = {**spec.params, "n_elite": lp["select"],
                               "history": [str(p) for p in history],
                               "history_sha256": history_digest(history),
                               "parents": str(last / "selected.parquet"),
                               "parents_sha256": file_sha256(last / "selected.parquet")}
            if spec.name == VHH_GEN:
                sel_csv = last / "selected.csv"
                spec.params = {**spec.params, "parents": str(sel_csv),
                               "parents_sha256": file_sha256(sel_csv)}
        rem = probemod.seconds_until(deadline)
        if rem is not None and (est is None or rem > 0.5 * est):
            log(f"=== round_{r:02d}: proposals for the next iteration (generate + filter only) ===")
            runner.run(cfg, stages=("generate", "filter"), log=lambda m: log(m))
            kept = pl.read_parquet(rd / "02_kept.parquet")
            scalar = [c for c, t in kept.schema.items() if not isinstance(t, (pl.List, pl.Struct, pl.Array))]
            kept.select(scalar).write_csv(rd / "proposals.csv")
            state_doc["next_iteration_proposals"] = {"round": r, "n": kept.height, "file": f"round_{r:02d}/proposals.csv"}
            state_path.write_text(json.dumps(state_doc, indent=2, default=str))
        else:
            log("skipping --propose-next: too close to the deadline")
    return state_doc
