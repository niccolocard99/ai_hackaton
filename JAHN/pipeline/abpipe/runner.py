"""Runs the four stages, caching each module call and recording a manifest.

generate -> filter -> score -> rank. Each stage reads the previous stage's file from the
run directory, so a stage can be re-run alone (`abpipe stage score config.yaml`).
A module call is cached by (module, version, params, content of its input), so changing
one module never recomputes another, and expensive remote calls are paid for once.
"""
from __future__ import annotations

import dataclasses
import hashlib
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from importlib import metadata

import polars as pl

from abpipe import registry, schema
from abpipe.base import INFO, SCORE, Context
from abpipe.config import Config

STAGES = ("generate", "filter", "score", "rank", "validate")
_KEY_COLS = ["sequence_id", "kind", "heavy_sequence", "light_sequence"]


class PipelineEmpty(RuntimeError):
    """Every candidate was rejected; later stages have nothing to work on."""


def _log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def build_context(cfg: Config) -> Context:
    registry.load_plugins(cfg.plugins)
    run_dir = cfg.out_dir / cfg.name
    run_dir.mkdir(parents=True, exist_ok=True)
    cfg.cache_dir.mkdir(parents=True, exist_ok=True)
    ctx = Context(data_dir=cfg.data_dir, repo_root=cfg.repo_root, run_dir=run_dir,
                  cache_dir=cfg.cache_dir, seed=cfg.seed)
    if cfg.parent:
        from abpipe.data import load_parent
        ctx.parent = load_parent(cfg.parent, cfg.data_dir)
    return ctx


def _parent_key(ctx: Context):
    return dataclasses.asdict(ctx.parent) if ctx.parent else None


def _cached(ctx: Context, parts: dict, fn):
    key = hashlib.sha256(json.dumps(parts, sort_keys=True, default=str).encode()).hexdigest()[:24]
    path = ctx.cache_dir / f"{key}.parquet"
    if path.exists():
        return pl.read_parquet(path), True
    out = fn()
    out.write_parquet(path)
    return out, False


def _entry(spec, module, rows_out: int, t0: float, cached: bool, **extra) -> dict:
    return {"name": spec.name, "version": module.version, "params": spec.params,
            "rows_out": rows_out, "seconds": round(time.time() - t0, 3), "cached": cached,
            "cost_usd": 0.0 if cached else float(getattr(module, "cost_usd", 0.0)), **extra}


def _git(repo_root) -> dict:
    def g(*a):
        r = subprocess.run(["git", "-C", str(repo_root), *a], capture_output=True, text=True)
        return r.stdout.strip() if r.returncode == 0 else None
    return {"commit": g("rev-parse", "HEAD"), "dirty": bool(g("status", "--porcelain"))}


def _versions() -> dict:
    out = {"python": platform.python_version()}
    for pkg in ("polars", "abguard", "numpy", "rapidfuzz", "pyyaml"):
        try:
            out[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            out[pkg] = None
    return out


def _new_manifest(cfg: Config) -> dict:
    from abpipe import __version__
    return {"abpipe_version": __version__, "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "stages": {}}


def _order(df: pl.DataFrame) -> pl.DataFrame:
    head = [c for c in schema.REQUIRED if c in df.columns]
    return df.select(head + [c for c in df.columns if c not in head])


# ----------------------------------------------------------------------------- stages

def stage_generate(cfg: Config, ctx: Context, log=_log):
    frames, entries = [], []
    for i, spec in enumerate(cfg.generators):
        gen = registry.get("generator", spec.name)(**spec.params)
        seed = cfg.seed + i
        t0 = time.time()
        parts = {"stage": "generate", "name": spec.name, "version": gen.version, "params": spec.params,
                 "n": spec.n, "seed": seed, "parent": _parent_key(ctx)}
        df, cached = _cached(ctx, parts, lambda: gen.generate(spec.n, seed, ctx))
        schema.validate(df, f"generator {spec.name}")
        frames.append(schema.with_optional(df))
        entries.append(_entry(spec, gen, df.height, t0, cached, n_requested=spec.n, seed=seed))
        log(f"  generate  {spec.name:<24} {df.height:>7} rows{'  (cached)' if cached else ''}")
    out = _order(pl.concat(frames, how="diagonal_relaxed"))
    schema.validate(out, "combined generator output")
    return out, entries


def _check_filter_result(res: pl.DataFrame, df: pl.DataFrame, name: str) -> None:
    if "sequence_id" not in res.columns or "reason" not in res.columns:
        raise ValueError(f"filter {name}: must return columns sequence_id and reason")
    if res["sequence_id"].is_duplicated().any() or set(res["sequence_id"].to_list()) != set(df["sequence_id"].to_list()):
        raise ValueError(f"filter {name}: must return exactly one row for every input sequence_id")


def stage_filter(cfg: Config, ctx: Context, df: pl.DataFrame, log=_log, specs=None, label="filter"):
    rejections, entries, stopped = [], [], None
    for spec in (cfg.filters if specs is None else specs):
        flt = registry.get("filter", spec.name)(**spec.params)
        t0, rows_in = time.time(), df.height
        parts = {"stage": label, "name": spec.name, "version": flt.version, "params": spec.params,
                 "input": schema.content_hash(df, _KEY_COLS), "parent": _parent_key(ctx)}
        res, cached = _cached(ctx, parts, lambda: flt.check(df, ctx))
        _check_filter_result(res, df, spec.name)
        info = [c for c in res.columns if c.startswith(INFO)]
        df = df.drop([c for c in info if c in df.columns])
        if spec.on_fail == "drop":
            failed = res.filter(pl.col("reason").is_not_null()).select("sequence_id", "reason")
            if failed.height:
                rejections.append(failed.join(df.select("sequence_id", "heavy_sequence", "light_sequence"),
                                              on="sequence_id", how="left")
                                  .with_columns(module=pl.lit(spec.name))
                                  .select("sequence_id", "module", "reason", "heavy_sequence", "light_sequence"))
            keep = res.filter(pl.col("reason").is_null())["sequence_id"].to_list()
            df = df.filter(pl.col("sequence_id").is_in(keep)).join(
                res.select("sequence_id", *info), on="sequence_id", how="left")
            n_fail = failed.height
        else:
            flag = f"flag__{spec.name}"
            df = df.drop([flag] if flag in df.columns else []).join(
                res.select("sequence_id", pl.col("reason").alias(flag), *info), on="sequence_id", how="left")
            n_fail = res.filter(pl.col("reason").is_not_null()).height
        entries.append(_entry(spec, flt, df.height, t0, cached, rows_in=rows_in, failed=n_fail, on_fail=spec.on_fail))
        log(f"  {label:<9} {spec.name:<24} {rows_in:>7} -> {df.height:>7}  ({n_fail} {'dropped' if spec.on_fail == 'drop' else 'flagged'}){'  (cached)' if cached else ''}")
        if df.height == 0:
            stopped = spec.name
            break
    rej = (pl.concat(rejections) if rejections else
           pl.DataFrame(schema={"sequence_id": pl.String, "module": pl.String, "reason": pl.String,
                                "heavy_sequence": pl.String, "light_sequence": pl.String}))
    return df, rej, entries, stopped


def stage_score(cfg: Config, ctx: Context, df: pl.DataFrame, log=_log):
    directions, entries = {}, []
    for spec in cfg.scorers:
        scr = registry.get("scorer", spec.name)(**spec.params)
        outs = dict(scr.outputs)
        for col, d in outs.items():
            if not col.startswith(SCORE) or d not in ("higher", "lower"):
                raise ValueError(f"scorer {spec.name}: outputs must map score__* columns to 'higher'/'lower', got {col!r}: {d!r}")
            if col in directions:
                raise ValueError(f"scorer {spec.name}: output column {col} is already produced by another scorer")
        t0 = time.time()
        parts = {"stage": "score", "name": spec.name, "version": scr.version, "params": spec.params,
                 "input": schema.content_hash(df, _KEY_COLS), "parent": _parent_key(ctx)}
        res, cached = _cached(ctx, parts, lambda: scr.score(df, ctx))
        missing = [c for c in ["sequence_id", *outs] if c not in res.columns]
        if missing or set(res["sequence_id"].to_list()) != set(df["sequence_id"].to_list()):
            raise ValueError(f"scorer {spec.name}: must return sequence_id and {sorted(outs)} for every row (missing {missing})")
        df = df.drop([c for c in outs if c in df.columns]).join(res.select("sequence_id", *outs), on="sequence_id", how="left")
        directions.update(outs)
        entries.append(_entry(spec, scr, df.height, t0, cached, outputs=outs))
        log(f"  score     {spec.name:<24} {', '.join(outs)}{'  (cached)' if cached else ''}")
    return df, directions, entries


def stage_rank(cfg: Config, ctx: Context, df: pl.DataFrame, directions: dict, log=_log):
    spec = cfg.aggregator
    agg = registry.get("aggregator", spec.name)(**spec.params)
    t0 = time.time()
    ranked = agg.combine(df, directions, ctx)
    for col in ("rank_score", "rank"):
        if col not in ranked.columns:
            raise ValueError(f"aggregator {spec.name}: must add a {col} column")
    ranked = ranked.sort("rank")
    log(f"  rank      {spec.name:<24} {ranked.height} candidates ranked")
    return ranked, [_entry(spec, agg, ranked.height, t0, False)]


# ----------------------------------------------------------------------------- driver

def _rej_summary(rej: pl.DataFrame) -> dict:
    out = {}
    for module in rej["module"].unique().to_list():
        sub = rej.filter(pl.col("module") == module)
        top = sub.group_by("reason").len().sort("len", descending=True).head(3)
        out[module] = {"count": sub.height, "top_reasons": [{"reason": r, "count": n} for r, n in top.iter_rows()]}
    return out


def _write_table(df: pl.DataFrame, rd, name: str) -> None:
    """A CSV companion for a stage output, dropping nested columns CSV cannot hold."""
    scalar = [c for c, t in df.schema.items() if not isinstance(t, (pl.List, pl.Struct, pl.Array))]
    df.select(scalar).write_csv(rd / f"{name}.csv")


def _write_ranked(ranked: pl.DataFrame, rd) -> None:
    ranked.write_parquet(rd / "04_ranked.parquet")
    scalar = [c for c, t in ranked.schema.items() if not isinstance(t, (pl.List, pl.Struct, pl.Array))]
    ranked.select(scalar).write_csv(rd / "ranked.csv")


def run(cfg: Config, stages=STAGES, log=_log) -> dict:
    ctx = build_context(cfg)
    rd = ctx.run_dir
    mpath = rd / "manifest.json"
    man = json.loads(mpath.read_text()) if mpath.exists() else _new_manifest(cfg)
    cfg_sha = hashlib.sha256(json.dumps(cfg.raw, sort_keys=True, default=str).encode()).hexdigest()[:16]
    man.update({"config": cfg.raw, "config_sha256": cfg_sha, "git": _git(cfg.repo_root),
                "versions": _versions(), "parent": _parent_key(ctx),
                "updated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")})
    log(f"run {cfg.name}  (config {cfg_sha})")
    try:
        if "generate" in stages:
            df, entries = stage_generate(cfg, ctx, log)
            df.write_parquet(rd / "01_generated.parquet")
            man["stages"]["generate"] = {"modules": entries, "rows_out": df.height, "config_sha256": cfg_sha}
        if "filter" in stages:
            df = pl.read_parquet(rd / "01_generated.parquet")
            kept, rej, entries, stopped = stage_filter(cfg, ctx, df, log)
            kept.write_parquet(rd / "02_kept.parquet")
            rej.write_parquet(rd / "02_rejections.parquet")
            man["stages"]["filter"] = {"modules": entries, "rows_in": df.height, "rows_out": kept.height,
                                       "rejections": _rej_summary(rej), "config_sha256": cfg_sha}
            if stopped:
                raise PipelineEmpty(f"every candidate was rejected by filter {stopped!r}")
        if "score" in stages:
            df = pl.read_parquet(rd / "02_kept.parquet")
            scored, directions, entries = stage_score(cfg, ctx, df, log)
            scored.write_parquet(rd / "03_scored.parquet")
            (rd / "03_directions.json").write_text(json.dumps(directions, indent=2))
            man["stages"]["score"] = {"modules": entries, "rows_out": scored.height, "config_sha256": cfg_sha}
        if "rank" in stages:
            df = pl.read_parquet(rd / "03_scored.parquet")
            directions = json.loads((rd / "03_directions.json").read_text())
            ranked, entries = stage_rank(cfg, ctx, df, directions, log)
            _write_ranked(ranked, rd)
            man["stages"]["rank"] = {"modules": entries, "rows_out": ranked.height, "config_sha256": cfg_sha}
        if "validate" in stages and cfg.validators:
            # The expensive gate runs last, on the top of the ranked list only: abguard
            # costs ~0.7 ms/design to number and the novelty cross join 26 ms/design
            # unpatched, and only the top of the list can ever be submitted.
            ranked = pl.read_parquet(rd / "04_ranked.parquet")
            top = ranked.head(cfg.validate_top) if cfg.validate_top else ranked
            kept, rej, entries, stopped = stage_filter(cfg, ctx, top, log,
                                                       specs=cfg.validators, label="validate")
            kept = kept.sort("rank") if "rank" in kept.columns else kept
            kept.write_parquet(rd / "05_validated.parquet")
            rej.write_parquet(rd / "05_validate_rejections.parquet")
            _write_table(kept, rd, "validated")
            man["stages"]["validate"] = {"modules": entries, "rows_in": top.height, "rows_out": kept.height,
                                         "validate_top": cfg.validate_top,
                                         "rejections": _rej_summary(rej), "config_sha256": cfg_sha}
            if stopped:
                raise PipelineEmpty(f"every validated candidate was rejected by {stopped!r}")
    except PipelineEmpty as e:
        man["stopped"] = str(e)
        log(f"  STOP      {e}")
    man["totals"] = {"cost_usd": round(sum(m.get("cost_usd", 0.0) for s in man["stages"].values() for m in s["modules"]), 6),
                     "seconds": round(sum(m.get("seconds", 0.0) for s in man["stages"].values() for m in s["modules"]), 2)}
    mpath.write_text(json.dumps(man, indent=2, default=str))
    log(f"manifest: {mpath}")
    return man
