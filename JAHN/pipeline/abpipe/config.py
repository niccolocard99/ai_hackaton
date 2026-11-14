"""YAML config -> validated Config. Swapping a module means editing this file, not code."""
from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import yaml

TOP_KEYS = {"run", "data_dir", "repo_root", "plugins", "parent", "generators", "generator",
            "filters", "scorers", "aggregator", "validators", "validate_top"}
RUN_KEYS = {"name", "out_dir", "cache_dir", "seed"}


class ConfigError(ValueError):
    pass


@dataclass
class ModuleSpec:
    name: str
    params: dict = field(default_factory=dict)
    on_fail: str = "drop"      # filters only: drop | flag
    n: int | None = None       # generators only


@dataclass
class Config:
    name: str
    out_dir: Path
    cache_dir: Path
    seed: int
    data_dir: Path
    repo_root: Path
    plugins: list
    parent: dict | None
    generators: list
    filters: list
    scorers: list
    aggregator: ModuleSpec
    validators: list          # filters that run AFTER ranking, on the top `validate_top`
    validate_top: int | None  # None: validate the whole ranked list
    raw: dict


def _expand(obj):
    if isinstance(obj, str):
        return os.path.expanduser(os.path.expandvars(obj))
    if isinstance(obj, list):
        return [_expand(x) for x in obj]
    if isinstance(obj, dict):
        return {k: _expand(v) for k, v in obj.items()}
    return obj


def _spec(item, where: str, *, kind: str) -> ModuleSpec:
    if isinstance(item, str):
        item = {"name": item}
    if not isinstance(item, dict) or "name" not in item:
        raise ConfigError(f"{where}: each module needs a `name`")
    allowed = {"name", "params"} | ({"on_fail"} if kind == "filter" else set()) | ({"n"} if kind == "generator" else set())
    extra = set(item) - allowed
    if extra:
        raise ConfigError(f"{where} ({item['name']}): unexpected keys {sorted(extra)}; allowed {sorted(allowed)}")
    spec = ModuleSpec(name=item["name"], params=dict(item.get("params") or {}),
                      on_fail=item.get("on_fail", "drop"), n=item.get("n"))
    if spec.on_fail not in ("drop", "flag"):
        raise ConfigError(f"{where} ({spec.name}): on_fail must be drop or flag")
    if kind == "generator" and (not isinstance(spec.n, int) or spec.n <= 0):
        raise ConfigError(f"{where} ({spec.name}): generators need a positive integer `n`")
    return spec


def load(source, base_dir: Path | None = None) -> Config:
    """`source` is a path to a YAML file or an already-parsed dict."""
    if isinstance(source, (str, Path)):
        path = Path(source)
        raw = yaml.safe_load(path.read_text())
        base_dir = base_dir or path.parent
    else:
        raw, base_dir = dict(source), (base_dir or Path.cwd())
    raw = _expand(raw or {})
    extra = set(raw) - TOP_KEYS
    if extra:
        raise ConfigError(f"unknown top-level keys {sorted(extra)}; allowed {sorted(TOP_KEYS)}")
    run = raw.get("run") or {}
    if set(run) - RUN_KEYS:
        raise ConfigError(f"unknown run keys {sorted(set(run) - RUN_KEYS)}; allowed {sorted(RUN_KEYS)}")
    if "data_dir" not in raw:
        raise ConfigError("`data_dir` is required (the organiser data folder)")
    gens = raw.get("generators") or ([raw["generator"]] if "generator" in raw else [])
    if not gens:
        raise ConfigError("at least one generator is required")
    if "aggregator" not in raw:
        raise ConfigError("`aggregator` is required")
    validate_top = raw.get("validate_top")
    if validate_top is not None and (not isinstance(validate_top, int) or validate_top <= 0):
        raise ConfigError("`validate_top` must be a positive integer (omit it to validate the whole list)")
    out_dir = Path(run.get("out_dir", "runs"))
    repo_root = raw.get("repo_root")
    if repo_root is None:
        top = subprocess.run(["git", "-C", str(base_dir), "rev-parse", "--show-toplevel"],
                             capture_output=True, text=True)
        repo_root = top.stdout.strip() if top.returncode == 0 else str(base_dir)
    return Config(
        name=run.get("name", "run"), out_dir=out_dir,
        cache_dir=Path(run["cache_dir"]) if run.get("cache_dir") else out_dir / "_cache",
        seed=int(run.get("seed", 0)), data_dir=Path(raw["data_dir"]), repo_root=Path(repo_root),
        plugins=list(raw.get("plugins") or []), parent=raw.get("parent"),
        generators=[_spec(g, "generators", kind="generator") for g in gens],
        filters=[_spec(f, "filters", kind="filter") for f in raw.get("filters") or []],
        scorers=[_spec(s, "scorers", kind="scorer") for s in raw.get("scorers") or []],
        aggregator=_spec(raw["aggregator"], "aggregator", kind="aggregator"),
        validators=[_spec(v, "validators", kind="filter") for v in raw.get("validators") or []],
        validate_top=validate_top, raw=raw)
