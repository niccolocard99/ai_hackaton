#!/usr/bin/env python3
"""Generation ledger: stop the team paying twice for the same generation.

Two append-only CSVs live next to this script:

  runs.csv      one row per run attempt: what was launched (model, inputs, seed, settings)
  variants.csv  one row per (run, sequence): what came out

A run is identified by its run_key: a hash of method + model + the CONTENT of its input
files + seed + settings. A variant is identified by its seq_key: a hash of the normalised
heavy and light sequences. Standard library only, Python 3.8+.

  check-run RUNFLAGS                      exit 0 = new, exit 3 = this exact run is already logged
  add-run   RUNFLAGS --run-id ID          log a run before launching it (status starts "started")
  set-status --run-id ID --status S       S is done, failed or empty
  check DESIGNS.csv [--new-only OUT.csv]  exit 3 = some sequence is already in the ledger
  add   DESIGNS.csv --run-id ID           log a run's variants; re-run it to update verdicts

RUNFLAGS (the same flags go to check-run and add-run, so they cannot drift apart):
  --method M       what kind of generation: design, proteinmpnn, mutate, manual, ...
  --model X        weight image tag or tool + version, e.g. the full quay.io/... tag
  --input FILE     every file the run reads (repeatable). A .yaml spec's `input:` file is followed.
  --seed N         the seed, if the method has one
  --setting k=v    every other flag that changes what is generated (repeatable), spelled as in
                   the command without the dashes: --setting num-replicates=100. Never the output path.

DESIGNS.csv is the BASE VALIDATION table: sequence_id, heavy_sequence, light_sequence (plus kind
and design_description). Optional columns cdrh3, track_assigned and base_validation (PASS or FAIL)
are stored too; without base_validation a variant is logged as "untested".

Exit codes: 0 fine, 3 duplicate found (check-run, check), 1 usage or ledger error.
Rows are never edited or deleted: a changed status or verdict is a new row and the latest wins.
To repeat a run on purpose, add --setting rerun=2 and say why in --notes.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import getpass
import hashlib
import io
import json
import os
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUN_COLS = ["run_id", "run_key", "status", "date_utc", "owner", "method", "model", "inputs",
            "input_sha256", "seed", "settings", "output_dir", "notes"]
VAR_COLS = ["seq_key", "run_id", "sequence_id", "kind", "heavy_sequence", "light_sequence",
            "cdrh3", "base_validation", "track_assigned", "date_utc", "notes"]
STATUSES = ("started", "done", "failed", "empty")
BLOCKING = ("started", "done", "empty")  # a failed run may be retried once its cause is fixed
SEEN_EXIT = 3
YAML_INPUT = re.compile(r"^\s*input:\s*['\"]?([^\s'\"#]+)", re.M)


def die(msg):
    sys.exit(f"ledger: {msg}")


def now():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ---- the two CSV files ------------------------------------------------------------------

def read_ledger(path, cols):
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames != cols:
            die(f"{path.name} has columns {reader.fieldnames}, expected {cols}. Do not edit the header.")
        return list(reader)


def append_ledger(path, cols, rows):
    if not rows:
        return
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=cols, lineterminator="\n", extrasaction="ignore")
    fresh = not path.exists() or path.stat().st_size == 0
    if fresh:
        writer.writeheader()
    writer.writerows(rows)
    if not fresh:  # a hand edit may have dropped the final newline; never glue onto the last row
        with path.open("rb") as f:
            f.seek(-1, os.SEEK_END)
            if f.read(1) != b"\n":
                buf = io.StringIO("\n" + buf.getvalue())
    with path.open("a", newline="", encoding="utf-8") as f:
        f.write(buf.getvalue())


def paths(a):
    d = Path(a.ledger_dir)
    return d / "runs.csv", d / "variants.csv"


def latest_runs(rows):
    latest = {}
    for r in rows:
        latest[r["run_id"]] = r
    return latest


# ---- run fingerprint --------------------------------------------------------------------

def file_sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def collect_inputs(raw_paths):
    """Hash every input file, and the file a .yaml spec points at with `input:`."""
    found = {}
    for raw in raw_paths or []:
        p = Path(raw)
        if not p.exists():
            die(f"input not found: {raw} (run from the folder the model command runs in, e.g. apheris/rfd4)")
        files = sorted(x for x in p.rglob("*") if x.is_file()) if p.is_dir() else [p]
        for f in files:
            found[f.as_posix()] = file_sha(f)
            if f.suffix.lower() in (".yaml", ".yml"):
                for ref in YAML_INPUT.findall(f.read_text(encoding="utf-8", errors="replace")):
                    dep = f.parent / ref
                    if dep.is_file():
                        found[dep.as_posix()] = file_sha(dep)
    return found


def parse_settings(items):
    out = {}
    for item in items or []:
        if "=" not in item:
            die(f"--setting needs key=value, got {item!r}")
        key, value = item.split("=", 1)
        key = key.strip().lstrip("-")
        if key in ("output", "out", "output-dir"):
            die("leave the output path out of --setting: it does not change the result and would "
                "make every run look new")
        out[key] = value.strip()
    return out


def norm_model(model):
    model = model.strip().split("/")[-1]
    return model.split(":", 1)[1] if model.startswith("oslo-hackathon:") else model


def fingerprint(a):
    inputs = collect_inputs(a.input)
    spec = {"method": a.method.strip().lower(), "model": norm_model(a.model),
            "inputs": sorted(set(inputs.values())), "seed": (a.seed or "").strip(),
            "settings": parse_settings(a.setting)}
    key = hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()[:12]
    return key, spec, inputs


def owner_name(a):
    if a.owner or os.environ.get("LEDGER_OWNER"):
        return a.owner or os.environ["LEDGER_OWNER"]
    try:
        name = subprocess.run(["git", "config", "user.name"], capture_output=True, text=True,
                              timeout=5).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        name = ""
    return name or getpass.getuser()


def describe(r):
    return (f"{r['run_id']}  status={r['status']}  {r['date_utc'][:10]}  owner={r['owner']}  "
            f"output={r['output_dir'] or '?'}")


# ---- run commands -----------------------------------------------------------------------

def cmd_check_run(a):
    runs_csv, _ = paths(a)
    key, _, inputs = fingerprint(a)
    print(f"run_key {key}")
    for path, sha in sorted(inputs.items()):
        print(f"  input {sha[:12]}  {path}")
    hits = [r for r in latest_runs(read_ledger(runs_csv, RUN_COLS)).values() if r["run_key"] == key]
    if not hits:
        print("NEW     this exact run is not in the ledger; register it with add-run, then launch")
        return 0
    blocking = [r for r in hits if r["status"] in BLOCKING]
    for r in hits:
        print(f"{'SEEN' if r in blocking else 'FAILED':<7} {describe(r)}")
    if blocking:
        print("Do not launch this run. Reuse the output above, or change the seed or a setting "
              "if you want a different result.")
        return SEEN_EXIT
    print("Every earlier attempt failed. Retry only once the cause is fixed, and say what you changed in --notes.")
    return 0


def cmd_add_run(a):
    runs_csv, _ = paths(a)
    key, spec, _ = fingerprint(a)
    runs = latest_runs(read_ledger(runs_csv, RUN_COLS))
    if a.run_id in runs:
        die(f"run_id {a.run_id} already exists; use set-status to update it, or pick a new id")
    for r in runs.values():
        if r["run_key"] == key and r["status"] in BLOCKING:
            print(f"SEEN    {describe(r)}")
            die("this exact run is already in the ledger, not adding it. Reuse that output, or change "
                "the seed or a setting (add --setting rerun=2 for a deliberate repeat)")
    append_ledger(runs_csv, RUN_COLS, [{
        "run_id": a.run_id, "run_key": key, "status": a.status, "date_utc": now(),
        "owner": owner_name(a), "method": spec["method"], "model": spec["model"],
        "inputs": ";".join(a.input or []), "input_sha256": ";".join(h[:12] for h in spec["inputs"]),
        "seed": spec["seed"], "settings": ";".join(f"{k}={v}" for k, v in sorted(spec["settings"].items())),
        "output_dir": a.output or "", "notes": a.notes or ""}])
    print(f"logged run {a.run_id}  run_key={key}  status={a.status}")
    return 0


def cmd_set_status(a):
    runs_csv, _ = paths(a)
    run = latest_runs(read_ledger(runs_csv, RUN_COLS)).get(a.run_id)
    if run is None:
        die(f"unknown run_id {a.run_id}; register it with add-run first")
    run = dict(run, status=a.status, date_utc=now())
    if a.notes:
        run["notes"] = (run["notes"] + " | " if run["notes"] else "") + a.notes
    append_ledger(runs_csv, RUN_COLS, [run])
    print(f"run {a.run_id} is now {a.status}")
    return 0


# ---- variant commands -------------------------------------------------------------------

def clean_seq(s):
    return re.sub(r"\s+", "", s or "").upper()


def read_designs(path):
    p = Path(path)
    if not p.is_file():
        die(f"designs table not found: {path}")
    delim = "\t" if p.suffix.lower() == ".tsv" else ","
    with p.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter=delim)
        rows, cols = list(reader), list(reader.fieldnames or [])
    for need in ("sequence_id", "heavy_sequence"):
        if need not in cols:
            die(f"{path} needs a {need} column (found: {cols})")
    for r in rows:
        heavy, light = clean_seq(r["heavy_sequence"]), clean_seq(r.get("light_sequence"))
        if not heavy or not re.fullmatch(r"[A-Z]+", heavy) or not re.fullmatch(r"[A-Z]*", light):
            die(f"{r['sequence_id']}: heavy_sequence must be letters only and not empty, light_sequence letters only")
        r["_heavy"], r["_light"] = heavy, light
        r["_key"] = hashlib.sha256(f"{heavy}|{light}".encode()).hexdigest()[:16]
    return rows, cols, delim


def classify(designs, variants, exclude_run=None):
    """Return (status, design, note) per design. SEEN and DUP-IN-FILE are the blocking ones."""
    by_key, by_cdrh3 = {}, {}
    for v in variants:
        if v["run_id"] == exclude_run:
            continue
        by_key.setdefault(v["seq_key"], []).append(v)
        if v["cdrh3"]:
            by_cdrh3.setdefault(v["cdrh3"], []).append(v)
    first_in_file, out = {}, []
    for d in designs:
        key = d["_key"]
        if key in first_in_file:
            out.append(("DUP-IN-FILE", d, f"same sequence as {first_in_file[key]} above"))
            continue
        first_in_file[key] = d["sequence_id"]
        hits = by_key.get(key)
        if hits:
            known = [h["base_validation"] for h in hits if h["base_validation"] != "untested"]
            out.append(("SEEN", d, f"already generated: {len({h['run_id'] for h in hits})} run(s), first "
                        f"{hits[0]['run_id']} as {hits[0]['sequence_id']} ({hits[0]['date_utc'][:10]}), "
                        f"validation={known[-1] if known else 'untested'}"))
            continue
        cdrh3 = clean_seq(d.get("cdrh3"))
        h = by_cdrh3.get(cdrh3, [None])[0] if cdrh3 else None
        out.append(("NEW", d, f"new sequence, but CDRH3 {cdrh3} is already in the ledger "
                    f"({h['run_id']}:{h['sequence_id']})" if h else ""))
    return out


def print_report(results):
    for status, d, note in results:
        print(f"{status:<11} {d['sequence_id']:<26} {note}".rstrip())
    counts = {s: sum(1 for r in results if r[0] == s) for s in ("NEW", "SEEN", "DUP-IN-FILE")}
    print(f"{len(results)} designs: {counts['NEW']} new, {counts['SEEN']} already in the ledger, "
          f"{counts['DUP-IN-FILE']} repeated within the file")
    return counts


def cmd_check(a):
    _, vars_csv = paths(a)
    designs, cols, delim = read_designs(a.designs)
    results = classify(designs, read_ledger(vars_csv, VAR_COLS))
    counts = print_report(results)
    if a.new_only:
        with open(a.new_only, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter=delim, lineterminator="\n", extrasaction="ignore")
            w.writeheader()
            w.writerows(d for status, d, _ in results if status == "NEW")
        print(f"wrote the {counts['NEW']} new designs to {a.new_only}")
    return SEEN_EXIT if counts["SEEN"] or counts["DUP-IN-FILE"] else 0


def cmd_add(a):
    runs_csv, vars_csv = paths(a)
    run = latest_runs(read_ledger(runs_csv, RUN_COLS)).get(a.run_id)
    if run is None:
        die(f"unknown run_id {a.run_id}; register the run with add-run first")
    designs, _, _ = read_designs(a.designs)
    variants = read_ledger(vars_csv, VAR_COLS)
    print_report(classify(designs, variants, exclude_run=a.run_id))
    logged = {(v["run_id"], v["seq_key"]): v for v in variants}  # file order: the latest row wins
    rows = []
    for d in designs:
        verdict = (d.get("base_validation") or "untested").strip()
        verdict = "untested" if verdict.lower() == "untested" else verdict.upper()
        if verdict not in ("PASS", "FAIL", "untested"):
            die(f"{d['sequence_id']}: base_validation must be PASS, FAIL or untested, got {verdict!r}")
        row = {"seq_key": d["_key"], "run_id": a.run_id, "sequence_id": d["sequence_id"],
               "kind": d.get("kind") or "", "heavy_sequence": d["_heavy"], "light_sequence": d["_light"],
               "cdrh3": clean_seq(d.get("cdrh3")), "base_validation": verdict,
               "track_assigned": (d.get("track_assigned") or "").strip(), "date_utc": now(),
               "notes": d.get("design_description") or ""}
        old = logged.get((a.run_id, d["_key"]))
        if old and all(old[c] == row[c] for c in ("cdrh3", "base_validation", "track_assigned")):
            continue
        logged[(a.run_id, d["_key"])] = row
        rows.append(row)
    append_ledger(vars_csv, VAR_COLS, rows)
    print(f"logged {len(rows)} rows for {a.run_id} ({len(designs) - len(rows)} unchanged) in variants.csv")
    if run["status"] == "started":
        print(f"note: run {a.run_id} is still 'started'; once it has finished run "
              f"set-status --run-id {a.run_id} --status done")
    return 0


# ---- command line -----------------------------------------------------------------------

def build_parser():
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--ledger-dir", default=str(HERE),
                        help="folder holding runs.csv and variants.csv (default: next to this script)")

    runflags = argparse.ArgumentParser(add_help=False)
    runflags.add_argument("--method", required=True)
    runflags.add_argument("--model", required=True)
    runflags.add_argument("--input", action="append")
    runflags.add_argument("--seed")
    runflags.add_argument("--setting", action="append")

    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("check-run", parents=[common, runflags], help="is this exact run already logged?")
    s.set_defaults(func=cmd_check_run)

    s = sub.add_parser("add-run", parents=[common, runflags], help="log a run before launching it")
    s.add_argument("--run-id", required=True, help="unique, e.g. leo-20261114-her2novel-s1")
    s.add_argument("--status", choices=STATUSES, default="started")
    s.add_argument("--output", help="where the outputs go (a pointer only, not part of the run_key)")
    s.add_argument("--owner")
    s.add_argument("--notes")
    s.set_defaults(func=cmd_add_run)

    s = sub.add_parser("set-status", parents=[common], help="update a logged run's status")
    s.add_argument("--run-id", required=True)
    s.add_argument("--status", choices=STATUSES, required=True)
    s.add_argument("--notes")
    s.set_defaults(func=cmd_set_status)

    s = sub.add_parser("check", parents=[common], help="which designs are already in the ledger?")
    s.add_argument("designs")
    s.add_argument("--new-only", metavar="OUT", help="write only the designs not yet in the ledger")
    s.set_defaults(func=cmd_check)

    s = sub.add_parser("add", parents=[common], help="log a run's variants")
    s.add_argument("designs")
    s.add_argument("--run-id", required=True)
    s.set_defaults(func=cmd_add)
    return p


def main(argv=None):
    a = build_parser().parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
