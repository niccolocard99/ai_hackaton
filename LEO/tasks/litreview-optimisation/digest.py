#!/usr/bin/env python3
"""Deduplicate the cached AMASS results and write a triage digest.

Usage: python digest.py  (writes DIGEST.md next to the cache)

Reads cache/*.json, dedupes records by pmid/amassId, tags each with the queries
that returned it, and sorts by how many distinct queries hit it (a crude but
useful relevance signal) then by citation count.
"""
import json
import glob
import os
import re
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "cache")

records = {}
hits = defaultdict(set)

for path in sorted(glob.glob(os.path.join(CACHE, "*.json"))):
    name = os.path.basename(path)[:-5]
    core, _, slug = name.partition("__")
    try:
        data = json.load(open(path, encoding="utf-8"))["data"]
    except Exception as exc:  # noqa: BLE001
        print(f"skip {name}: {exc}")
        continue
    for rec in data:
        key = rec.get("pmid") or rec.get("amassId") or rec.get("publicationNumber")
        if not key:
            continue
        rec["_core"] = core
        records.setdefault(key, rec)
        hits[key].add(slug)


def year(rec):
    d = rec.get("publicationDate") or rec.get("createDate") or ""
    m = re.match(r"(\d{4})", str(d))
    return m.group(1) if m else "?"


def cites(rec):
    try:
        return int(rec.get("citationCount") or 0)
    except (TypeError, ValueError):
        return 0


ordered = sorted(records.items(), key=lambda kv: (-len(hits[kv[0]]), -cites(kv[1])))

by_core = defaultdict(list)
for key, rec in ordered:
    by_core[rec["_core"]].append((key, rec))

lines = [
    "# AMASS digest — trastuzumab optimisation literature",
    "",
    f"{len(records)} unique records from {len(glob.glob(os.path.join(CACHE, '*.json')))} queries.",
    "Sorted by number of distinct queries that returned the record, then citations.",
    "",
]

for core in sorted(by_core):
    lines += [f"## {core} ({len(by_core[core])} records)", ""]
    lines += ["| # queries | year | cites | title | id |", "|---|---|---|---|---|"]
    for key, rec in by_core[core]:
        title = (rec.get("title") or "").replace("|", "/").strip()
        retr = " **[RETRACTED]**" if rec.get("isRetracted") else ""
        lines.append(
            f"| {len(hits[key])} | {year(rec)} | {cites(rec)} | {title}{retr} | {key} |"
        )
    lines.append("")

out = os.path.join(HERE, "DIGEST.md")
open(out, "w", encoding="utf-8").write("\n".join(lines))
print(f"{len(records)} unique records -> {out}")
for core in sorted(by_core):
    print(f"  {core}: {len(by_core[core])}")
