#!/usr/bin/env python3
"""Integrity checks for the de novo VHH / HER2 experiment ledger.

ADAPTED COPY of LEO/tasks/experiment-ledger/validate.py. It points at THIS
folder's files (HERE-relative, so nothing of Leo's is read or written) and adds
three checks that matter for this campaign:

  * an imported theory or idea must carry `imported_from` (provenance back to
    Leo's ledger, since we preserve his ids);
  * every theory must carry an `implies` field -- a theory that does not say
    what it changes about a generation arm or a funnel threshold is not
    actionable and does not belong here;
  * LIT_DIGEST.md must exist.

Run before every commit:

  python -I validate.py            # check, exit 1 on any error
  python -I validate.py --quiet    # errors only

Inherited checks: the CSV header is exact; exp_ids are unique and well formed;
verdicts are legal; every theory_ids / idea_id / supersedes / lit_checked
reference resolves; no row has a verdict but an empty lit_checked (the
literature step is not optional); every experiment has its long-form file;
theory and idea blocks parse and carry the required fields.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(HERE, "experiments.csv")
THEORIES = os.path.join(HERE, "THEORIES.md")
IDEAS = os.path.join(HERE, "IDEAS.md")
RESEARCH = os.path.join(HERE, "research")
DIGEST = os.path.join(HERE, "LIT_DIGEST.md")

HEADER = [
    "exp_id", "date", "theory_ids", "idea_id", "lit_checked", "hypothesis",
    "what_we_did", "inputs", "outputs", "metrics", "verdict", "learned",
    "cost", "supersedes",
]
VERDICTS = {"worked", "failed", "inconclusive", "aborted", "closed-by-literature"}
CONFIDENCES = {"high", "medium", "low", "retired"}
IDEA_STATUSES = {"open", "lit-checked", "running", "done", "dropped"}

errors: list[str] = []
warnings: list[str] = []


# [ \t]* not \s*: \s matches newlines, so an empty field
# ("- **evidence_against:**") would swallow the line below it.
FIELD_RE = re.compile(r"^[ \t]*-[ \t]+\*\*([a-z_]+):\*\*[ \t]*(.*)$", re.M)


def parse_fields(body: str) -> dict[str, str]:
    """Field lines, each absorbing any indented continuation beneath it.

    `history` and `notes` hold their value as indented sub-bullets, so a
    single-line match would read them as empty.
    """
    fields: dict[str, str] = {}
    matches = list(FIELD_RE.finditer(body))
    for i, fm in enumerate(matches):
        stop = matches[i + 1].start() if i + 1 < len(matches) else len(body)
        tail = body[fm.end():stop]
        extra = " ".join(
            line.strip().lstrip("-").strip()
            for line in tail.splitlines() if line.strip()
        )
        value = fm.group(2).strip()
        fields[fm.group(1)] = f"{value} {extra}".strip() if extra else value
    return fields


def parse_blocks(path: str, prefix: str) -> dict[str, dict]:
    if not os.path.exists(path):
        errors.append(f"missing {os.path.basename(path)}")
        return {}
    text = open(path, encoding="utf-8").read()
    pattern = re.compile(rf"^###\s+({prefix}-\d+)\s*(?:--|—|-)\s*(.+?)$", re.M)
    matches = list(pattern.finditer(text))
    out: dict[str, dict] = {}
    for i, m in enumerate(matches):
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        body = text[m.end():end]
        fields = parse_fields(body)
        if m.group(1) in out:
            errors.append(f"duplicate id {m.group(1)} in {os.path.basename(path)}")
        out[m.group(1)] = {"claim": m.group(2).strip(), **fields}
    if not out:
        errors.append(f"no {prefix}-* blocks parsed from {os.path.basename(path)}")
    return out


def refs(value: str) -> list[str]:
    return [t.strip() for t in re.split(r"[;,]", value or "") if t.strip()]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    theories = parse_blocks(THEORIES, "TH")
    ideas = parse_blocks(IDEAS, "IDEA")

    for tid, t in theories.items():
        if t.get("confidence") not in CONFIDENCES:
            errors.append(f"{tid}: confidence {t.get('confidence')!r} not in {sorted(CONFIDENCES)}")
        for required in ("falsified_by", "history", "query_terms"):
            if not t.get(required):
                errors.append(f"{tid}: empty required field {required}")
        # Campaign-specific: a theory must say what it changes.
        if not t.get("implies"):
            errors.append(f"{tid}: empty required field implies — a theory must "
                          "name the generation-arm parameter or funnel threshold "
                          "it sets, or it is not actionable")
        # Campaign-specific: preserved ids must name where they came from.
        if not t.get("imported_from") and not t.get("new_in"):
            errors.append(f"{tid}: needs imported_from (an id preserved from "
                          "Leo's ledger) or new_in (first raised in this campaign)")
    for iid, i in ideas.items():
        if i.get("status") not in IDEA_STATUSES:
            errors.append(f"{iid}: status {i.get('status')!r} not in {sorted(IDEA_STATUSES)}")
        if i.get("status") == "dropped" and not i.get("notes"):
            errors.append(f"{iid}: dropped without a reason in notes")
        for ref in refs(i.get("tests", "")):
            if ref not in theories:
                errors.append(f"{iid}: tests unknown theory {ref}")
        if not i.get("imported_from") and not i.get("new_in"):
            errors.append(f"{iid}: needs imported_from or new_in")

    if not os.path.exists(DIGEST):
        errors.append("missing LIT_DIGEST.md — the literature digest is the "
                      "evidence behind every threshold in THEORIES.md")

    if not os.path.exists(CSV_PATH):
        errors.append("missing experiments.csv")
        return report()

    with open(CSV_PATH, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))
    with open(CSV_PATH, encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh)
        actual_header = next(reader, [])
    if actual_header != HEADER:
        errors.append(f"CSV header mismatch.\n  expected: {HEADER}\n  actual:   {actual_header}")

    pass_ids = {
        re.search(r"pass-(\d+)\.md$", p).group(1)
        for p in (os.listdir(RESEARCH) if os.path.isdir(RESEARCH) else [])
        if re.search(r"pass-(\d+)\.md$", p)
    }

    seen: set[str] = set()
    for n, row in enumerate(rows, start=2):
        eid = (row.get("exp_id") or "").strip()
        where = f"row {n} ({eid or 'no id'})"
        if not re.fullmatch(r"EXP-\d{3}", eid):
            errors.append(f"{where}: exp_id must look like EXP-001")
        if eid in seen:
            errors.append(f"{where}: duplicate exp_id")
        seen.add(eid)

        verdict = (row.get("verdict") or "").strip()
        if verdict and verdict not in VERDICTS:
            errors.append(f"{where}: verdict {verdict!r} not in {sorted(VERDICTS)}")
        if verdict and not (row.get("lit_checked") or "").strip():
            errors.append(f"{where}: has a verdict but lit_checked is empty — "
                          "the literature step is not optional")

        lit = (row.get("lit_checked") or "").strip()
        for token in refs(lit):
            if token.startswith("skipped:"):
                continue
            m = re.fullmatch(r"pass-(\d+)", token)
            if not m:
                errors.append(f"{where}: lit_checked token {token!r} is neither "
                              "pass-NNN nor skipped:<reason>")
            elif m.group(1) not in pass_ids:
                warnings.append(f"{where}: lit_checked cites {token}, which has no file yet")

        for ref in refs(row.get("theory_ids", "")):
            if ref not in theories:
                errors.append(f"{where}: unknown theory {ref}")
        idea = (row.get("idea_id") or "").strip()
        if idea and idea not in ideas:
            errors.append(f"{where}: unknown idea {idea}")
        sup = (row.get("supersedes") or "").strip()
        if sup and sup not in seen and sup not in {r.get("exp_id") for r in rows}:
            errors.append(f"{where}: supersedes unknown experiment {sup}")

        long_form = os.path.join(HERE, "experiments", f"{eid}.md")
        if eid and not os.path.exists(long_form):
            warnings.append(f"{where}: no experiments/{eid}.md")

    return report()


def report() -> int:
    for w in warnings:
        print(f"warn  {w}")
    for e in errors:
        print(f"ERROR {e}")
    if errors:
        print(f"\n{len(errors)} error(s), {len(warnings)} warning(s) — ledger is invalid")
        return 1
    print(f"ledger OK ({len(warnings)} warning(s))")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
