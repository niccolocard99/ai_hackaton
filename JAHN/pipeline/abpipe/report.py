"""Reports for the orchestrated loop: one small report per step, one major report per iteration.

STEPS.md (one file per round) grows while the round runs: every module completion the runner logs is
appended with a timestamp, so the funnel is readable live. REPORT.md (one file per round) is written
when the round finishes: the decision that shaped it, the funnel with rejection reasons, the score
statistics, and the selected list — everything the pitch's "clear reasoning" criterion wants, one
file per iteration.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import polars as pl


class StepLog:
    """log= replacement that prints, and appends every runner message to <round>/STEPS.md."""

    def __init__(self, rd: Path, t0: float, log=print):
        self.rd, self.t0, self.log = rd, t0, log
        self.path = rd / "STEPS.md"
        if not self.path.exists():
            self.path.write_text(f"# Steps — {rd.name}\n\n")
            self._f(f"round starts {time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime())} UTC")

    def _f(self, line: str) -> None:
        with open(self.path, "a") as fh:
            fh.write(line + "\n")

    def __call__(self, msg: str) -> None:
        self.log(msg)
        m = msg.strip()
        if m and not m.startswith(("run ", "manifest:", "=== ")):
            self._f(f"- [{time.time() - self.t0:7.1f}s] {m}")


def _funnel(man: dict) -> list[str]:
    lines = ["| stage | module | rows in | rows out | dropped | seconds |",
             "|---|---|---|---|---|---|"]
    gen = man["stages"].get("generate", {})
    for m in gen.get("modules", []):
        lines.append(f"| generate | {m['name']} | {m.get('n_requested', '')} | {m['rows_out']} | — | "
                     f"{m['seconds']}{' (cached)' if m['cached'] else ''} |")
    flt = man["stages"].get("filter", {})
    for m in flt.get("modules", []):
        lines.append(f"| filter | {m['name']} | {m.get('rows_in', '')} | {m['rows_out']} | "
                     f"{m.get('failed', 0)} | {m['seconds']}{' (cached)' if m['cached'] else ''} |")
    for m in man["stages"].get("score", {}).get("modules", []):
        lines.append(f"| score | {m['name']} | {m.get('rows_out', '')} | {m.get('rows_out', '')} | — | "
                     f"{m['seconds']}{' (cached)' if m['cached'] else ''} |")
    for m in man["stages"].get("rank", {}).get("modules", []):
        lines.append(f"| rank | {m['name']} | {m['rows_out']} | {m['rows_out']} | — | {m['seconds']} |")
    return lines


def _rejections(man: dict, top: int = 4) -> list[str]:
    rej = man["stages"].get("filter", {}).get("rejections") or {}
    if not rej:
        return []
    lines = ["", "Why candidates were dropped (per filter, top reasons):"]
    for name, r in rej.items():
        if not r.get("count"):
            continue
        reasons = ", ".join(f"{x['reason']} ×{x['count']}" for x in r["top_reasons"][:top])
        lines.append(f"- {name}: {r['count']} dropped — {reasons}")
    return lines


def _selected_table(rd: Path, k: int = 12) -> list[str]:
    sel = pl.read_parquet(rd / "selected.parquet")
    score_cols = [c for c in sel.columns if c.startswith("score__")]
    info_cols = [c for c in ("info__track", "info__cdrh3_min_distance", "info__strategy",
                             "info__parent_score") if c in sel.columns]
    keep = ["sequence_id", "origin"] + score_cols + info_cols if "origin" in sel.columns \
        else ["sequence_id"] + score_cols + info_cols
    head = sel.head(k).select([c for c in keep if c in sel.columns])
    lines = ["", f"Top {min(k, sel.height)} of {sel.height} selected:", ""]
    lines.append("| " + " | ".join(c.replace("score__", "").replace("info__", "")
                                   for c in head.columns) + " |")
    lines.append("|" + "---|" * head.height.__class__(len(head.columns)))
    for r in head.iter_rows(named=True):
        cells = []
        for c in head.columns:
            v = r[c]
            cells.append(f"{v:.3f}" if isinstance(v, float) else str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return lines


def iteration_report(rd: Path, entry: dict, man: dict | None = None) -> None:
    """The major report for one finished iteration: <round>/REPORT.md."""
    r = entry.get("round", "?")
    L = [f"# Iteration report — round {int(r):02d}", "", f"- finished: {entry.get('utc', '')}"]
    if entry.get("wall_seconds") is not None:
        L.append(f"- wall time: {entry['wall_seconds']} s")
    gpu = (entry.get("probe") or {}).get("gpu")
    if gpu:
        L.append(f"- machine while deciding: GPU {gpu['util_pct']}% ({gpu['mem_used_mb']}/"
                 f"{gpu['mem_total_mb']} MB), disk {entry['probe'].get('disk_free_gb')} GB free")
    b = entry.get("budget") or {}
    if b.get("deadline_remaining_s") is not None:
        L.append(f"- budget: {b['deadline_remaining_s'] / 60:.0f} min left; a round costs about "
                 f"{(b.get('est_round_seconds') or 0) / 60:.0f} min")
    L += ["", "## The decision that shaped this round", f"- policy: {entry.get('source', '')}",
          f"- patch: {json.dumps(entry.get('patch') or {}, sort_keys=True)}"]
    for why in entry.get("reasoning") or []:
        L.append(f"- why: {why}")
    if man is not None:
        L += ["", "## Funnel", *_funnel(man), *_rejections(man)]
    if entry.get("generated") is not None:
        L += ["", f"## Outcome: {entry.get('generated')} generated -> {entry.get('survivors')} survive filters "
              f"-> {entry.get('selected')} selected ({entry.get('carried_over', 0)} carried over, "
              f"{entry.get('new_in_selected', 0)} new)"]
        v = entry.get("validate")
        if v:
            L.append(f"- BASE VALIDATION on the shortlist: {v['in']} in -> {v['out']} out "
                     f"({v['dropped']} dropped by the official checkers)")
            for rej in (v.get("top_rejections") or [])[:5]:
                L.append(f"    - {rej['sequence_id']}: {rej['reason']}")
        if entry.get("selected_medians"):
            L.append(f"- selected medians: {json.dumps(entry['selected_medians'], sort_keys=True)}")
        if entry.get("mean_pairwise_hamming") is not None:
            L.append(f"- mean pairwise Hamming of the selection: {entry['mean_pairwise_hamming']}")
        L += _selected_table(rd)
    if entry.get("stop"):
        L += ["", f"**STOP after this round: {entry['stop_reason']}**"]
    (rd / "REPORT.md").write_text("\n".join(L) + "\n")
