# Ledger — de novo VHH binders to HER2 domain IV

This campaign's experiment ledger. Same schema as
[`../../../LEO/tasks/experiment-ledger/`](../../../LEO/tasks/experiment-ledger/),
whose cycle is defined in that folder's `PROMPT.md` and still governs: pick an
idea, literature-check it, **write the hypothesis row before the run**, run,
re-check, update.

**Run `python -I validate.py` before every commit.** It exits 1 on any error.

## Files

| File | What it holds |
|---|---|
| [`experiments.csv`](experiments.csv) | One append-only row per thing we tried. Header only right now — hypotheses are written before runs. Correct history by adding a superseding row, never by editing. |
| [`THEORIES.md`](THEORIES.md) | What we believe, each with an `implies` field naming the generation-arm parameter or funnel threshold it sets |
| [`IDEAS.md`](IDEAS.md) | The backlog |
| [`LIT_DIGEST.md`](LIT_DIGEST.md) | The evidence behind every threshold, one section per design decision, each number marked CITED or OUR JUDGEMENT |
| [`research/pass-001.md`](research/pass-001.md) + `.json` | The seeding literature pass. Cache-only, $0.00, 7 queries queued. |
| [`research/spend.json`](research/spend.json) | Cumulative AMASS spend. Currently 0 searches, $0.00. |
| [`cache/`](cache/) | The 32 cached AMASS responses (585 unique records), copied from Leo's ledger so this one is self-contained and reproducible offline |
| `experiments/EXP-nnn.md` | The long form per experiment: expected, what the literature said, what happened, **why** |

## Id policy

Append-only. Ids imported from Leo's ledger **keep their original ids** so a
claim has one address across both ledgers, and carry `imported_from`. Ids new
to this campaign start past Leo's highest: **TH-009+**, **IDEA-008+**.
`validate.py` requires every theory and idea to carry either `imported_from` or
`new_in`.

## Appending an experiment row

The hypothesis goes in **before** the run. Two steps:

```bash
cd /home/hackathon/hackathon-team-3/NICCOLO/denovo/ledger

# 1. BEFORE the run: append the row with hypothesis filled and verdict EMPTY.
#    An empty verdict is legal; a verdict with an empty lit_checked is not.
python3 -I - <<'PY'
import csv
row = {
    "exp_id": "EXP-001", "date": "2026-11-14",
    "theory_ids": "TH-006;TH-009;TH-011", "idea_id": "IDEA-005",
    "lit_checked": "pass-001",
    "hypothesis": "Falsifiable sentence, written before the run.",
    "what_we_did": "", "inputs": "", "outputs": "", "metrics": "",
    "verdict": "", "learned": "", "cost": "", "supersedes": "",
}
with open("experiments.csv", "a", newline="", encoding="utf-8") as fh:
    csv.DictWriter(fh, fieldnames=[
        "exp_id","date","theory_ids","idea_id","lit_checked","hypothesis",
        "what_we_did","inputs","outputs","metrics","verdict","learned",
        "cost","supersedes"]).writerow(row)
PY
python -I validate.py

# 2. AFTER the run: supersede it. Never edit the row in place.
#    Append a new EXP-002 row with supersedes=EXP-001, metrics and verdict filled,
#    and write experiments/EXP-002.md.
```

Legal `verdict` values: `worked`, `failed`, `inconclusive`, `aborted`,
`closed-by-literature`. Legal `lit_checked`: `pass-NNN` (the file must exist)
or `skipped:<reason>`.

## Constraints this campaign runs under

- **Cache-only literature.** `AMASS_API_KEY` is unset and credits are
  unclaimed, so no new AMASS call has been made. Wanted queries are **queued**
  in `research/pass-001.md`, not silently dropped.
- **De novo track**, so the novelty rule that binds is `abguard`'s "< 70%
  sequence identity to any known antibody" — not the CDRH3-vs-trastuzumab rule,
  which belongs to the optimisation track.
- Hard `abguard` gates: no unpaired cysteines, <= 250 aa total, `immunum`
  numbering confidence >= 0.1.
- Measured and predicted never share a cell. Every number says CITED or
  OUR JUDGEMENT.
