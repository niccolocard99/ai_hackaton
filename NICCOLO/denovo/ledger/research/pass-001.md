# Research pass 001 — 2026-11-14

The literature pass that seeded this ledger. **Cache-only.**

- **scope:** campaign seed — import the de novo VHH priors from Leo's ledger and map each to a generation-arm parameter or a funnel threshold
- **targets:** TH-003, TH-004, TH-005, TH-006, TH-008 (imported); TH-009, TH-010, TH-011 (new); IDEA-001 … IDEA-010
- **searches:** 32 cache files read (585 unique records), **0 new**, 0 failed, 7 queued (see below)
- **cost this pass:** **$0.00**
- **cumulative:** 0 searches, $0.00, 0% of the ~2,000-search team budget

## Why zero spend

`AMASS_API_KEY` is **unset** and the team's AMASS credits are **unclaimed**, so
this pass made **zero network calls** by construction, not by thrift. Every
record cited anywhere in this ledger was read from the 32 cached responses
copied into [`../cache/`](../cache/) from
`LEO/tasks/experiment-ledger/cache/` (read-only for us; nothing was written
there). Queries we wanted and could not run are queued below rather than
quietly dropped.

## Sources read

| Source | What it gave |
|---|---|
| `LEO/tasks/experiment-ledger/` — PROMPT.md, THEORIES.md, IDEAS.md, literature/REVIEW.md, STATE.md, experiments.csv, validate.py, research_pass.py | The schema this ledger reuses, and the 8 theories / 7 ideas we imported from |
| `LEO/tasks/experiment-ledger/cache/` — 32 JSON responses, 585 unique records | Every PMID cited in this ledger |
| `apheris/rfd4/*.md`, `apheris/rfd4/Data/`, `apheris/rfd4/Settings/`, `apheris/challenge/` | Platform-side numbers: CDR length ranges, hotspot atoms, pass rates, `abguard` rules, per-design GPU cost |
| `github.com/Jerome-Alvarez/TEMPRO` (installed and run) | TEMPRO is obtainable; a real VHH returned Tm 75.55 C |

## What was imported, and what changed on import

Ids are preserved so a claim has one address across both ledgers. Three imports
needed **re-scoping**, which is recorded in each block's `history` and `notes`
rather than hidden:

| Id | Imported as | Change on import |
|---|---|---|
| TH-006 / IDEA-005 | The feasibility prior for the VHH route (PMID:37632647) | IDEA-005 promoted from "a hedge alongside the scFv" to **the whole campaign** |
| IDEA-003 | Thermostability proxy decision | **Closed favourably.** Leo's open question was format-dependent; we submit a VHH, so TEMPRO applies. Verified by installing and running it. |
| TH-003 / IDEA-002 | Arg cap in CDRs | Numeric cap added as OUR JUDGEMENT — the source has no number |
| TH-004 | No rigidity / B-factor filter | Scope limit added: this does not cover pLDDT or PAE, which we do rank on |
| TH-005 / IDEA-001 | Binding risk in CDRs, stability elsewhere | **Re-scoped.** Framework redesign is undocumented on RFD4-Proteina, so "framework design" became "framework **choice**" (caplacizumab, kept fixed). The stability half of Leo's idea is not available to us. |
| TH-008 / IDEA-004 | A contact is not a hot spot | **Re-scoped.** We take the 4 cited `C_HOT` atoms to aim generation; we drop the per-residue-energetics plan, which we have no tooling for. |
| IDEA-007 | Patent / freedom-to-operate | Blocked on unclaimed credits. One zero-cost action salvaged: add PMID:37632647's published VHH to the novelty-comparison set. |
| PMID:41193805 | New: TH-009 (CDR3 length), TH-011 (design counts), IDEA-008 (mine the fulltext) | See the blocker below |

## Blocker found: the cache is abstract-only

The brief asked us to mine **PMID:41193805** ("Atomically accurate de novo
design of antibodies with RFdiffusion", Nature 2025) for CDR length ranges,
hotspot counts, designs per target, filters and reported success rates.

**It cannot be done from the cache.** A biomedcore record carries metadata plus
`abstract` and nothing else — there is a `hasFulltext: True` flag but no
fulltext body on the record. So the paper's methods numbers are simply not
present locally. What the abstract does give is recorded in IDEA-008, and the
two theories that needed length and count numbers (TH-009, TH-011) were built
from other cited sources instead: PMID:19062174 for H3 prediction accuracy by
loop length, and the Apheris docs for pass rates and per-design cost. The paper
is open access as **PMC12727541**; reading its methods by hand is IDEA-008.

## Queued queries — wanted, not run, zero spent

Run these when `AMASS_API_KEY` is set and credits are claimed. Written down
here so the gap is visible instead of forgotten.

| # | Core | Query | Which decision it would settle |
|---|---|---|---|
| Q1 | biomedcore | `RFdiffusion de novo antibody VHH CDR length hotspot designs tested success rate` | TH-009, TH-011 — the PMID:41193805 recipe numbers the abstract withholds |
| Q2 | biomedcore | `nanobody CDR3 length distribution camelid repertoire residues` | TH-009 — a natural VHH CDR3 length prior, which we do not have; our 12-16 ceiling is currently judgement |
| Q3 | biomedcore | `antibody developability sequence liability threshold deamidation isomerization glycosylation CDR` | TH-010 — would replace four OUR JUDGEMENT thresholds with cited ones |
| Q4 | biomedcore | `nanobody aggregation hydrophobic patch surface prediction threshold` | TH-010 — the hydrophobic-patch cap, currently our crudest number |
| Q5 | biomedcore | `de novo binder design experimental validation hit rate nanomolar fraction` | TH-011 — an independent hit rate not borrowed from upstream |
| Q6 | patentcore | `anti-HER2 VHH single domain antibody sequence claim` | IDEA-007 — freedom to operate on the VHH route specifically |
| Q7 | patentcore | `nanobody HER2 domain IV epitope claim` | IDEA-007 |

## Two real bugs in Leo's `research_pass.py`

Both claims in the brief were checked against the source and **both are real**.
Not fixed in Leo's folder (read-only for us). Detail, reproduction and the fix
shape are in [`../LIT_DIGEST.md`](../LIT_DIGEST.md) under "Bugs in the upstream
research_pass.py".

1. **trialcore records are scored on fields that do not exist.** Confirmed: all
   20 records in `cache/trialcore__*.json` have `briefTitle` (20/20) and
   `briefSummary` (17/20), and **zero** have `title` or `abstract`.
   `score_record` builds its blob from `title` + `abstract`, so every trialcore
   record scores 0 and can never surface. This is the "empty titles" blocker in
   Leo's `STATE.md`.
2. **The cache key collides.** `cache_path` is `{core}__{slug(query)[:60]}`.
   Confirmed: **12 of the 32** existing cache filenames are already at exactly
   60 characters, so the truncation is live, and two queries differing only
   after character 60 provably resolve to one file. Worse than the brief
   stated: `limit` is **already** a request-determining parameter that is not
   in the key, so the bug does not wait for filters to be added.

## Proposed updates — applied

| Target | Change | Justified by |
|---|---|---|
| TH-003, TH-004, TH-005, TH-006, TH-008 | Imported at Leo's confidence, ids preserved, `implies` added | pass-001 |
| TH-009, TH-010, TH-011 | Created | PMID:19062174; PMID:24492299; Apheris docs |
| IDEA-001, 002, 004, 005, 007 | Imported, three re-scoped | pass-001 |
| IDEA-003 | `lit-checked` → **`done`** | PMID:39154093 + PMID:41564239 + TEMPRO installed and run |
| IDEA-008, 009, 010 | Created | pass-001 |
| `LIT_DIGEST.md` | Written | pass-001 |

## Experiment rows

**None, deliberately.** `experiments.csv` carries its header and no rows.
Hypotheses are written **before** runs, so the first `EXP-001` row gets written
when the first generation arm is launched, not now.
