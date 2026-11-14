# Literature review — trastuzumab optimisation

**Branch:** `LEO-litreview-optimisation` · **Started:** 2026-11-14 · **Agent:** Claude Opus 5

**Read both agents files:** `../../../AGENTS.md` (root, shared) ☑ · `../../AGENTS.md` (Leo's) ☑

## What this task is

Gather everything published on engineering trastuzumab and turn it into actionable,
position-level design guidance for the three scored properties (binding affinity,
epitope preservation, thermostability). Feeds the 40% "judgment" criterion — the
reasoning behind the molecules we pick — and closes the open item "no thermostability
predictor chosen".

Scope is the **optimisation** track. A de novo counterpart pass is not yet written.

## Right now

Five synthesis sub-agents running in parallel over the cached AMASS results, one per
topic area (affinity, epitope, thermostability, prediction tools, patents/clinical).
Started 2026-11-14 ~11:10 UTC.

## Done so far

| What | Verified how | Commit |
|---|---|---|
| Review prompt written | `PROMPT.md`, reviewed against root AGENTS.md constraints | (this commit) |
| AMASS API reachable, key valid | `HTTP 200` + `x-ratelimit-remaining: 59` on a 5-result test query | (this commit) |
| 21 AMASS queries fetched and cached | `amass_search.sh` reported `fetched=21 cached=0 failed=0`, cost ~$1.05 | (this commit) |
| 385 unique records deduped | `digest.py` -> `DIGEST.md`; biomedcore 327, patentcore 38, trialcore 20 | (this commit) |

## Blocked

| Item | Blocker | Needed to unblock |
|---|---|---|
| Branch naming per `LEO/AGENTS.md` | Git cannot create `refs/heads/LEO/<task>` while `refs/heads/LEO` exists — a ref cannot be both file and directory. Used `LEO-litreview-optimisation` instead | Leo to decide: rename the `LEO` branch (e.g. `LEO-main`) so `LEO/<task>` works, or change `LEO/AGENTS.md` to the hyphen convention |
| Full-text reading | Only abstracts are cached; `hasFulltext` is set on some records | Fetch full text for the shortlisted papers via AMASS `include=` or publisher |
| Trastuzumab's own measured Tm baseline | Not yet confirmed from literature | Check the shortlisted thermostability papers; the organisers' measured baseline is authoritative |

## Sub-agents

| Sub-agent | Folder | What it is doing |
|---|---|---|
| affinity | (returns to this folder) | Binding affinity + CDR/CDRH3 affinity maturation |
| epitope | (returns to this folder) | Trastuzumab paratope / HER2 domain IV epitope, preservation risk |
| thermostability | (returns to this folder) | Stabilising mutations, affinity-stability trade-off |
| tools | (returns to this folder) | Thermostability predictors, solubility/expression filters, CDRH3 loop rules |
| priorart | (returns to this folder) | Patent crowding + clinical unmet need for the pitch |

Sub-agents here synthesise from the cache and return text; they were not given write
scope, so all their output is consolidated by the parent into `FINDINGS.md`.

## Next, in order

1. Consolidate the five syntheses into `FINDINGS.md` with the findings table, ranked
   recommendations, contradictions and sources.
2. Pick a thermostability predictor and record the reasoning (closes an open item).
3. Pull full text for the shortlisted papers and verify every number used in a design
   decision against the paper, not the abstract.
4. Write the de novo counterpart prompt and run that pass.
5. Translate the findings into concrete CDRH3 design specs for RFD4-Proteina.

## Log

| Date | Event |
|---|---|
| 2026-11-14 | Branch created (as `LEO-litreview-optimisation`; `LEO/<task>` impossible, see Blocked) |
| 2026-11-14 | `LEO/amass_key.txt` found untracked but unignored; added `LEO/.gitignore` rules on the `LEO` branch before any other work |
| 2026-11-14 | AMASS API spec retrieved from the public skill repo; base `https://api.amass.tech/api/v1`, bearer auth, 60 req/60 s, $0.05/20 results |
| 2026-11-14 | 21 queries across biomedcore/patentcore/trialcore cached, ~$1.05 of the $100 credit. Re-runs are free (cache hits) |
| 2026-11-14 | 385 unique records; top hits confirm relevance (cryo-EM HER2-trastuzumab-pertuzumab 31042744, CDR3 randomisation for HER2 26225765, framework-mutation stability 40653772) |
| 2026-11-14 | Five synthesis sub-agents launched in parallel |
