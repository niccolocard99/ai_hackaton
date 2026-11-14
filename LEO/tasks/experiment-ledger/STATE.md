# experiment-ledger

**Branch:** `LEO-experiment-ledger` · **Started:** 2026-11-14 · **Agent:** Claude Opus 5

**Read both agents files:** [`../../../AGENTS.md`](../../../AGENTS.md) (root, shared) ☑ · [`../../AGENTS.md`](../../AGENTS.md) (Leo's) ☑

## Counts

| | |
|---|---|
| Experiments | 9 (`worked` 7, `failed` 1, `inconclusive` 1) — **3 of them ours** |
| Theories | 18 active + 2 retired — `low` 1, `medium` 6, `high` 11 |
| Ideas | 14 open/lit-checked, 3 done |
| Research passes | 7 (step B blocked: no credits) |
| AMASS spend | 14 searches, $0.70 — **team credits now exhausted**, see Blocked |
| Published | PR #2 **merged** into `main` at `27b6982`; later work goes out as a new PR |

## What this task is

The bookkeeping layer for the design work: an append-only record of every
experiment, a falsifiable theory set that literature and results update, an idea
backlog, and an AMASS autoresearch loop that runs as a required step of the
cycle rather than a one-off review. It exists because 40% of the pitch score is
judgment and 20% is a reproducible repo — both are scored on written reasoning,
not on the designs alone. See [`PROMPT.md`](PROMPT.md) for the full cycle.

## Right now

**Autoresearch loop running** — one pass every 15 minutes (at :08, :23, :38,
:53) until 18:00 local, pushing to `fork` (`leossgh/ai_hackaton`), which feeds
**open PR #8** against `niccolocard99/ai_hackaton:main`. Pushes go to `fork`,
never `origin`: `leossgh` has no write access upstream.

**Step B (AMASS) is blocked — credits exhausted.** Step A, the repo scan, is
free and is where every substantive result has come from: EXP-002 through
EXP-005 were all found by reading teammates' commits. Passes continue on Step A
alone until credits are topped up.

## Done so far

| What | Verified how | Commit |
|---|---|---|
| Ledger skeleton: CSV, theories, ideas, standing review, state | `python -I validate.py` passes | initial |
| 8 theories seeded from 30 cached AMASS queries, every claim citing a PMID | abstracts read directly from `cache/`, not from titles | initial |
| EXP-001 backfilled from `LEO/STATE.md` with its five reasons for being unrankable | cross-checked against `runs/her2/scored-boltz2` figures in that file | initial |
| `research_pass.py` runs a scoped pass, caches, tracks spend, writes pass files | pass-001 executed at `--max-new 0`, zero spend | initial |
| Field parser bug fixed: `\s*` swallowed the line after an empty field | 4 of 8 theories had silently lost `query_terms`; re-ran and all 8 parse | initial |
| **EXP-009: cross-checked this ledger against the full SPR-768 table** — overturned TH-019's hit rate and inverted IDEA-008's recipe | `python -I validate.py` passes; scripts in `NICCOLO/eda/p15`–`p17` | this pass |

## Blocked

| Item | Blocker | Needed to unblock |
|---|---|---|
| **AMASS research (Step B)** | HTTP 403 `"You are out of API credits"`. This ledger spent 14 searches / $0.70, so the ~2,000-search team budget went elsewhere — the sibling literature task ran ~30 queries plus a full-text corpus fetch, which costs far more per record | Top up at <https://platform.amass.tech/settings/credits>, or accept repo-scan-only passes. `research_pass.py` now stops on the first credits 403 instead of retrying every query |
| Theory confidence above `medium` | No measured result of our own yet; everything rests on literature | A design run with a recorded hypothesis (IDEA-006 is cheapest) |
| TH-007 / IDEA-003 (thermostability predictor) | All candidates are trained on proxies (TS50), not the Tm Adaptyv measures | A dedicated pass, then a decision |
| EXP-001 promotion from `inconclusive` | Single Boltz-2 seed, 2 of 4 designs refolded, novelty unchecked, `abguard` never run | IDEA-006 plus base validation |
| trialcore evidence | 20 cached records returned empty titles | Re-query or drop the core |

## Next, in order

0. **Top up AMASS credits** (or accept repo-scan-only passes). Step B cannot run at all until then — though note Step A has produced every substantive result so far.
1. **IDEA-011 — switch the ranker.** Adopt the sequence pre-filter (`NICCOLO/eda/prefilter.py`, branch `origin/nick`) as the primary ranking for optimisation-track designs. Costs no GPU, runs on a laptop, and replaces a metric EXP-003 shows is near-useless (Boltz-2 ipTM AUROC 0.645) with one at Spearman +0.622 against measured KD. Highest-value action available.
2. ~~IDEA-017 — pull the full SPR-768 table~~ **DONE (EXP-009).** It cost us two claims: the legal-band hit rate is **1 of 131 (0.8%)**, not 1 of 21, and **0 of 131** against trastuzumab's best replicate; and the per-position pattern behind IDEA-008 was a base-rate artefact.
3. **IDEA-008 — hand the design to whoever has the GPU**, with the recipe **revised by EXP-009**: `G100 + G101 + A106 + D108`, 63.6% identity. **Do not run the old recipe** — M107 (q=0.040) and F104 (q=0.028) measurably *hurt* affinity on 131 legal variants, A106/D102/D108 are inert, and G100/G101 are the only two positions where mutating helps. M107 never contacts HER2 yet is the costliest to mutate, so non-contact means tolerated, not free.
3. **IDEA-005 — add a VHH** alongside the scFv designs. TH-012: the only format where both no-proxy properties (Tm, structure accuracy) are covered.
4. Each pass: scan the repo for new committed results, `--scope` the weakest theory, apply, commit.
5. IDEA-010 — energetics for Y109, the one assumption IDEA-008 rests on and the one position SPR-768 cannot speak to. pass-003 proved the literature will not answer it.
6. IDEA-012 — make the decoy control a standing rule before any new metric is trusted.
4. IDEA-001 (CDRH3 for binding + framework for Tm) — still the only idea attacking two scored properties at once, and TH-005 is unaffected by EXP-002.
5. Resolve IDEA-003 (thermostability predictor); answer now depends on whether we submit a VHH.

## Log

| Date | Event |
|---|---|
| 2026-11-14 | Branch created from `LEO` in an isolated git worktree, because the litreview agent holds the primary working tree and a branch switch there would have redirected its commits |
| 2026-11-14 | Seeded 8 theories and 7 ideas from the 30 cached AMASS queries; abstracts read directly, since `DIGEST.md` carries titles only |
| 2026-11-14 | pass-001 run at `--max-new 0` (zero spend). Found and fixed a parser bug that had dropped `query_terms` from every theory with no counter-evidence |
| 2026-11-14 | Loop started: 30-minute cadence until 18:00 local, full-apply mode, output to this branch with a PR to `main` |
| 2026-11-14 | pass-002 (live, $0.10, 2 searches): found PMID:41564239, absent from the seeded queries. TH-007 `low` → `medium`; IDEA-003 split by format (`open` → `lit-checked`); `literature/REVIEW.md` §3 rewritten |
| 2026-11-14 | Rebased onto `origin/main` (3cbf4a7), which had moved 5 commits ahead. Noted that `main` lacks `LEO/.gitignore`, so the AMASS key is unprotected on a fresh clone — this branch carries the fix |
| 2026-11-14 | Push blocker resolved: `leossgh` has no write access to `niccolocard99/ai_hackaton`, so forked to `leossgh/ai_hackaton` and opened **PR #2** from the fork. Loop now pushes to `fork` each pass |
| 2026-11-14 | Repo scan of sibling branch `LEO-litreview-optimisation` (commits 9b45cdb, ae2cdef) found measured structural work. Recorded as **EXP-002**: TH-008 **retired**, TH-009 and TH-010 created at `high`, IDEA-004 `done`, IDEA-008/009/010 opened. First theory falsified by evidence rather than downgraded |
| 2026-11-14 | Cadence changed 30 min → 15 min, and a repo-scan step added to each pass so teammates' committed results feed the theories |
| 2026-11-14 | **PR #2 merged** into `main` at `27b6982` (46 files). Branch rebased onto the merge; subsequent passes go out as a new PR |
| 2026-11-14 | pass-003 ($0.15, 3 searches, scoped TH-009): found **nothing** on Y109 energetics — a negative result that escalates IDEA-010 to high. Side effects: TH-001 `medium` → `high` (PMID:29702641, 10-50 fold measured), and new **TH-011** on whether refold scores can separate close designs at all |
| 2026-11-14 | pass-004 ($0.15, 3 searches, scoped TH-011): repo scan found nothing new. TH-011 `medium` → **`high`** — AlphaFold3 needs **20 seeds** for a 2.9 Å median CDR-H3 (PMID:40814020), and four methods 2009-2025 converge on 2.4-2.9 Å. IDEA-006 rescoped from "second seed" to "measure the spread". New **TH-012**: a VHH is better supported by our tooling on two legs, so IDEA-005 → high |
| 2026-11-14 | **Repo scan found `origin/nick` `aec87f4`** (NICCOLO HER2 CDRH3 pre-filter). Recorded as **EXP-003**, the most consequential result so far: Boltz-2 ipTM separates real from shuffled HER2 pairs at AUROC 0.645 with 94.3% of decoys scoring ≥0.7, while a one-hot sequence model reaches Spearman +0.622 against measured log10 KD. New **TH-013/014/015**; TH-011's last caveat closed by direct measurement; EXP-001's scores shown to sit inside the decoy distribution. IDEA-011/012 opened, IDEA-006 downgraded to medium, IDEA-008 given measured substitutions |
| 2026-11-14 | pass-005 ($0.15, 3 searches, scoped TH-013): **no published validation of ipTM/ipAE as a binder screen exists** in what the searches returned, so our own decoy control is the best evidence on it. New IDEA-013 (ESM-IF1 vs ProteinMPNN, 6.5% vs 1.5% predicted binder success). Recorded the de novo ranking gap as open |
| 2026-11-14 | **Repo scan found three more teammate commits.** `3dc1f24` recorded as **EXP-004**: fold confidence calibrated against *measured KD* on 60 SPR variants — ipTM AUROC 0.502, min_ipae 0.497, sequence-only 0.849, and ipTM ≥ 0.7 passes **100%** of true non-binders. TH-013's "keep it as a coarse filter" clause **retracted**. `8743169`+`e542324` recorded as **EXP-005**: NanoMelt is the only usable Tm predictor, and **CDRH3 is thermostability-inert** (0 of 40 beyond MAE, framework control p=0.002). TH-005 → `high`, TH-007 **retired**, new TH-016/TH-017, IDEA-003 **done**, IDEA-014/015 opened |
| 2026-11-14 | Merged `origin/main` (track checkers moved `scripts/` → `NICCOLO/scripts/`, so base-validation commands in root `AGENTS.md` now point at the new path) |
| 2026-11-14 | pass-006 (3 searches, auto scope): TH-002 `medium` → **`low`** — EXP-005 shows CDRH3 cannot move Tm, so the trade-off it asserts does not exist in our design space; two papers show affinity and stability improved together. New IDEA-016 (expression/self-association scoring, a gap nothing covered). **Then AMASS credits ran out mid-pass**; `research_pass.py` fixed to stop on the first credits 403 rather than make seven doomed calls. pass-007 confirmed the stop works |
| 2026-11-14 | **PR #3 merged** into `main` at `da498b8`; PR #4 (nick) merged too, so the teammate work is on main. Merged `origin/main`, moved to fork branch `LEO-experiment-ledger-3` and opened a fresh PR |
| 2026-11-14 | **EXP-006, our own first experiment.** Hypothesis stated before the run: Protenix ranks binders better than Boltz-2. **Wrong** — all ten Protenix metrics span AUROC 0.426-0.526 on the same 60 measured variants, and `px_iptm_ab_her2`, the antibody↔HER2 chain-pair ipTM, is *below* chance at 0.448. TH-013 broadened from one model to the class. Closes a GPU-hungry line of enquiry for two minutes of analysis |
| 2026-11-14 | Repo scan also found `dac6935` (ranking power of fold metrics vs random) and `fd1b002` (`AGENT_GUIDE.md` for the pre-filter). `dac6935` committed the script but **not its output**, so its conclusion is recorded from the commit message, not from numbers we can cite. The guide supplied the out-of-fold figures now in TH-014: top 10 at 14 nM against 74 nM for a random 10 |
| 2026-11-14 | **PR #6 merged** at `661ddcc`. Merged `origin/main`, moved to fork branch `LEO-experiment-ledger-4` |
| 2026-11-14 | **EXP-007 corrects two of our own numbers.** The design band is **>=4** substitutions, not >=3 (3/10 = exactly 70%, and the rule is *strictly* under) — a teammate caught it independently in `fcb238e`. And the **1.42 nM top-20 median we published was never achievable**: that top 20 is full of 1-2 substitution variants the identity rule forbids. In the legal band it is **6.5 nM**. TH-014's numbers corrected to 0.913 / +0.628; new **TH-018** on the ~4x cost of the constraint |
| 2026-11-14 | **EXP-008 (ours) found the bar and the pattern.** Trastuzumab's own CDRH3 measures **0.152 nM**; exactly **1 of 21** legal variants beats it, by 1.6x. Across the 8 tightest legal variants M107 is mutated 8/8, D102 8/8, A106 7/8 — and A106/M107 are precisely the never-contacting positions from EXP-002. **IDEA-008's recipe revised** to D102 + A106 + M107 + one of F104/D108. New **TH-019**: binding alone is a weak basis for an S-tier claim |
