# design-loop

**Branch:** `LEO/design-loop` (VM worktree `~/wt/leo-design-loop`, off `origin/main` 5d639b1; nothing committed) · **Started:** 2026-11-14 · **Agent:** Claude Opus 5, then Claude Sonnet 5.5 (model switched mid-session)

**Read both agents files:** `../../../AGENTS.md` (root, shared) ☑ · `../../AGENTS.md` (Leo's) ☑. The root file's GENERATION LEDGER section is uncommitted on the laptop and is not on the VM; it was read there.

## What this task is

A loop on the VM that enumerates optimisation-track CDRH3 designs, gates them, scores them with the team pre-filter (with bootstrap uncertainty), validates them with the REAL abguard checkers, and folds and scores survivors on the shared H100. It adds what the abpipe loop lacks (epitope contact, Tm, structure, literature levers) and builds a hedged ranked list, since only the top 4 are guaranteed a test. See `README.md`.

## Right now

- The CPU loop **converged after round 5** and sits in its `finalize` step, waiting up to 20 min for the GPU worker (PID in `out/loop.pid`).
- The GPU worker (`out/gpu_worker.pid`) is waiting for the H100: a teammate's RFD4 job (`antibody-2026-10-30`) has held it for 12+ minutes. The team rule is one heavy job at a time, so it does not start beside it. 20 exact final sequences are queued in `out/fold_queue.csv`.
- The ranked list is built: `out/FINAL.md`, `out/FINAL.csv`.

## Done so far

Verified by running it; nothing is committed.

| What | Verified how |
|---|---|
| abguard rules measured, not assumed: 13-mer CDRH3 needs >= 4 substitutions; novelty (>= 3 from all 525k reference loops) is the binding gate; no measured SPR variant is submittable | `probe_abguard.py`, `feasible.py`; real checker on a 7-rung ladder |
| Fast novelty index vs real abguard | 204/204, 344/346, 527/529; the misses are indel neighbours (abguard counts Levenshtein), so real abguard screens every shortlist |
| 5 rounds: 2.3 billion candidates enumerated, 668 shortlisted, 664 PASS real abguard; round 2 added only 34 new designs | `out/round_*/funnel.json`, `out/pool.csv` |
| GENERATION LEDGER protocol every round (check-run, add-run, set-status, check, add); restart reused 529 verdicts and 527 Tm predictions instead of redoing them | `ledger/runs.csv`, `out/loop.log` |
| Held-out backtest: bootstrap lower bound and anchored estimate do NOT beat the plain ranker; legal-band top pick beats trastuzumab's KD 1.0%, lands within 3x 14% | `backtest.py`, 20 repeats of cluster-grouped 5-fold CV |
| Framework Tm arm (R87S, K76L, structure-filtered): 12 variants, all PASS | `tm_arm.py`, `out/tm_arm.csv` |
| Literature arm (D102F/W CDRH3, VL D28K/R from the ADAPT papers): 16-18 variants, all PASS | `adapt_arm.py`, `out/adapt_arm.csv` |
| Final list: 20 designs, 20/20 PASS optimisation track, 0/20 de novo (expected); every sequence diffed residue by residue against trastuzumab, 0 mismatches | `final_portfolio.py`, `out/FINAL.csv` |

## Blocked

| Item | Blocker | Needed to unblock |
|---|---|---|
| Boltz-2 folds of the final 20 | The H100 is held by a teammate's job; team rule is one heavy job at a time | A gap, or Leo's go-ahead to fold alongside (about 25 structures, 4 min, ~12 GB beside their 3.6 GB) |
| Canonical ledger | `LEO/tasks/generation-ledger/` is uncommitted on the laptop; the VM cannot see it | Merge `ledger/runs.csv` and `ledger/variants.csv` from here (the CSVs are `merge=union`) |
| Ledger smoke-test row | Its input hashes came from the laptop and differ from the VM (`668d383aec56` / `8136a9abb040` here) | Re-hash on the VM, probably CRLF vs LF |
| Developability scorers | Not installed | None chosen |

## Next, in order

1. Fold the final 20 as soon as the card is free (worker does it); rerun `final_portfolio.py` to merge them.
2. Decide the lead four with the team; copy `out/` and `ledger/` off the VM before it shuts down.
3. Fold or drop designs whose H3 is flagged as displaced from the control.

## Log

| Date (UTC) | Event |
|---|---|
| 2026-11-14 14:05 | Pilot generator run, then discarded: its epitope-risk column used guessed contact counts, not the 8-structure table |
| 2026-11-14 14:17 | Real abguard on 204 designs: all PASS, assigned `novel` |
| 2026-11-14 14:33 | First end-to-end test on scratch folders: NanoMelt failed, `hmmscan` not on PATH. Fixed |
| 2026-11-14 14:38 | Real loop launched; first attempt blocked inside round 1 waiting for the GPU. Killed, GPU work moved to a separate worker, restarted |
| 2026-11-14 14:47 | Round 2: 1.65 billion enumerated, 34 new designs. The additive surrogate has one optimum, as the abpipe doc says |
| 2026-11-14 14:49 | Loop converged at round 5 |
| 2026-11-14 14:50 | Backtest: lower bound and anchored estimate are neutral, not better. Legal-band hit rate 1% |
| 2026-11-14 14:50 | Tm arm. Only two single mutations clear +1 C once ->C/P/E/D are excluded; both inside NanoMelt's error |
| 2026-11-14 14:54 | Literature arm from `findings/08-hard-numbers-verified.md` |
| 2026-11-14 14:59 | Final list built; repeated `final` runs in the ledger are annotated as superseded by `leo-dl-final-1459` |
