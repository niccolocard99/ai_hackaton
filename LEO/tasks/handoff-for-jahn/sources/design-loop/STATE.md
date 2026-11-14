# design-loop

**Branch:** `LEO/design-loop` (VM worktree `~/wt/leo-design-loop`, off `origin/main` 5d639b1; nothing committed) · **Started:** 2026-11-14 · **Agent:** Claude Opus 5, then Claude Sonnet 5.5 (model switched mid-session)

**Read both agents files:** `../../../AGENTS.md` (root, shared) ☑ · `../../AGENTS.md` (Leo's) ☑. The root file's GENERATION LEDGER section is uncommitted on the laptop, not on the VM; it was read there.

## What this task is

A pipeline on the VM that generates optimisation-track CDRH3 designs (enumeration near measured binders, and RFD4-Proteina her2-novel), gates them (identity, novelty, the team classifier, the REAL abguard), predicts VH Tm with NanoMelt, and assembles a hedged ranked list, since only the top 4 designs are guaranteed a test. The handoff document (Claude Docs, "Handoff: HER2 optimisation-track design loop") has the full picture. See `README.md`.

## Right now

16:02 UTC. RFD4 chunk `novel_s302` has generated its 40 samples and is finishing; my chain (`out/rfd4_chain.sh`) ends after it and nothing else of mine is running. hackathon-3e's waiter (`wait_and_run.sh`) is queued for the card behind it, so launch nothing further until their job has run. Next is the harvest of both chunks, then `final_portfolio.py`.

## Done so far

Verified by running it; nothing is committed.

| What | Verified how |
|---|---|
| abguard rules measured: 13-mer CDRH3 needs 4 substitutions; novelty (3+ from all 525k reference loops) is the binding gate; no measured SPR variant is submittable | `probe_abguard.py`, `feasible.py` |
| Fast novelty index vs real abguard | 204/204, 344/346, 527/529; misses are indel neighbours, so the real abguard screens everything |
| Enumeration: 5 rounds, about 2.2 billion candidates, 668 shortlisted, 664 pass real abguard; converged (round 2 added 34 designs) | `out/round_*/funnel.json`, `out/pool.csv` |
| Ledger protocol on every generation run; a restart reused 529 verdicts instead of redoing them | `ledger/runs.csv`, `out/loop.log` |
| Held-out backtest: lower bound and anchored estimate are neutral; legal-band top pick beats 0.152 nM 1.0%, within 3x 14% | `backtest.py` |
| Framework Tm arm (R87S, K76L) and literature arm (D102F/W loops, VL D28K/R); D28 faces Glu598 at 4.5 A in 1N8Z | `tm_arm.py`, `adapt_arm.py` |
| RFD4 her2-novel: `novel_s301` (40) generated and ledgered; `novel_s302` (40) generated. 12-design harvest test: 6 of 12 loops fail real abguard; the 6 that pass have classifier p10 0.03 to 0.49 | `rfd4_arm.py`, `rfd4_harvest.py` |
| CPU rule: `run_limited.sh`, abguard capped at 2 shards. Fairness: `rfd4_arm.py` yields 30 s to a queued `wait_and_run.sh` waiter and waits for its job | code, peer requests |

## Blocked

| Item | Blocker | Needed to unblock |
|---|---|---|
| Final list v3 | `final_portfolio.py` is written and syntax-checked but has never run; it needs the RFD4 harvests | Harvest both chunks, then run it |
| Canonical ledger | `LEO/tasks/generation-ledger/` is uncommitted on the laptop; the VM cannot see it | Merge `ledger/runs.csv` and `variants.csv` (they are `merge=union`) |
| Ledger smoke-test row | Its input hashes came from the laptop and differ from the VM's | Re-hash on the VM, probably CRLF vs LF |
| Developability scorers | Aggrescan3D, SoluProt, NetSolP do not install; ProteinMPNN and Boltz are banned | None chosen |

## Next, in order

1. Harvest `novel_s301`, then `novel_s302`: `./run_limited.sh ~/.venvs/abpipe/bin/python rfd4_arm.py harvest --tag <tag>`.
2. Run `./run_limited.sh ~/.venvs/abpipe/bin/python final_portfolio.py`; check every row PASSes the optimisation track and FAILs de novo.
3. Settle the lead four with the team (the decisions are listed in the handoff document).
4. Copy `out/` and `ledger/` to the laptop before the VM shuts down at the end of today; never copy `ref/` or `data/`.

## Log

| Date (UTC) | Event |
|---|---|
| 2026-11-14 14:05 | Pilot generator run, discarded: its epitope-risk column used guessed contact counts |
| 2026-11-14 14:38 | Loop launched; first attempt blocked waiting for the GPU. GPU work moved to a separate worker, restarted |
| 2026-11-14 14:49 | Loop converged at round 5; backtest and add-on arms followed |
| 2026-11-14 15:18 | 20 exact designs folded with Boltz-2 (later banned; results unused) |
| 2026-11-14 15:27 | ProteinMPNN run launched; later found INVALID (a comment swallowed the soluble flag) and killed by a peer. Retired |
| 2026-11-14 15:31 | Team: no Boltz anywhere. Worker and my one Boltz container stopped; generation moved to RFD4 |
| 2026-11-14 15:34 | RFD4 `novel_s301` launched (ledger run `leo-rfd4-novel_s301`); `novel_s302` followed at 15:48 |
| 2026-11-14 15:48 | hackathon-3e's waiter was starved by the immediate chunk 2 launch; fairness yield added at 16:00 |
| 2026-11-14 16:02 | Handoff document written; README and STATE brought up to date |
