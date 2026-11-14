# generation-ledger

**Branch:** `LEO-litreview-optimisation` (uncommitted; not yet on its own branch, see Next) · **Started:** 2026-11-14 · **Agent:** Claude Sonnet 5.5

**Read both agents files:** [`../../../AGENTS.md`](../../../AGENTS.md) (root, shared) ☑ · [`../../AGENTS.md`](../../AGENTS.md) (Leo's) ☑

## What this task is

A CSV ledger of every generation run and every sequence it produced, plus a checker, so nobody regenerates, refolds or re-validates a result the team already has. The GPU is shared and only 4 designs get tested, so duplicate runs cost real time. The rule that makes agents use it is the **GENERATION LEDGER** section of the root `AGENTS.md`.

This is not the experiment ledger on `LEO-experiment-ledger` (`experiments.csv`: one row per research experiment, with theories and verdicts). This one is per run and per sequence, and its only job is deduplication.

## Right now

idle

## Done so far

| What | Verified how | Commit |
|---|---|---|
| `ledger.py`: `check-run`, `add-run`, `set-status`, `check`, `add`; stdlib only | 18 unit tests in `test_ledger.py` pass; two injected bugs (stop following the spec's `input:` file, wrong blocking statuses) were each caught by the matching test | uncommitted |
| Run fingerprint covers method, model, input file contents (including the structure a `.yaml` spec points at), seed, settings | real CLI on `apheris/rfd4/Data/12-her2/her2_cdrh3.yaml`: both files hashed, hashes match `sha256sum` | uncommitted |
| `runs.csv` backfilled with the HER2 smoke test (`unknown-20261114-her2-smoke`) | `check-run` with the smoke-test flags now exits 3 | uncommitted |
| Root `AGENTS.md`: GENERATION LEDGER section plus a Rules-for-agents bullet | read back; anchored on lines `main` has not changed | uncommitted |

## Blocked

| Item | Blocker | Needed to unblock |
|---|---|---|
| Backfill of the smoke test's sequences | `runs/her2/design` exists only on the VM (git-ignored) | Run `ledger.py add` on the VM against its designs table with `--run-id unknown-20261114-her2-smoke` |
| Verifying the smoke-test input hashes | They were taken from this checkout's `Data/`, not the VM's | `sha256sum` the two files on the VM and compare with `input_sha256` in `runs.csv` (`333d6213f4e8`, `2920ebcf369a`) |

## Next, in order

1. Decide where the ledger lives. It sits under `LEO/` to follow the folder-scope rule, but the root rule applies to the whole team, so a teammate's agent has to write into Leo's folder. Moving it means changing one path in `AGENTS.md` and `ledger.py` needs no change (`--ledger-dir`).
2. Give it its own branch and commit it (the work is uncommitted; no commit was asked for).
3. On the VM, backfill any other runs already in `apheris/rfd4/runs/`.

## Log

| Date | Event |
|---|---|
| 2026-11-14 | Task started. Existing `LEO-experiment-ledger` work reviewed and found to be a different thing (research experiments, not runs), so this is separate |
| 2026-11-14 | `ledger.py`, tests, header-only `variants.csv`, `.gitattributes` (`*.csv merge=union`) written; 18 tests pass |
| 2026-11-14 | HER2 smoke test backfilled as a run row; `check-run` on it exits 3 |
| 2026-11-14 | Caught up after `origin/main` was merged in (f037ab9). My uncommitted section and bullet survived in `AGENTS.md`, tests still pass, CSVs untouched. Re-read `NICCOLO/PIPELINE.md` and `NICCOLO/eda/AGENT_GUIDE.md`. Added "log the shortlist, not the screen" because the pre-filter scores millions of CDRH3s |
| 2026-11-14 | Not wired in, in other people's folders (not touched): `ANDERS/vhh_variants.py` (a generator) and `NICCOLO/eda/prefilter.py`. Their agents have to call `ledger.py` by hand, as the `AGENTS.md` rule says |
