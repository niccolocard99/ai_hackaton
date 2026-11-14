# handoff-for-jahn

**Branch:** `LEO-handoff-for-jahn` (from `origin/main` `9ac33d0`) · **Started:** 2026-11-14 · **Agent:** Claude Sonnet 5.5

## What this task is

A brief for Jahn, summarising the handoff documents of three agents (JapanFold RFD4 loop, design-loop, de novo VHH) and listing what to port into the `abpipe` loop. It carries copies of those documents in `sources/` and no code: the code is uncommitted or on the VM only. Matters for the pitch because it records the reasoning behind the loop and keeps the pipeline reproducible.

## Right now

Branch pushed to the fork; pull request page opened for Leo to create.

## Done so far

- Read the three handoffs: `LEO/tasks/rfd4-japanfold/HANDOFF.md` (plus RUNBOOK and FINDINGS), `denovo-vhh-epitope/HANDOFF.md`, and the design-loop handoff, which was a Claude Doc and is exported to `sources/design-loop/HANDOFF.md`.
- Wrote `README.md`: engines, eight ports ranked by value, findings, instructions in force, code locations, caveats.
- Checked live facts on the VM (load, containers, ledger sizes). Did not re-verify any scientific claim.

## Blocked

- The code behind the ports is not in git. Committing it is Leo's call (`AGENTS.md`: do not commit unless asked).
- Pull request creation: the GitHub CLI is installed but not logged in, so the page is opened in the browser for Leo to confirm.

## Next

1. Leo creates the pull request.
2. Jahn picks ports 1 to 3 first; port 1 (design S97) needs no new infrastructure.
3. Re-sync `sources/` if an agent updates its handoff before the VM shuts down.

## Log

- 2026-11-14 16:10 UTC: brief compiled from handoffs written 16:00-16:05 UTC.
