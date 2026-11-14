# denovo-vhh-epitope

**Branch:** `LEO-denovo-vhh-epitope` (the template's `LEO/<name>` form is impossible: a branch named `LEO` already exists, so git cannot create `LEO/...`) · **Started:** 2026-11-14 · **Agent:** Claude Sonnet 5.5

<<<<<<< HEAD
**Read both agents files:** `../AGENTS.md` (root, shared) ☑ · `LEO/AGENTS.md` (Leo's) ☑ · `LEO/README.md` ☑ · after the 14 Nov updates I read the changed parts of the root `AGENTS.md` (base validation, GENERATION LEDGER, rules), `NICCOLO/eda/AGENT_GUIDE.md` to section 6, `NICCOLO/PIPELINE.md` to section 3, `ledger.py`'s CLI, and, from `origin/nick` (read only), `EXP-011.md` and `NICCOLO/SCORED_DENOVO_VHH.md`. Not read in full.
=======
**Read both agents files:** `../AGENTS.md` (root, shared) ☑ · `LEO/AGENTS.md` (Leo's) ☑ · `LEO/README.md` ☑ · after the 14 Nov update I read the changed parts of the root `AGENTS.md` (base validation, GENERATION LEDGER, rules, from the main checkout where the ledger section is uncommitted), `NICCOLO/eda/AGENT_GUIDE.md` to section 6, `NICCOLO/PIPELINE.md` to section 3, and `ledger.py`'s CLI. Not read in full.
>>>>>>> origin/main

## What this task is

Generate completely novel VHH binders to HER2 domain IV that sit on trastuzumab's epitope, targeting the de novo track. The epitope constraint is what the competitive-binding assay scores. Reasoning and assumptions are in [`NOTES.md`](NOTES.md).

## Right now

<<<<<<< HEAD
Round 2 is queued on the VM (`team-2`, repo `~/hackathon-team-3`), detached, strictly one GPU job at a time, waiting for another agent's `her2-novel` run (`optimistic_wu` at the last check). RFD4-Proteina only, per Leo. `queue_round2.sh` runs: (1) RFD4 partial-diffusion refinement of the 5 strong-epitope parents, 16 children each; (2) a fresh batch with CDR3 17 to 21, 100 designs, seed 3. No refold step. Log: `runs/round2.log` on the VM. Nothing from round 2 has been generated yet.

## Done so far

Commits on this branch: `68f0ce3` (initial setup) and the merge `215f1ed` of `origin/main` (clean, 70 commits). Everything else listed here is uncommitted. The folder is copied to the VM at `~/hackathon-team-3/LEO/tasks/denovo-vhh-epitope`.

| What | Verified how |
|---|---|
| Structure: 1N8Z domain IV crop; 5MY6 rejected (its nanobody binds domain I, 43.7 Å from trastuzumab) | `compare_5my6.py`, `compare_5my6.txt` |
| Hotspots: the organisers' four atoms, all within 4.5 Å of trastuzumab CDRH3, max 12.4 Å apart. **Unchanged in round 2** | `check_spec.py` |
| Smoke test (2 per variant) and batch 1 (`fixedfw`, 100 designs, seed 2): spec accepted, all conditions resolved as designed, ledger-registered | runs `vhh-fixedfw-n2-s1`, `vhh-freefw-n2-s1`, `vhh-fixedfw-n100-s2` |
| Base validation on 104 designs: 99 pass the de novo track; the 5 failures have 12-residue CDR3s and were assigned to optimisation | `validate_and_log.sh` on the VM; all 104 logged to the ledger |
| Epitope check on all 104 (HER2 superposed onto 1N8Z): 8 of 99 valid designs reach the team's bar (cov >= 0.6, in_ep >= 0.7) | `epitope_check.py`, `results/batch1-epitope.csv` |
| Epitope-first ranking reproduces Nick's four best (`fixedfw_48`, `_89`, `freefw_0`, `_42`) | `score_designs.py`, `results/batch1-scored.csv` |
| Boltz-2 refold of 12 (kept; no more will run): 2 of 12 reproduce the designed pose | `refold_each.py`, `results/batch1-shortlist-refold-boltz2.csv` |
| `freefw` gains no identity margin (83% recovery of freed positions) | `NOTES.md` section 5 |
| NanoMelt on all 104 finished on the VM (`reports/batch1_nanomelt.csv`); not analysed here: tie-break only, no thermostability claim | see EXP-011 |
| Round-2 tooling: `make_refine_spec.py`, refine mode in `run_design.sh`, `wait_and_run.sh`, `queue_round2.sh`, retry in `refold_each.py`, CPU caps (`--cpus 3`, threads 2) | `bash -n`, syntax checks; not yet run |
=======
Generating. `./run_design.sh fixedfw 100 2` is running on the VM's H100 (`team-2`, repo `~/hackathon-team-3`, VM branch `main-loop`), ledger run `vhh-fixedfw-n100-s2`, container `dreamy_bassi`. 87 of 100 designs were written at the last check. When it ends, the scoring pipeline runs on all designs (section Next).

## Done so far

First commit `68f0ce3` on this branch (initial setup); everything marked "uncommitted" below is in the working tree only. The whole folder is copied to the VM at `~/hackathon-team-3/LEO/tasks/denovo-vhh-epitope`.

| What | Verified how |
|---|---|
| Structure: 1N8Z domain IV crop. 5MY6 rejected: its nanobody (2Rs15d) binds domain I, 0 shared epitope residues with trastuzumab, 43.7 Å apart | `compare_5my6.py`, output in `compare_5my6.txt` |
| Hotspots = the organisers' four atoms (D593 NZ, F573 CD1, I591 CG2, P572 CD), all within 4.5 Å of trastuzumab CDRH3, max separation 12.4 Å | `check_spec.py` passes; deliberately broken specs fail with exit 1 |
| Smoke test (N=2 per variant, seed 1): the spec is accepted and every condition resolves as designed (192 sequence residues, 436 coordinate atoms, 4 hotspots, 332 internal-coordinate atoms); `hotspot_obedience` 1.0, no clashes | runs `vhh-fixedfw-n2-s1`, `vhh-freefw-n2-s1`, ledger-registered, outputs under `runs/` on the VM |
| The 4 smoke designs pass `NICCOLO/scripts/check_denovo.py` as `de_novo` (identity to trastuzumab 65–68%) | run on the VM with the installed `abguard` |
| `freefw` redesigns only 3 of its 18 freed positions (83% recovery) and gains no identity margin, so no 100-design `freefw` run was queued | `NOTES.md` section 5, item 2 |
| HER2 does not keep the 1N8Z frame in the output (about 150 Å raw) but is rigid (0.05 Å after superposition) | `epitope_check.py` |
| Early read of 54 finished `fixedfw` designs: 29 have 4 of 4 hotspots in contact, 16 have 3; every design overlaps trastuzumab's Fv | `epitope_check.py` (uncommitted) |
| Tools, tested on the 4 smoke designs only: `make_designs_table.py`, `epitope_check.py` (needs `gemmi`, which the VM lacks, so it runs on a laptop), `score_designs.py` | uncommitted |
>>>>>>> origin/main

## Blocked

| Item | Blocker | Needed to unblock |
|---|---|---|
<<<<<<< HEAD
| Round 2 launching | The H100 is held by `heuristic_curran`; the queue waits (rule: one RFD4 job at a time) | That run finishing |
| Committing the new tools, results and notes | Leo asked for the first commit and for the `origin/main` merge; nobody has asked for more commits | Leo's go-ahead |
| EXP-010 / EXP-011 and IDEAs 018 to 020 | They are on `origin/nick` (PR #17), not on main; I read them without merging | Leo says whether to merge that branch, or PR #17 lands on main |
| Official track assignment | The installed `abguard` is not confirmed as the organisers' release; the local checker measures identity to trastuzumab only, not to caplacizumab (the scaffold) | The official `abguard` run |
| Untested: multi-parent refine spec | First use of `input_as_annotated_cif` with a partial-trajectory fraction here | Step 1 of round 2; fallback is one parent per run |

## Next, in order

1. Watch `runs/round2.log`. When step 1 ends: `validate_and_log.sh`, fetch the CIFs, `epitope_check.py`, `score_designs.py`. Compare the share reaching the strong-epitope bar with batch 1's 8%.
2. Same for step 2 (the fresh batch).
3. Final shortlist and the reasoning per design. Leo decides the ranking order for the 4 guaranteed slots.
4. Leo decides on commits and on merging `origin/nick`.
=======
| Committing the new tools and notes | Leo asked for the first commit explicitly; nobody has asked for more | Leo's go-ahead |
| Catching up with `origin/main` (a peer session asked, merge `f037ab9`) | A merge creates a merge commit, and nobody has asked for one | Leo's go-ahead |
| Official track assignment | The installed `abguard` is not confirmed as the organisers' release, and the local checker measures identity to trastuzumab only, not to caplacizumab (the scaffold) | The official `abguard` run |

Resolved: the ledger lives on the VM at `/home/hackathon/leo-work/generation-ledger/ledger.py` (run with `LEDGER=` set to it). The GPU was shared with Leo's design loop; Leo said to use the H100, so the runs went alongside it with `FORCE=1`.

## Next, in order

1. When the batch ends: `make_designs_table.py`, ledger `check --new-only`, `check_denovo.py` on the VM, ledger `add` (with `base_validation`, `track_assigned`, `cdrh3`).
2. Fetch the CIFs, `epitope_check.py` on all designs, `score_designs.py`, read the shortlist.
3. Refold the shortlist with Boltz-2 (`refold.sh`) as a pose sanity check and for provenance, not as a screen. NanoMelt for thermostability.
4. Record the reasoning per shortlisted design. Decide whether another seed or a `freefw` batch is worth the GPU time.
5. Leo decides on commit and the `origin/main` catch-up.
>>>>>>> origin/main

## Log

| Date | Event |
|---|---|
<<<<<<< HEAD
| 2026-11-14 | Branch and worktree created off `LEO`; epitope measured (32.4 Å, too wide for one VHH); identity edge found; `freefw` added |
| 2026-11-14 | 5MY6 rejected (domain I epitope); hotspots swapped to the organisers' four CDRH3-nearest atoms; compactness bound relaxed 12 to 20 Å |
| 2026-11-14 | Committed `68f0ce3`, copied to the VM, smoke runs and batch 1 launched alongside the design loop (`FORCE=1`, at Leo's "use h100") |
| 2026-11-14 | Wrong assumptions corrected: output does not keep the 1N8Z frame; `freefw` recovery is 83%, not 30% |
| 2026-11-14 | Bug fixed in `score_designs.py`: free-gap CDRH3 identity merged unrelated loops; now edit distance |
| 2026-11-14 | Boltz-2 refold of the 12 ran before the team said no Boltz; stopped nothing (it had finished); Boltz option removed from `refold.sh`; results kept |
| 2026-11-14 | Merged `origin/main` at Leo's request (`215f1ed`, clean). Read EXP-011 from `origin/nick`: epitope coverage is the only usable axis; Protenix, not Boltz |
| 2026-11-14 | Scoring re-ranked epitope-first. Round 2 designed from that: refine the 5 strong parents, fresh CDR3 17 to 21 batch, Protenix refold. Queued behind `heuristic_curran` |
| 2026-11-14 | Slip: `pkill -f` matched my own ssh command line and killed my earlier waiter and that ssh session. Nothing else was affected; the round-2 queue was launched afterwards and checked |
| 2026-11-14 | Leo: RFD4-Proteina ONLY. Dropped the Protenix v2 refold step; relaunched the queue (steps 1 and 2 only) by PID, before anything of mine had started. Still waiting for the GPU |
=======
| 2026-11-14 | Branch and worktree created off `LEO` (`.worktrees/denovo-vhh`) |
| 2026-11-14 | Measured the epitope: too wide (32.4 Å) for one VHH; hotspots must sit in one patch |
| 2026-11-14 | Found the identity edge: fixed caplacizumab framework gives a 65–71% floor, on the 70% line. Added the `freefw` variant (18 exposed positions redesigned) |
| 2026-11-14 | Dead end: 5MY6 (HER2 + nanobody 2Rs15d) was proposed as the source of hotspots. Its epitope is HER2 87–195 (domain I), 85.8 Å from trastuzumab's; unusable for a VHH that must compete with trastuzumab |
| 2026-11-14 | Swapped hotspots to the organisers' four CDRH3-nearest atoms (my first pick, E558/D560/P572/F573, is dropped). Relaxed the compactness check from 12 to 20 Å (the paratope span) because the set measures 12.4 Å |
| 2026-11-14 | `run_design.sh` follows the generation ledger; read-only SSH to `team-2` confirmed it is our VM |
| 2026-11-14 | Committed `68f0ce3` and copied the folder to the VM at Leo's request |
| 2026-11-14 | Smoke runs `vhh-fixedfw-n2-s1` and `vhh-freefw-n2-s1` done on the H100 (run alongside the design loop, `FORCE=1`, at Leo's "use h100"). The first launch hung the ssh call (a `&` wrapped the whole `cd && …` chain); the job itself was unaffected |
| 2026-11-14 | Launched `vhh-fixedfw-n100-s2` (100 designs, seed 2), queued behind the smoke runs |
| 2026-11-14 | Wrong assumption corrected: the output does not keep the 1N8Z frame (about 150 Å raw RMSD). `epitope_check.py` superposes first. Wrong assumption corrected: `freefw` recovery is 83%, not 30% |
| 2026-11-14 | Bug found and fixed in `score_designs.py`: free-gap CDRH3 identity merged unrelated loops (about 59% identical by eye); now edit-distance identity |
>>>>>>> origin/main
