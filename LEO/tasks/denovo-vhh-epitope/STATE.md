# denovo-vhh-epitope

**Branch:** `LEO-denovo-vhh-epitope` (the template's `LEO/<name>` form is impossible: a branch named `LEO` already exists, so git cannot create `LEO/...`) · **Started:** 2026-11-14 · **Agent:** Claude Sonnet 5.5

**Read both agents files:** `../AGENTS.md` (root, shared) ☑ · `LEO/AGENTS.md` (Leo's) ☑ · `LEO/README.md` ☑ · after the 14 Nov update I read the changed parts of the root `AGENTS.md` (base validation, GENERATION LEDGER, rules, from the main checkout where the ledger section is uncommitted), `NICCOLO/eda/AGENT_GUIDE.md` to section 6, `NICCOLO/PIPELINE.md` to section 3, and `ledger.py`'s CLI. Not read in full.

## What this task is

Generate completely novel VHH binders to HER2 domain IV that sit on trastuzumab's epitope, targeting the de novo track. The epitope constraint is what the competitive-binding assay scores. Reasoning and assumptions are in [`NOTES.md`](NOTES.md).

## Right now

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

## Blocked

| Item | Blocker | Needed to unblock |
|---|---|---|
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

## Log

| Date | Event |
|---|---|
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
