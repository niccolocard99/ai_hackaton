# Handoff for Jahn: what three agents learned, and what to port into the abpipe loop

**As of:** Saturday 14 Nov 2026, about 16:10 UTC · **Compiled from:** the handoff documents of the three Claude sessions that ran today (copies in [`sources/`](sources/)) · **Branch:** `LEO-handoff-for-jahn`, based on `origin/main` `9ac33d0`
**Written for:** Jahn, extending `JAHN/pipeline` (`abpipe iterate`). Read this page; open the sources only for the item you port.

**Deadline.** The user told one agent their deadline is 18:00 Oslo (17:00 UTC), and the VM shuts down at the end of 14 Nov. **Almost all of the code below is on the VM, not in git.** If you want a script, copy it off `~/leo-work/` or `~/wt/leo-design-loop/` now.

## The three engines, in one paragraph each

| Engine | What it does | Numbers (about 16:05 UTC) | Compute |
|---|---|---|---|
| **A. JapanFold RFD4 loop** (`sources/jf-loop/`) | Designs trastuzumab CDRH3 variants with RFD4-Proteina through the JapanFold REST API, gates them (ledger, abguard-equivalent check), keeps a ranked shortlist. Every arm except the `jfree*` controls pins Y105 | v3 runs 38 jobs per cycle (24 fresh + 14 refine), about 304 designs per 6.5 min, **$0.087/min, no spend cap**. Pool: 1,432 designs, 393 distinct, 262 eligible; merged shortlist 296 eligible | JapanFold API only; local CPU limited to about 2 cores |
| **B. design-loop** (`sources/design-loop/`) | Enumerates about 2.2 billion CDRH3 candidates near measured binders, plus RFD4 `her2-novel` on the H100, then gates (fast novelty index, your teammate's pre-filter, real abguard), adds NanoMelt Tm, assembles a hedged top-20 | Enumeration converged at round 5: 664 designs pass the real abguard. RFD4: two chunks of 40 generated, not yet harvested at 16:02 | H100 for RFD4 only; CPU capped by `run_limited.sh` |
| **C. de novo VHH** (`sources/denovo-vhh/`) | Designs new nanobodies on trastuzumab's epitope with RFD4 antibody weights; scores epitope geometry against 1N8Z | Batch 1: 99 of 104 pass the de novo track; 8 of 99 reach the team's epitope bar. Round 2 queued behind the H100 | H100, one RFD4 job at a time |

## What to port into your loop, in order of value

Your `LOOP.md` says thermostability, epitope and expression are not scored and that extra rounds mostly concentrate around the additive optimum. The agents' evidence agrees (engine B found the same optimum). So the gain is in **new kinds of module**, not more rounds.

| # | Port | Why | Where | Your module |
|---|---|---|---|---|
| 1 | **Design position S97 (and 98, 109), not only the 10-mer 99-108** | abguard's CDRH3 is the 13-mer 97-109 (`SRWGGDGFYAMDY`). The 525k reference never varies position 97, so any S97X change is +1 novelty distance from every known variant. S97 makes 0 heavy-atom contacts to HER2 in 1N8Z. An S97 substitution lifted 23 of 24 top parents to margin 4 or more while staying legal. **Effect of S97X on binding is unmeasured.** | A: FINDINGS section 5, HANDOFF "three things that matter most"; arms `jw97`, `jy97`, `jfree97`, `jyf97` | `scored_evolution` window, `window_novelty` |
| 2 | **Add a non-additive scorer: NanoMelt VH Tm** | Fills the thermostability gap in `LOOP.md`. Error 4.1 degC, trained on nanobodies (out of domain for a VH): tie-break only, never a gate. Needs `hmmscan` from `~/nick-work/nmenv/bin` on `PATH`; run CPU only | B: `gpu_stage.py`, `NICCOLO/thermo/` | new scorer next to `niccolo_prefilter` |
| 3 | **RFD4-Proteina as a generator**, not only mutagenesis | Produces CDRH3 loops far from the measured binders, which is where novelty margin is. Cost: the loops the classifier likes best fail novelty; the six that passed in a 12-design test had classifier p10 0.03-0.49. Novelty and classifier confidence pull in opposite directions | A: `jf_loop.py` (API, parallel). B: `rfd4_arm.py launch/harvest` (H100, one job at a time, yields 30 s to a queued waiter) | new generator |
| 4 | **Ledger gate before every generation and validation** | `check-run` before each run, `check --new-only` before validating, as `AGENTS.md` requires. Today three separate ledgers hold the same kind of rows (see below) | `LEO/tasks/generation-ledger/ledger.py` | a filter step before `abguard_hard` |
| 5 | **One fast validation path, not four** | Candidates: A's `fast_check.py` (rapidfuzz, identical to the team checker on 60/60 verdicts), B's `validate_abguard.py` (real checkers, rapidfuzz swapped in, 2 shards), `origin/main` #18 "fast exact novelty gate (6,300x)", and Niccolo's `~/nick-1m/pipeline/abpipe/abguard_fast.py` (not benchmarked by A). A Hamming index misses about 0.5% of indel neighbours, so **finalists still need the real abguard** | A: RUNBOOK section 6; B: `validate_abguard.py` | `abguard_hard`, `window_novelty` |
| 6 | **CPU etiquette** | The VM has 16 cores shared by the team. Unthrottled abguard spawns 86 threads. At 15:31 the load was about 70; a process with working directory `~/wt/jahn` was using about 380% with 86 threads and no thread limit. Use `OMP/MKL/OPENBLAS/POLARS/RAYON_NUM_THREADS=2`, `nice -n 10`, at most 2 validations at once | A: RUNBOOK; B: `run_limited.sh` | wrap your runner |
| 7 | **Keep Y fixed in some designs** | Standing user instruction. A read it as Y105 (largest buried surface, alanine-scan hotspot); if the user meant Y109, add a pin | A: HANDOFF "Instructions in force" | generator constraint |
| 8 | Epitope geometry check (de novo only) | Superpose HER2 first (RFD4 output does not keep the 1N8Z frame); add 488 to chain C and 590 to chain D to get 1N8Z numbering; needs `gemmi`, which the VM lacks | C: `epitope_check.py` | optional scorer |

## Findings that change the search space (measured, from the sources)

- A 13-mer CDRH3 needs **4 substitutions** to fall under 70% identity. The gate that actually binds is **novelty**: at least 3 mutations from every one of the 525k reference loops. Every measured SPR variant is itself in that reference, so **none can be submitted**. (B)
- abguard counts **indels**; a Hamming index misses about 0.5% of neighbours. (B)
- **All 20 of A's top designs sit at novelty distance exactly 3**: they pass at the threshold with no margin, and keeping Y105 and having margin are anti-associated (8% of 240 Y105-keeping eligible designs have margin, against 61% of 56 Y105-losing ones). The organisers' reference could differ from the local checker. (A)
- Held-out backtest in the legal band: the ranker's top pick beats trastuzumab's 0.152 nM **1.0%** of the time and lands within 3x **14%** of the time. A bootstrap lower bound and an anchored estimate did not beat the plain ranker. (B)
- **Boltz-2 does not separate binders**: ipTM of at least 0.7 for 94% of known non-binders (EXP-011). The team banned it. (C)
- `freefw` (18 framework positions freed) was not worth scaling: the model changed 3 of 18. (C)

## Instructions in force (from the user)

No Boltz anywhere. Generate with RFD4-Proteina. No ProteinMPNN (retired). One RFD4 job on the H100 at a time. About 3 CPU cores per agent. Ledger check before every generation run. JapanFold API with no spend cap. Some designs keep the Y fixed. Do not commit, push or merge into `main` unless asked.

## Where the code is

| Engine | VM | In git |
|---|---|---|
| A | `~/leo-work/` (`h100/jf_loop.py`, `fast_check.py`, `finalize.py`, `compare_lists.py`; log `jf/loop.log`; shortlist `final/candidates.csv`); copies of the docs in `~/leo-work/handoff/` | `LEO/tasks/rfd4-japanfold/` on branch `LEO-litreview-optimisation`; everything after commit `1f49287` is uncommitted |
| B | `~/wt/leo-design-loop/LEO/tasks/design-loop/` (branch `LEO/design-loop`, nothing committed) | laptop copy only, uncommitted |
| C | `~/hackathon-team-3/LEO/tasks/denovo-vhh-epitope/` | branch `LEO-denovo-vhh-epitope`, 18 files uncommitted |

**Ledgers.** Three copies exist: B's `ledger/` (about 1.4k variants), A's `~/leo-work/generation-ledger/` (about 1.6k), and the repo copy (nearly empty). The CSVs are `merge=union` and append-only. Merging them is Leo's call.

## Decisions that belong to Leo, not to you

Which four designs lead the list and in what order; whether RFD4 loops take lead slots (their classifier p10 is low); VL D28K (unmeasured, and D28 contacts HER2); whether to apply the S97 hedge to the lead designs; ledger merge; when to commit; whether the de novo VHH queue may use `FORCE=1`.

## Caveats on this brief

- Compiled from the four documents, read at about 16:10 UTC. The numbers move every few minutes. I checked only live facts (processes, load, ledger sizes), **not the science**: the claims above are each agent's, and each cites its own evidence file.
- Local checker verdicts are predictions. The real track assignment is abguard's at submission time, and the installed abguard is not confirmed as the official release.
- Two source documents say the "no spend cap" instruction came from a peer session claiming to relay the user. **The user confirmed it directly.**
- Engine A's comparison of its list with B's (FINDINGS section 6) predates its S97 arms; re-run `compare_lists.py` before quoting it.

## Do not

Start Boltz or ProteinMPNN. Run unthrottled validation. Use `pkill -f` or `pgrep -f` with a pattern that also appears in your own command line (it killed two agents' own shells today); use PIDs. Edit or delete ledger rows by hand. Commit keys, the 130 MB `ref/` table or organiser `data/`. Switch branches in `~/hackathon-team-3` on the VM: Niccolo works there directly.
