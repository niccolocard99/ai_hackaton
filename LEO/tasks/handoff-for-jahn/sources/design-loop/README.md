# design-loop: optimisation-track generate, gate and score loop

Runs on the team VM; everything here was executed there. The handoff document (Claude Docs, "Handoff: HER2
optimisation-track design loop") has the current status and the decisions that are the user's.

```
generate: enumeration near measured binders (run_loop.py)  or  RFD4-Proteina her2-novel (rfd4_arm.py)
  -> cheap gates: >= 4 substitutions, no cysteine, novelty (fast index)
  -> NICCOLO pre-filter (classifier + ranker) with a 100-fit cluster bootstrap
  -> generation ledger check, then the REAL abguard track checkers (BASE VALIDATION)
  -> NanoMelt Tm of the VH
  -> add-on arms (tm_arm.py, adapt_arm.py) and final_portfolio.py -> ranked list
```

## Run it (VM, from this folder; run_limited.sh = nice 10 plus a 2-thread cap, the team's CPU rule)

```bash
export LEDGER_OWNER="<who you are>"            # else the shared VM git identity is recorded as owner
mkdir -p data ref
cp ~/hackathon-team-3/apheris/hackathon-data/her2-trastuzumab/trastuzumab_her2_spr_all_768.csv data/
unzip ~/hackathon-team-3/apheris/hackathon-data/her2-trastuzumab/tz_her2_facs_3classes_524346.csv.zip -d ref/
cp <generation-ledger>/{ledger.py,runs.csv,variants.csv} ledger/    # stdlib-only tool; never copy data/ or ref/ into git
PY=~/.venvs/abpipe/bin/python                  # abguard, rapidfuzz, sklearn, pandas, polars
./run_limited.sh $PY run_loop.py --until 15:50 --max-rounds 8 --no-boltz     # enumeration (it converges by round 5)
./run_limited.sh $PY rfd4_arm.py launch  --tag novel_s303 --n 40 --seed 303  # RFD4 on the H100, ledger-checked
./run_limited.sh $PY rfd4_arm.py harvest --tag novel_s303                    # extract, abguard, classifier, NanoMelt
./run_limited.sh $PY tm_arm.py ; ./run_limited.sh $PY adapt_arm.py           # add-on arms (need out/portfolio.csv)
./run_limited.sh $PY final_portfolio.py                                      # out/FINAL.md, FINAL.csv (v3)
./run_limited.sh $PY backtest.py                                             # does selection beat the plain ranker?
```

`rfd4_arm.py launch` registers the run in the ledger, yields at least 30 s to any queued `wait_and_run.sh` waiter,
waits for an idle card (one heavy job at a time), then runs the container with `--batch-size 1`. Use a new seed.

## Files

`her2lib.py` constants and contact table | `novelty.py` exact novelty index | `models.py` pre-filter as tables and
bootstrap | `designloop.py` one enumeration round | `run_loop.py` the loop with the ledger protocol |
`validate_abguard.py` the real checkers, rapidfuzz swapped in, 2 shards | `rfd4_arm.py`, `rfd4_harvest.py` RFD4 |
`tm_arm.py`, `adapt_arm.py` add-on arms | `final_portfolio.py` ranked list | `developability.py` liabilities, charge, pI |
`backtest.py`, `probe_abguard.py`, `feasible.py` measurements | `gpu_stage.py` NanoMelt only.

## Retired on the team's instruction: do not call

Boltz-2 (the Boltz functions in `gpu_stage.py`, `gpu_worker.py`) and ProteinMPNN (`developability_retired_mpnn.py.txt`;
its one run was also invalid, a trailing comment dropped the soluble-weights flag).

## What the measurements say

- A 13-mer CDRH3 needs 4 substitutions to fall under 70% identity. The binding gate is novelty: 3 or more mutations from
  every reference loop, and every measured SPR variant is in that reference, so none can be submitted.
- abguard counts indels; a Hamming index misses about 0.5%, so the real abguard screens every shortlist.
- Held-out backtest, legal band: the ranker's top pick beats trastuzumab's 0.152 nM 1% of the time and lands within 3x
  14% of the time. A bootstrap lower bound and an anchored estimate did not beat the plain ranker.
- Enumeration converges: round 2 enumerated 1.65 billion and added 34 designs.
- RFD4 her2-novel: the organisers report 4% known binders and 86% passing novelty, with no measured binding rate. In a
  12-design test, loops the classifier likes best failed abguard novelty and those that passed had p10 0.03 to 0.49.
