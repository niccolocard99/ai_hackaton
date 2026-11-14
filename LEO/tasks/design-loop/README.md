# design-loop: optimisation-track generate, gate, score and fold loop

Runs on the team VM. Everything here was executed there; nothing was run on a laptop.

```
enumerate (hundreds of millions of CDRH3 10-mers per round, arms)
  -> hard gates: >= 4 substitutions, novelty >= 3 from all 525,042 reference loops, within reach of measured data
  -> NICCOLO pre-filter (classifier + ranker) with cluster-bootstrap uncertainty
  -> GENERATION LEDGER check, then REAL abguard track checkers (BASE VALIDATION)
  -> H100: NanoMelt Tm, Boltz-2 pose check (separate worker, only when the card is free)
  -> pool -> ranked list
```

## Run it

```bash
cd LEO/tasks/design-loop                      # venv: ~/.venvs/abpipe (abguard, rapidfuzz, sklearn, pandas, polars)
mkdir -p data ref
cp ~/hackathon-team-3/apheris/hackathon-data/her2-trastuzumab/trastuzumab_her2_spr_all_768.csv data/
unzip ~/hackathon-team-3/apheris/hackathon-data/her2-trastuzumab/tz_her2_facs_3classes_524346.csv.zip -d ref/
cp <generation-ledger>/{ledger.py,runs.csv,variants.csv} ledger/    # stdlib-only tool from LEO/tasks/generation-ledger
export LEDGER_OWNER="<who you are>"           # the shared VM git identity would be recorded otherwise
python run_loop.py --until 15:50 --max-rounds 8 &      # CPU loop; resumes from out/; create out/STOP to end it
python gpu_worker.py &                                   # folds only when nvidia-smi shows the card idle
python tm_arm.py            # framework Tm variants (needs out/portfolio.csv)
python adapt_arm.py         # literature arm: D102F/W loops, VL D28K/R
python final_portfolio.py   # ranked list -> out/FINAL.csv, out/FINAL.md, out/fold_queue.csv
python backtest.py          # does the selection beat the plain ranker on held-out variants?
```

Paths to NICCOLO/eda, NICCOLO/scripts and the Boltz reference inputs are absolute VM paths at the top of each
module; they are the one thing to change to run elsewhere.

## Files

| File | What it is |
|---|---|
| `her2lib.py` | Trastuzumab sequences, the 10-mer window, the 8-structure contact table |
| `novelty.py` | Exact novelty index (Hamming, 56 mask tables) over the reference loops |
| `models.py` | The pre-filter as additive tables, and the cluster bootstrap |
| `designloop.py` | One generation round: arms, gates, scoring, diverse selection |
| `validate_abguard.py` | Runs the team's real `check_*.py` in shards, rapidfuzz swapped for the Levenshtein |
| `gpu_stage.py`, `gpu_worker.py` | NanoMelt and Boltz-2 on the H100, with the free-card wait |
| `run_loop.py` | The loop, with the ledger protocol per round |
| `tm_arm.py`, `adapt_arm.py`, `final_portfolio.py` | The add-on arms and the final list |
| `backtest.py`, `probe_abguard.py`, `feasible.py` | Measurements that set the design rules |
| `ledger/` | Copy of the ledger tool and the rows this task added |

## What the measurements say (so nobody re-derives them)

- A 13-mer CDRH3 `SRWGGDGFYAMDY` needs 4 substitutions to fall under 70% identity (4/13 = 69.2%).
- The binding gate is novelty: 3 or more mutations from **every** reference loop. Every measured SPR variant is
  itself in the reference, so none can be submitted. abguard counts indels (Levenshtein); a Hamming index misses
  about 0.5% of neighbours, so real abguard screens every shortlist.
- Held-out backtest, legal band: the ranker's top pick beats trastuzumab's 0.152 nM 1% of the time, lands within 3x
  14% of the time. The bootstrap lower bound and an anchored estimate did not beat the plain ranker.
- A search past the first round finds nothing new: round 2 enumerated 1.65 billion candidates and added 34 designs.
