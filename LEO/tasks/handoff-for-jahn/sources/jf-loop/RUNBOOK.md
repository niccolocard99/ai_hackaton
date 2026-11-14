# Runbook: operating the JapanFold design loop

**As of:** 2026-11-14, about 16:00 UTC · **Written for:** an agent or teammate who must run, stop, restart, verify or shut down this system. Overview and open decisions are in [HANDOFF.md](HANDOFF.md); evidence is in [FINDINGS.md](FINDINGS.md). All commands run **on the VM** (`ssh team-2`) unless marked laptop.

## 1. Where things are (VM, under `/home/hackathon`)

| Path | What |
|---|---|
| `leo-work/h100/jf_loop.py` | **The loop (v3).** Imports `h100_loop.py` for the gates, scoring and ledger wrapper, and points its globals at `leo-work/jf/` |
| `leo-work/h100/h100_loop.py` | Shared library (harvest gates, `score_new`, `table`, `ledger()` with its lock, `pick_parents`). Also the **stopped** H100 loop: its `main()` must not be run |
| `leo-work/h100/{specs.py,lib.py}` | `lib.py`: mmCIF parsing, CDRH3-HER2 contacts, epitope footprint. `specs.py`: the H100 Docker spec generator (not used by the API loop) |
| `leo-work/h100/fast_check.py` | Team optimisation checker with a rapidfuzz Levenshtein. Verified identical (HANDOFF open item 1) |
| `leo-work/jf/` | The loop's working dir: `runs/cNN/<job>/` (`job.json`, `in/`, `design/*.cif`, `designs.csv`, `fresh.csv`), `pool.csv`, `shortlist.csv`, `state.json`, `abguard/<tag>/opt.csv`, `loop.log`, `Data/her2_fv_1n8z.cif` |
| `leo-work/finalize.py`, `leo-work/compare_lists.py` | Shortlist + both-track BASE VALIDATION; three-list comparison + S97 hedge |
| `leo-work/final/` | Outputs: `candidates.csv`, `designs.csv`, `reports/`, `comparison.csv`, `s97_hedge.csv` |
| `leo-work/generation-ledger/` | Ledger tool, `runs.csv`, `variants.csv`, `.lock`. **Shared with the `denovo-vhh-epitope` session** |
| `leo-work/prefilter/` | `data/` (SPR-768 table), `scripts_main/` (checkers copied from `NICCOLO/scripts/`), NICCOLO's `prefilter.py` |
| `leo-work/loop/prefilter.pkl` | The pickled fitted `PreFilter` the loops load (needs `leo-work/prefilter` on `sys.path`) |
| `leo-work/loop/`, `leo-work/japanfold/` | The first (CLI, unledgered) JapanFold runs: kept for provenance, backfilled into the ledger as RETROACTIVE |
| `~/.venvs/abmodels` | Python for the loops (pandas, sklearn, numpy) |
| `~/.venvs/abguard` | abguard 0.1.0 (plain team checker) |
| `~/.venvs/abpipe` | abguard + rapidfuzz. **Another session's venv: execute it, never modify it** |
| `~/.config/japanfold/api_key` | The API key, mode 600. Never print it |

Repo (laptop): `LEO/tasks/rfd4-japanfold/` holds the documents, `snapshot/` (a copy of the above taken at about 15:58 UTC) and `tools/`.

## 2. The loop: start, stop, check

```bash
# START (safe to re-run: state is in leo-work/jf/state.json, and reconcile() closes any run left 'started')
cd ~/leo-work/h100 && rm -f ~/leo-work/jf/STOP
OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 POLARS_MAX_THREADS=2 RAYON_NUM_THREADS=2 \
  nohup nice -n 10 ~/.venvs/abmodels/bin/python -I jf_loop.py >> ~/leo-work/jf/loop.log 2>&1 &

# STOP cleanly: finishes the cycle in flight, then does a final ingest and exits
touch ~/leo-work/jf/STOP

# CHECK
ps -eo pid,ni,etime,pcpu,args | awk '/jf_loop/ && !/awk/ && !/ps -eo/'   # alive? (see the pgrep trap in section 9)
grep -E "JF CYCLE|plan:|submitted|jobs finished" ~/leo-work/jf/loop.log | tail
```

Options: `--max-cycles N`, `--until HH:MM` (VM clock; default none), `--ingest-only`, `--dry-run` (prints a plan and submits nothing; **the user asked not to use it**).

**One cycle** (about 6.5 min): `reconcile()` closes stale runs, then `plan()` reads `/v1/cluster` and sizes up to 40 jobs (refine jobs first, then fresh arms allocated by recent yield), `submit_all()` does ledger `check-run` + `add-run` and a `POST /v1/designs` per job, then **ingests the previous cycle's results while these run**, then waits and downloads. Ingest is one batched validation: `ledger check --new-only` per job, the fast checker once, `ledger add` per job. The loop stops on `STOP`, `--until`, or an HTTP 402.

**Spend.** There is no cap (the user's instruction, relayed). The balance is logged each cycle with `$/min`. To re-add a guard, in `main()` after `now = balance()` add `if bal0 is not None and now is not None and bal0 - now > CAP: log("budget guard"); break`. The v2 loop had exactly this.

## 3. Arms

Pins are encoded in the **contig**: a fixed segment (`B105-105`) keeps that residue's identity **and coordinates**; a designed segment (`3-3`) fixes only the length, so the loop stays 10 long and abguard's 13-residue CDRH3 keeps its length. R98 is always fixed (it salt-bridges HER2 D585). Arms are generated from `ARM_SPECS` in `jf_loop.py`.

| Arm | Pinned | S97 designed | Note |
|---|---|---|---|
| `jy` | Y105 | no | the base Y-fixed arm; collapses to few distinct sequences |
| `jy97` | Y105 | yes | |
| `jyf`, `jyf97` | F104, Y105 | no / yes | |
| `jyg`, `jyg97` | G103, Y105 | no / yes | |
| `jym` | Y105, M107 | no | |
| `jygf`, `jygm` | G103+F104+Y105 / G103+Y105+M107 | no | |
| `jw97` | W99, Y105 | yes | the most productive anchor-rich arm so far |
| `jwy`, `jan4` | W99+Y105 / W99+G103+F104+Y105 | no | probes only, every 4th cycle: **0 legal in 201 designs** from this generic RFD4 |
| `jfree`, `jfree97` | nothing (R98 only) | no / yes | controls, capped at 1 and 2 jobs; they lose Y105 |

Example contigs (chain A, B, C, D; JapanFold uses the 1N8Z numbering of the input file):
```
jy    A1-107,/0,B1-98,6-6,B105-105,3-3,B109-120,/0,C489-580,/0,D591-607
jy97  A1-107,/0,B1-96,1-1,B98-98,6-6,B105-105,3-3,B109-120,/0,C489-580,/0,D591-607
```
**Refinement is inpainting**: take a top design's own CIF as the structure, keep every residue except 3-5 random non-pinned loop positions, and redesign those (`refine_contig`). JapanFold renumbers HER2 chains from 1 in its outputs, so a refine contig ends `C1-92,/0,D1-17`, and `harvest_jf` adds +488 / +590 to compare with 1N8Z numbering.

To add an arm: add a line to `ARM_SPECS` (`name: ({position: residue}, design97)`). `arm_contig()` builds the contig, `plan()` asserts chain B is 120 residues, `allocate()` gives it an exploration job automatically.

## 4. JapanFold API cheat sheet

Base `https://api.japanfold.aiand.com`, header `Authorization: Bearer <key>`.

| Call | Use |
|---|---|
| `POST /v1/designs` | body `{"protocol":"rfd4-scaffold","structure":"<mmCIF text>","contig":"...","name":"...","params":{"num_designs":8,"recipe":"high/front_loaded","nsteps":60,"seed":N}}`, header `Idempotency-Key`. 202 = job. 402 = no balance. 429 = rate limit or busy queue. 503 = no galaxy serves the model. 422 = key reused with a different body |
| `GET /v1/jobs/{id}?wait=30` | long-poll until the job ends |
| `GET /v1/jobs/{id}/archive` | zip of `design_k_generated.cif` + metrics. Validate with `zipfile.testzip()` before extracting: a truncated body raised `IncompleteRead` once |
| `GET /v1/cluster`, `GET /v1/credits` | fleet load; balance |

`hotspots` is PXDesign-only in this API. `num_designs` is at most 8. The ledger `run_key` is used as the `Idempotency-Key`, so a retried POST cannot create a second paid job. Measured cost about $0.015 per 8-design job (38 jobs = $0.58). The Boltz-2 prediction endpoint takes no template, so its folds are not comparable to the organisers' templated Docker folds.

## 5. Gates and the pool

`eligible` = abguard optimisation PASS **and** no Cys **and** no N-glycosylation sequon **and** pins held **and** framework intact **and** docked (at least 30 CDRH3-HER2 heavy-atom contacts in the design's own complex; the parent has 53). Tiers (`finalize.py`): tier 1 = Y105 kept, p_binder >= 0.7, contacts >= 0.8x parent, epitope overlap >= 0.7, **and** a measured binder within 4 substitutions; tier 2 = the same without that local evidence; tier 3 = Y105 kept, docked, epitope ok, p_binder < 0.7; tier 4 = the rest. Within a tier: p_binder, then novelty margin. **None of this predicts affinity** (FINDINGS section 9).

`pool.csv` columns that matter: `design_id`, `run_id` (the ledger run), `arm`, `mode`, `parent_id`, `cdrh3` (10-mer, 99-108), `cdrh3_13` (the pool key, 97-109), `n_mut`, `min_dist` (abguard novelty distance), `p_binder`, `nn_dist`, `n_bind_le4`, `best_kd_le4`, `contact_ratio`, `epi_overlap`, `eligible`, `cif`, `heavy_sequence`.

## 6. Validation, CPU etiquette, and the patched abguard

The VM is shared (16 cores). **Always** run validation throttled and niced, one at a time:
```bash
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 POLARS_MAX_THREADS=2 RAYON_NUM_THREADS=2
nice -n 10 ~/.venvs/abpipe/bin/python ~/leo-work/h100/fast_check.py designs.csv -o report.csv \
    --provided ~/leo-work/prefilter/data/trastuzumab_her2_spr_all_768.csv
```
`fast_check.py` is the team's `check_optimisation.py` with one change (rapidfuzz Levenshtein); on 60 designs (32 of them rejected) it matched the standard checker's stored verdicts on passes, distance and identity, 60/60, and took 2.65 s against about 90 s. The plain checker spawns 86 threads: do not use it unthrottled.

**Patched abguard (open; the user asked for it).** `~/nick-1m/pipeline/abpipe/abguard_fast.py` replaces abguard's novelty cross-join with an index. To integrate, in `fast_check.py` before `import check_optimisation`:
```python
sys.path.insert(0, str(Path.home() / "nick-1m/pipeline"))
from abpipe import abguard_fast
abguard_fast.install()        # idempotent; uninstall() restores the original
```
Then (1) run it on `final/compare/fc_in.csv` (60 designs, 32 rejected) and compare with the stored verdicts; expect identical `passes`, and identical `min_dist` for rejects, with comfortable passes reported as `cap+1` (default 4); (2) **time it including the index build**: the index is built per process, so for about 300 designs per cycle it may not beat the 26 ms/design cross-join unless the process is kept alive; (3) `ABGUARD_FAST_MAX_DISTANCE=N` raises the reported cap at a memory/time cost. Margin tiering only needs `>= 4`, but the pool would lose true distances of 5 and 6.

## 7. The generation ledger

The root `AGENTS.md` rule applies to every run. The loops do it automatically; for a manual run (always from the folder holding the input files):
```bash
L="python3 ~/leo-work/generation-ledger/ledger.py"; D="--ledger-dir $HOME/leo-work/generation-ledger"; export LEDGER_OWNER=LEO-claude
$L check-run $D --method design --model japanfold:rfd4 --input Data/her2_fv_1n8z.cif --seed 7 --setting protocol=rfd4-scaffold --setting contig=... --setting num-designs=8 --setting recipe=high/front_loaded --setting nsteps=60   # exit 3 = already run
$L add-run   $D ...same flags... --run-id <unique> --output <dir>
$L set-status $D --run-id <unique> --status done|failed|empty
$L check $D designs.csv --new-only fresh.csv ; $L add $D designs.csv --run-id <unique>
```
Wrap manual calls in `flock ~/leo-work/generation-ledger/.lock` (the loops do), because the VM copy is shared with another session. Rows are append-only: a changed status is a new row.

**Append the VM ledger back to the repo ledger** (laptop; only rows not already present, existing rows untouched): `python LEO/tasks/rfd4-japanfold/tools/sync_ledger_to_repo.py`.

## 8. Shutdown checklist (the VM dies at the end of 14 Nov)

1. `touch ~/leo-work/jf/STOP`; wait for `loop ended` in the log (final ingest included).
2. `cd ~/leo-work && ~/.venvs/abmodels/bin/python -I finalize.py 20 && ~/.venvs/abmodels/bin/python -I compare_lists.py` (both throttle their own subprocesses; the system `python` lacks pandas).
3. Build and pull the snapshot: a copy of `tools/make_snapshot.py` is already on the VM at `~/leo-work/make_snapshot.py` (otherwise `scp` it there). Run `~/.venvs/abmodels/bin/python -I ~/leo-work/make_snapshot.py` (it writes `/tmp/leo_snapshot.tgz`), then on the laptop `scp team-2:/tmp/leo_snapshot.tgz` into `LEO/tasks/rfd4-japanfold/`, `rm -rf snapshot`, and extract (it creates `snapshot/`). Scan for credentials before committing: `grep -rIl -E "jf_live_|amass_[A-Za-z0-9]{8}|Bearer [A-Za-z0-9]{12}" snapshot`.
4. Run `tools/sync_ledger_to_repo.py` (laptop).
5. Also copy the sibling's `out/` (`~/wt/leo-design-loop/LEO/tasks/design-loop/out/`) if its owner has not: it is not mine, so ask first.
6. Commit only if asked.

## 9. Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `IncompleteRead` in a download | truncated chunked body; `HTTPException` is not an `OSError` | fixed in v3 (retry + zip test); `reconcile()` re-downloads |
| ledger runs stuck at `started` | a crash or an error that escaped | `reconcile()` runs at every cycle start and closes them; mine are only ever `started` while in flight |
| `pgrep -f` / `ps | grep` says the loop "is still running" after it exited, or `pkill -f` kills your own shell | the pattern also appears in **your own SSH command line** | use `ps -eo pid,args | awk '/jf_loop/ && !/awk/ && !/ps -eo/'`, kill by PID |
| VM load high | unthrottled validation (86 threads) | `ps -eo pid,nlwp,pcpu,args --sort=-pcpu | head`; use section 6 |
| HTTP 429 / 5xx on submit | shared fleet is busy | built-in exponential backoff; `plan()` also shrinks the batch from `/v1/cluster` |
| HTTP 402 | no balance | the loop stops by design |
| a large file written by SSH heredoc fails to parse | the command string is too big for the shell | write the file locally and `scp` it |
| `LF will be replaced by CRLF` from git | Windows autocrlf | harmless |
| an arm gets no jobs | `allocate()` weights by recent yield; probes run every 4th cycle | check `st["gen"]` in `state.json` |
