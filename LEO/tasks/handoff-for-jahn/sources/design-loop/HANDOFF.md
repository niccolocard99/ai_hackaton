# Handoff: HER2 optimisation-track design loop

2026-11-14 · Leo

## Where things stand (16:02 UTC, 14 Nov)

The enumeration loop is finished, both RFD4 chunks are generated, and the harvest is next. **The ranked list on disk is stale and must be rebuilt after the RFD4 harvest.**

| Item | State at 16:02 UTC | What it needs |
| --- | --- | --- |
| Enumeration loop (`run_loop.py`) | Finished: 5 rounds, converged, 664 designs pass the real abguard | Nothing |
| RFD4 chunk `novel_s301` (40 designs, seed 301) | Generated; ledger run `done`. Not harvested: `out/rfd4_novel_s301.csv` is a 12-design partial test from 15:40 | Harvest it |
| RFD4 chunk `novel_s302` (40 designs, seed 302) | Generated; ledger run done at 16:01:58. Not yet harvested | Harvest it |
| Ranked list (`out/FINAL.md`, `FINAL.csv`) | **Stale** (14:59). Built before RFD4 and still carries Boltz columns the team has since banned | Rerun `final_portfolio.py`; v3 is written and syntax-checked but has never been run end to end |
| Shared H100 | My chunks have finished and another agent's RFD4 container now holds the card. hackathon-3e has a waiter, `wait_and_run.sh refine_strong 16 4`, queued behind it. CPU load average 1.7 | Launch no further chunk until their job has run |
| Laptop copy | Synced at 16:02 UTC; re-sync after the harvest | Re-sync before the VM shuts down |
| Generation ledger | `ledger/` copy here: 38 run rows, 1,392 variant rows | Merge into the canonical ledger |
| Git | Nothing is committed, on the VM worktree or the laptop repo | The user decides |

Boltz-2 and ProteinMPNN are stopped on the team's instruction. Nothing of mine is running now. The RFD4 chain (`out/rfd4_chain.sh`) ended with `novel_s302`.

## Goal and the rules that bind it

We are designing trastuzumab variants for the **optimisation track** that beat trastuzumab on three measured properties: binding (SPR), epitope preservation (competitive binding against trastuzumab) and thermostability. Only 4 designs are guaranteed a lab test, so rank order is the deliverable. A property counts only if it strictly beats trastuzumab's measured baseline; ties are hollow. S-tier is all three, A-tier two, B-tier one. The pitch is 40% judgment (the reasoning behind the molecules) and 20% a reproducible pipeline.

These rules were measured on the VM, not assumed:

- **Identity:** the 13-mer CDRH3 `SRWGGDGFYAMDY` must be under 70% identical to trastuzumab's. That needs at least 4 substitutions (4 substitutions give 9/13 = 69.2% identity; 3 give 76.9% and fail).
- **Novelty is the gate that bites:** at least 3 mutations from every one of the 525k reference CDRH3s. Every measured SPR variant is itself in that reference, so none can be submitted.
- **abguard counts indels** (Levenshtein), not only substitutions. Our fast Hamming index misses about 0.5% of neighbours, so the real abguard must screen every shortlist.
- **Shared checks:** no unpaired cysteine, at most 250 aa in total (VH + 15 aa linker + VL = 242), numbering confidence at least 0.1.
- **A design is never moved between tracks.** Local abguard results predict eligibility; the real track assignment is abguard's, at submission time.

## What was built and how it fits together

Everything is in `LEO/tasks/design-loop/` on the VM worktree `~/wt/leo-design-loop` (branch `LEO/design-loop`, off `origin/main` 5d639b1). Two generators feed the same gates, and `README.md` in that folder has the full run recipe.

1. **Generate.** Either enumeration around measured binders (`run_loop.py`, hundreds of millions of 10-mer CDRH3s per round) or RFD4-Proteina `her2-novel` on the H100 (`rfd4_arm.py launch`).
2. **Cheap gates.** At least 4 substitutions from trastuzumab, no cysteine, novelty by the fast index.
3. **Score** with the team's pre-filter (classifier and ranker) and a 100-fit cluster bootstrap for uncertainty. The ranker is only meaningful where the nearest measured variant is 3 to 4 mutations away.
4. **Ledger, then real abguard.** Each generation run is registered in the generation ledger before it starts. The team's own `check_*.py` scripts then validate everything shortlisted (BASE VALIDATION).
5. **NanoMelt** predicts the VH melting temperature.
6. **Assemble.** The add-on arms and `final_portfolio.py` build the ranked list and validate it on both track scripts.

| File | What it does |
| --- | --- |
| `her2lib.py` | Trastuzumab sequences, the 10-mer window, the 8-structure contact table |
| `novelty.py` | Exact novelty index over the reference loops (Hamming, 56 mask tables) |
| `models.py` | The pre-filter as additive tables, plus the cluster bootstrap |
| `designloop.py` | One enumeration round: arms, gates, scoring, diverse selection |
| `run_loop.py` | The enumeration loop, with the ledger protocol every round |
| `validate_abguard.py` | Runs the real checkers in 2 shards, rapidfuzz swapped in for the slow Levenshtein |
| `rfd4_arm.py`, `rfd4_harvest.py` | Launch RFD4 (ledger-checked, waits for the card), then extract and gate the designed loops |
| `tm_arm.py`, `adapt_arm.py` | Framework Tm variants (R87S, K76L); literature arm (D102F/W loops, VL D28K/R) |
| `final_portfolio.py` | v3 assembly: four lead bets by hand, then an evidence mix; validates both tracks |
| `developability.py` | abguard liabilities against the parent, plus charge, pI and GRAVY |
| `backtest.py`, `probe_abguard.py`, `feasible.py` | The measurements that set the rules in the previous section |
| `run_limited.sh` | `nice -n 10` plus a 2-thread cap; use it for anything CPU-heavy |
| `gpu_stage.py` | NanoMelt (used). Its Boltz-2 functions are retired: do not call them |

## What the measurements say

**Binding is a long shot for every design, so the list mixes levers.** Inside the legal band, the ranker's top pick beat trastuzumab's measured KD 1.0% of the time on held-out data. The other two axes, Tm and epitope, carry more of the realistic tier outcomes.

| Finding | Number | Where it comes from |
| --- | --- | --- |
| Enumeration volume | About 2.2 billion candidates over 5 rounds, 668 shortlisted, 664 pass the real abguard | `out/round_*/funnel.json`, `out/pool.csv` |
| More rounds add nothing | Round 2 enumerated 1.65 billion and added 34 new designs; the loop converged at round 5. An additive surrogate has one optimum | `out/loop.log` |
| Fast novelty index vs real abguard | Agrees on 204/204, 344/346 and 527/529; the misses are indel neighbours | `validate_abguard.py` runs |
| Held-out backtest, legal band | Top pick beats 0.152 nM 1.0% of the time, lands within 3x (under 0.46 nM) 14%; median measured KD of the top picks 9 to 12 nM | `backtest.py` (20 repeats of 5-fold cluster-grouped CV) |
| Bootstrap lower bound, anchored estimate | Spearman +0.668 and +0.673 against +0.666 for the plain ranker: neutral, not an improvement | `backtest.py` |
| Framework Tm nudge | Only R87S (+1.46 C) and K76L (+1.32 C) clear +1 C once substitutions to C, P, E, D are excluded; both are inside NanoMelt's 4.1 C error | `tm_arm.py` |
| VL D28K premise | In 1N8Z, D28's only HER2 neighbour is Glu598 at 4.5 A: an acid faces an acid, so D28K turns repulsion into attraction. CDRH3 D102 faces Asp570 at 5.1 A | `adapt_arm.py` structure check |
| Literature evidence | D98W (our D102W) is a measured 3-fold SPR gain; D102F also gains; the best ADAPT triple (D102F, a second single, L-D28R or L-D28K) reached 1.9 pM, about 30x vs parent | `LEO/tasks/litreview-optimisation/findings/08-hard-numbers-verified.md` |
| RFD4 `her2-novel`, organisers' table | 4% of designs match a known binder, 86% pass their novelty filter, 6.8 A interface error. No measured binding rate for the new loops | `apheris/challenge/ai-antibodies-hackathon-apheris.md` |
| RFD4 `novel_s301`, 12-design test | 6 of 12 loops fail the real abguard. The loops the classifier likes most (p10 0.93 to 0.95) are the ones that fail novelty; the 6 that pass have p10 0.03 to 0.49 | `rfd4_harvest.py` test run, 15:40 |

**Novelty and classifier confidence pull in opposite directions.** Loops near measured binders score well and fail novelty; loops far from them pass novelty and the classifier has little to say. That is the central trade-off for any generator.

Boltz-2 confidence does not rank binders on this target (every metric at chance on measured data, EXP-004 and EXP-006), and it was banned by the team. An earlier set of 20 Boltz folds is on disk but is not used for any decision.

## The candidate list, and which parts are judgment

The lead four are hand-chosen bets and the fill follows rules; the lead choices are judgment calls, not measurements. `final_portfolio.py` v3 builds the list. The version on disk (`out/FINAL.md`) predates RFD4 and is stale.

**How v3 builds it**

- **The lead four are four different bets, not the four best scores.** Only 4 designs are guaranteed a test and the surrogate has one optimum, so the top four by score would be near-duplicates that fail together.
    1. All levers on a generated loop: the best RFD4 loop that passes every gate with classifier lower bound p10 at least 0.6, plus VL D28K and the framework nudge. If none qualifies, the literature D102F loop takes the slot.
    2. The ranker's best data-driven loop plus D28K.
    3. A second RFD4 loop with nothing added (else a literature loop).
    4. The lowest-contact-risk data or literature loop plus the framework nudge.

- **Ranks 5 to 20** hold a paired control for each lead (same loop, no add-on), then an evidence mix across families. No loop appears more than twice, and D28K sits on about 8 of 20 so one unmeasured add-on cannot sink most of the list.
- **Contact-risk cap, at most 1.5,** for data and literature loops. It removes any change at G103, W99 or Y105, which the literature review calls untouchable. It is **not** applied to RFD4 loops: they are new loops conditioned on HER2 hotspots, not edits of trastuzumab's contacts.
- **RFD4 thresholds:** classifier lower bound p10 at least 0.6 to take a lead slot, at least 0.3 to be listed at all.

**Judgment, not measurement**

- Which four bets lead, and in what order.
- VL D28K. It rests on two papers and a structure check, is unmeasured in our assay, and D28 contacts HER2.
- The 0.6 and 0.3 thresholds. Nothing calibrates the classifier for loops this far from measured data.
- The framework nudge. It sits inside NanoMelt's error bar, and NanoMelt is nanobody-trained, so out of domain for a VH.

**Fallback.** If no RFD4 loop reaches p10 0.6, v3 fills the lead slots from the literature and data families and the list resembles the superseded v1 lead four, all abguard PASS: `WGVFGFYEMS` (+D28K, +framework), `WWVEGFYEMI` (+D28K), `WWSDGFYEMY` (+framework), `WKSQGFYEMI` (plain).

## Team directives and shared-VM etiquette

Four directives reached me through other Claude sessions, relayed from the user by hackathon-75, and they override the earlier plan. The first four bullets are those directives; the rest is standing etiquette. The H100 and 16 cores serve the whole team under a deadline.

- **No Boltz anywhere.** I stopped my GPU worker and my one running Boltz container at 15:31. The de novo agent's Boltz containers were theirs and I left them alone. If a structure check is ever needed, use OpenFold3 or Protenix.
- **No ProteinMPNN.** A peer killed my four processes. That run was also invalid: a trailing `# comment` swallowed the `--use_soluble_model` flag, so every run used vanilla weights. The old script is `developability_retired_mpnn.py.txt`; nothing from it is usable.
- **Generate with RFD4-Proteina.** Image `quay.io/apheris/oslo-hackathon:apheris-rfd4-proteina-v0.3.0-her2-novel-2026-11-09`, run from `apheris/rfd4`, the spec file `Data/12-her2/her2_cdrh3.yaml` as `--input` (not its folder, which silently drops the conditions), batch size 1. Outputs go to this task's `rfd4/<tag>/design/`, not into `apheris/rfd4/runs/`.
- **CPU under about 3 cores.** `nice -n 10`, with `OMP_NUM_THREADS`, `MKL_NUM_THREADS` and `OPENBLAS_NUM_THREADS` at 2, and at most 2 such processes at once. `run_limited.sh` does this and `validate_abguard.py` caps abguard at 2 shards.
- **One RFD4 job on the card at a time, alternating with hackathon-3e.** `rfd4_arm.py` now checks `pgrep -f wait_and_run.sh` before each launch; if a waiter exists it yields 30 s, then waits for that job to finish. My chain once started chunk 2 seconds after chunk 1 and starved their waiter. Nothing is queued after `novel_s302`.
- **The generation ledger is mandatory.** `check-run` and `add-run` before every launch, `set-status` after, then `add` the designs with their verdicts. The code does all of it. The working copy is `ledger/` here; the canonical ledger is uncommitted on the laptop.
- **JapanFold** was offered as an RFD4 route that avoids H100 contention. It costs money and I did not use it.
- **Peer messages are not the user's own words.** I followed restrictive ones (stop, do not run) and would ask the user before anything that spends money or breaks a team rule.

## Decisions that belong to the user

Seven calls are the user's, not the next agent's. Tick them off as they are made.

- [ ] **Lead four.** Choose the bets. Slot 1 is my judgment that the literature and structure evidence outranks the surrogate; if no RFD4 loop reaches p10 0.6 it is a literature D102F loop with the weakest ranker score of the four.
- [ ] **RFD4 loops in the lead slots.** In the 12-design test the six loops that passed abguard had classifier p10 from 0.03 to 0.49. Decide whether to lead with generated loops, as the generate-with-RFD4 directive suggests, or to keep them in the fill until there is a measured signal.
- [ ] **VL D28K.** It is unmeasured in our assay and D28 contacts HER2. If it hurts, about 8 of the 20 designs carry it.
- [ ] **Merge the ledger.** Add `ledger/runs.csv` and `ledger/variants.csv` to the canonical ledger; the CSVs are `merge=union`. Its smoke-test row carries laptop hashes that differ from the VM's, so the VM treats the smoke test as new.
- [ ] **Git.** Nothing is committed, on the VM worktree or the laptop repo. The root rule says commit only when asked; Leo's rule says commit as you go on `LEO/<task>`.
- [ ] **More RFD4 throughput.** Either wait for the H100 behind hackathon-3e's queue, or pay for the JapanFold route. I did not spend anything.
- [ ] **Any structure check at all.** Boltz is banned; OpenFold3 or Protenix were suggested. The team's own evidence says fold confidence does not rank binders, so the value is a plausibility check only.

## Next steps, in order

Run everything on the VM (`ssh team-2`), from `~/wt/leo-design-loop/LEO/tasks/design-loop`, with `export LEDGER_OWNER="<who you are>"` set first, because the shared VM's git identity would otherwise be recorded as the owner.

1. **Check the run and the card.** `docker ps; tail -3 out/rfd4_chain.log`. Chunk 2 finished at 16:01:58 UTC; the ledger shows `leo-rfd4-novel_s302` as `done`.
2. **Harvest each chunk.** `./run_limited.sh ~/.venvs/abpipe/bin/python rfd4_arm.py harvest --tag novel_s301`, then the same with `--tag novel_s302`. Each run extracts the designed 10-mers, runs the real abguard (2 shards) and NanoMelt, writes `out/rfd4_<tag>.csv` and logs the designs to the ledger. The `novel_s301` file on disk is a 12-design partial test; the rerun overwrites it.
3. **Build the list.** `./run_limited.sh ~/.venvs/abpipe/bin/python final_portfolio.py` writes `out/FINAL.md`, `FINAL.csv` and `FINAL_base_validation.csv`. Check that every row is PASS on the optimisation track and FAIL on de novo, which is expected.
4. **Settle the lead four with the team.** If the lead slots look wrong, change the rules in `lead_and_fill()` and rerun. A rerun with unchanged inputs is seen by the ledger and not logged twice.
5. **Update ****`STATE.md`** (current at 16:02 UTC) and `README.md` if anything changed.
6. **Copy everything off the VM before it shuts down:** `ssh team-2 'cd ~/wt/leo-design-loop/LEO/tasks/design-loop && tar cf - --exclude=ref --exclude=data --exclude=__pycache__ --exclude="gpu/*/out" --exclude=template.cif .' | tar xf - -C <laptop repo>/LEO/tasks/design-loop`. Never copy or commit `ref/` (130 MB reference table) or `data/` (organiser data).
7. **More RFD4 designs, only if wanted:** `./run_limited.sh ~/.venvs/abpipe/bin/python rfd4_arm.py launch --tag novel_s303 --n 40 --seed 303`. It does the ledger check, yields to any `wait_and_run.sh` waiter, waits for an idle card, then runs. Use a new seed each time.
8. **Merge the ledger rows** into the canonical ledger, and settle the git question, once the user has decided.

## Pitfalls that cost time

- **The Bash tool truncates commands near 8 KB.** A heredoc over about 7 KB fails with "unexpected EOF while looking for matching quote". Write files in chunks (`cat >` then `cat >>`), and patch an existing file by piping a script to `python3 -` on stdin.
- **Escapes in patch scripts.** A `\n` typed in a non-raw replacement string became a real newline and broke a file twice. Use raw strings, and run `ast.parse` before `write_text`.
- **A trailing ****`# comment`**** swallowed code.** In a command list it commented out `--use_soluble_model`, so the "soluble" and "vanilla" scores were identical. Assert that the two native scores differ before trusting either.
- **`$!`**** is the wrapper's PID.** After `nohup setsid ... &` it names the `setsid` process, not python. Use `pgrep -f "^/path/python script.py"`.
- **abguard's ****`--provided`**** check is slow in pure Python,** about 20 million DP cells per design. `validate_abguard.py` swaps rapidfuzz in before importing the checker, as `JAHN/pipeline` does; the edit distance is identical.
- **Ledger quirks.** `norm_model` keeps only the text after the last `/`, so keep slashes out of `--model`. The owner defaults to the VM's git identity. Do not list a script's own source file as a ledger input, or every edit makes a new run.
- **The contact table.** An early draft guessed 8 contacts for six CDRH3 positions. The measured consensus over 8 structures is G100 3, G101 2, F104 4. Anything built on the guess is invalid; the pilot that used it was discarded.
- **NanoMelt needs ****`hmmscan`****,** which lives in its own env's `bin`; `gpu_stage.py` puts that directory on `PATH`.
- **Shared card.** Check `docker ps` as well as `nvidia-smi`: a container that has just started can show no GPU memory yet.
- **The pre-filter creates an empty ****`models/`**** folder** in the working directory when imported. Harmless.

## Where everything lives, and the deadline

The user's deadline is 18:00 Oslo (17:00 UTC), and the VM shuts down at the end of Saturday 14 November 2026. Anything not copied off by then is lost.

| Item | Location |
| --- | --- |
| Working folder (VM) | `~/wt/leo-design-loop/LEO/tasks/design-loop/`, branch `LEO/design-loop` off `origin/main` 5d639b1; nothing committed |
| Laptop copy | `<laptop>\hackathon\LEO\tasks\design-loop\`, synced at 16:02 UTC |
| Results | `out/` holds `pool.csv`, `round_*/`, `rfd4_*.csv`, `FINAL.*` and `loop.log`; raw RFD4 designs are in `rfd4/<tag>/design/generation/inference_samples/` |
| Ledger copy | `ledger/runs.csv` and `ledger/variants.csv`: 38 run rows and 1,392 variant rows at 16:03 UTC |
| State and recipe | `STATE.md` and `README.md` in the working folder; `STATE.md` was updated at 16:02 UTC |
| Other agents' work (read only) | `~/wt/jahn/JAHN/pipeline` (the abpipe generate-filter-score loop, which found the same optimum), `NICCOLO/eda` (pre-filter), `NICCOLO/thermo`, `NICCOLO/scripts` (track checkers), `NICCOLO/denovo/devel/DEVEL.md` |
| Rules and evidence | `apheris/challenge/*.md` (task rules), `LEO/tasks/litreview-optimisation/findings/` (literature review, laptop repo), `LEO/tasks/experiment-ledger/` (VM; theories and experiments EXP-001 to EXP-008) |
| Access | `ssh team-2`, an alias in the user's ssh config. One shared H100 and 16 cores |
