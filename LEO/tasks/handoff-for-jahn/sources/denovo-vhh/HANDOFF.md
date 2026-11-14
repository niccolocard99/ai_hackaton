# Handoff: de novo VHH binders on the trastuzumab epitope

Written 2026-11-14, about 16:00 UTC (the VM clock). Everything below was verified at that time unless it says
otherwise. Read this first, then [`STATE.md`](STATE.md) (live state) and [`NOTES.md`](NOTES.md) (the reasoning
behind every design choice).

## 1. The task in one paragraph

Design new VHH (nanobody) binders to HER2 domain IV that sit on **trastuzumab's epitope**, so they compete with
trastuzumab in Adaptyv's epitope-preservation assay. They are entries for the **de novo track** (under 70%
identity to any known antibody, at least 2 CDRH3 mutations from trastuzumab). Generation uses
**RFD4-Proteina, antibody weights** on the team VM's H100. Only 4 designs per team are guaranteed a lab slot,
so the real output is a short, well-argued ranked list.

## 2. Standing instructions. Follow these.

They come from Leo, directly or relayed by teammates' Claude sessions.

| Rule | Detail |
|---|---|
| **RFD4-Proteina only** | No Boltz, no Protenix or OpenFold3 refolds, no other model on the GPU. Leo's latest instruction. |
| **One RFD4 job on the H100 at a time** | The card is shared with other agents' `her2-novel` runs. Wait for it; do not set `FORCE=1` unless Leo says so. |
| **Stay under about 3 CPU cores** | `common.sh` already sets `docker --cpus 3` and thread caps (2). Run `abguard`/validation with `OMP_NUM_THREADS=2`, `nice -n 10`, at most 2 at once. The VM hit load 70 on 16 cores earlier. |
| **The epitope is fixed** | Same four HER2 hotspot atoms (D593 NZ, F573 CD1, I591 CG2, P572 CD) and the same HER2 domain IV target. Do not move the VHH to another epitope. |
| **Ledger before every generation run** | `run_design.sh` does `check-run`, `add-run` and `set-status`. Ledger: `/home/hackathon/leo-work/generation-ledger/ledger.py` (pass `LEDGER=`, `LEDGER_OWNER=`). |
| **Base validation before ranking** | `NICCOLO/scripts/check_denovo.py` via `validate_and_log.sh`. A local PASS is necessary, not sufficient. |
| **No thermostability claim from NanoMelt** | A VHH has no comparable baseline against trastuzumab. Tie-break among our own designs only. |
| **Do not commit, push or merge unless Leo asks** | Leo asked for exactly two things so far: the first commit (`68f0ce3`) and merging `origin/main` (`215f1ed`). Never print keys. Nick's folder (`~/nick-work`) is read only. |

## 3. Where things stand

| Item | State |
|---|---|
| Batch 1 (`fixedfw`, 100 designs) plus 4 smoke designs | Generated, validated, logged to the ledger, scored |
| **Round 2** | **Queued, not started.** `queue_round2.sh` is waiting on the VM for the GPU, which another agent's `her2-novel` container held at 16:00 (8 minutes into its run). Nothing from round 2 exists yet |
| Branch `LEO-denovo-vhh-epitope` | `68f0ce3` initial setup, `215f1ed` merge of `origin/main` (clean). **18 files are uncommitted**, including all scripts added after the first commit, `results/` and the notes |
| PR #17 (EXP-010, EXP-011, IDEAs 018 to 020) | On `origin/nick`, **not on main**. I read it without merging |

Round 2 is two steps, in order:
1. **Refine** 5 strong-epitope parents by RFD4 partial diffusion (fraction 0.4, 16 children each, 80 designs).
   The parent's annotations carry the fixed HER2 and the hotspots, so the epitope stays fixed by construction.
   Parents: `fixedfw_48`, `fixedfw_89`, `freefw_0`, `fixedfw_42`, `fixedfw_18`.
2. **Fresh batch**, `fixedfw_long` (CDR3 constrained to 17 to 21), 100 designs, seed 3.

## 4. Where everything is

Laptop worktree: `<laptop>\hackathon\.worktrees\denovo-vhh\LEO\tasks\denovo-vhh-epitope\`.
VM: `ssh team-2`, then `~/hackathon-team-3/LEO/tasks/denovo-vhh-epitope/` (the VM checkout is on branch
`main-loop` and holds this folder as an **untracked directory**: remove it before any `git pull` there).

| File | Purpose |
|---|---|
| `build_inputs.py` | Builds `spec/` from organiser files: target CIF, framework CIF, three specs (`fixedfw`, `freefw`, `fixedfw_long`). Deterministic |
| `check_spec.py` | Validates a spec against its structures and against 1N8Z contacts; run before any GPU time |
| `run_design.sh`, `common.sh` | One design run: spec check, ledger, GPU guard, docker. `PARTIAL=0.4` switches to refinement |
| `wait_and_run.sh`, `queue_round2.sh` | Queue behind other jobs; the round-2 plan |
| `make_refine_spec.py` | Builds `spec/design_refine_strong.json` and `spec/refine/*.cif` from the parents list |
| `make_designs_table.py`, `validate_and_log.sh` | Base-validation table, ledger duplicate check, `check_denovo.py`, ledger logging |
| `epitope_check.py` | Per-design geometry against 1N8Z (**needs `gemmi`, which the VM lacks: run on a laptop**) |
| `score_designs.py` | Joins everything, ranks epitope-first, picks a diverse shortlist |
| `results/` | `batch1-scored.csv` (all 104, with sequences), `batch1-shortlist-v2.csv` (epitope-first), `batch1-shortlist.csv` (the original 12), `batch1-refine-parents.csv`, `batch1-epitope.csv`, `batch1-shortlist-refold-boltz2.csv` |
| `compare_5my6.*`, `epitope_patch.*` | Evidence for the structure and epitope decisions |
| `refold.sh`, `refold_each.py`, `select_for_refold.py` | Refold tooling, Boltz removed, **unused under the RFD4-only rule** |
| On the VM only | `runs/<run>/design/...` (generated CIFs, git-ignored), `reports/` (tables, checker reports, NanoMelt output), and the round-2 inputs built there from `runs/`: `spec/design_refine_strong.json` and `spec/refine/*.cif` (re-create with `make_refine_spec.py`) |

## 5. How to continue

**Check the queue** (read only):
```bash
ssh team-2 'cd ~/hackathon-team-3/LEO/tasks/denovo-vhh-epitope; tail -5 runs/round2.log; docker ps --format "{{.Names}} | {{.RunningFor}}"'
```
**When a run finishes** (run name, for example `vhh-refine_strong-n16-s4` or `vhh-fixedfw_long-n100-s3`):
1. On the VM: `LEDGER=/home/hackathon/leo-work/generation-ledger/ledger.py ./validate_and_log.sh TAG RUN`
2. Fetch the CIFs into a local folder **named exactly like the run**:
   `ssh team-2 'cd ~/hackathon-team-3/LEO/tasks/denovo-vhh-epitope/runs/RUN/design/generation/inference_samples && find . -name "*_generated.cif" -print0 | xargs -0 tar -cf -' | tar -xf - -C RUN`
3. Locally: `python -I -X utf8 epitope_check.py RUN -o epitope_RUN.csv`, then `score_designs.py` with the table, metrics,
   checker report and epitope CSV from step 1 and 2 (`--help` shows the flags).
4. The number to look at: the **share of designs reaching the strong-epitope bar**, against batch 1's 8%.

**If RFD4 rejects the multi-parent refine spec** (untested): run one parent per job instead,
`--input spec/refine/<parent>.cif --partial-trajectory-fraction 0.4 --num-replicates 16` with the same image and flags as `run_design.sh`.
**Unverified:** that refine children keep the chain layout `epitope_check.py` expects (HER2 chains C and D, VHH chain H). Check the first child.

## 6. What we have learned (batch 1: 100 designs plus 4 smoke)

- **Base validation:** 99 of 104 pass the de novo track. The 5 failures all have 12-residue CDR3s, 70 to 71% identity to trastuzumab, and were assigned to optimisation (never reclassified).
- **CDR3 length and track risk:** my estimate says any design with a CDR3 of 15 residues or fewer has at least 70% framework identity to its caplacizumab scaffold (arithmetic, not `abguard`'s number). Only 17 and longer have a comfortable margin.
- **Epitope signal:** 8 of 99 valid designs reach the team's bar (`epitope_cov` at least 0.6 and `contacts_in_ep` at least 0.7). Mean coverage does not depend on CDR3 length (0.51 to 0.53).
- **RFD4 behaviour:** the output does **not** keep the 1N8Z frame (about 150 Å raw) but HER2 stays rigid (0.05 Å after superposition), so `epitope_check.py` superposes first. Output numbering restarts at 1 per chain: add 488 to chain C and 590 to chain D to get 1N8Z numbering.
- **`freefw` (18 framework positions freed) was not worth scaling:** the model changed only 3 of 18 (83% recovery) and the identity margin was about 2.5 points.
- **Boltz-2 (12 designs, kept, none further):** confident VHH folds (median pLDDT 93), but only `fixedfw_6` (2.2 Å) and `freefw_0` (1.7 Å) reproduced the designed pose. EXP-011 shows Boltz-2 cannot separate binders from non-binders (ipTM at least 0.7 for 94% of known non-binders). Do not treat it as endorsement.
- **NanoMelt** (Nick's run on the 12, mine agrees): Tm 65.6 to 73.5 °C against a 4.1 °C error. Tie-break only.

## 7. Candidates (epitope-first ranking, `results/batch1-shortlist-v2.csv`)

| # | Design | CDRH3 | Hotspots in contact | In-epitope share | Epitope coverage | Track risk | Identity to trastuzumab |
|---|---|---|---|---|---|---|---|
| 1 | `n100-s2` `fixedfw_48` | ARWDPAMRDGSHWDVSGWAY | 4 of 4 | 0.94 | **0.79** | MEDIUM | 0.648 |
| 2 | `n100-s2` `fixedfw_89` | AMASYPYLPTSLWTLREAWSY | 4 of 4 | 0.87 | 0.68 | LOW | 0.630 |
| 3 | `n2-s1` `freefw_0` | ARYNGTGSSTAPSSYDY | 4 of 4 | 0.92 | 0.63 | n/a | **0.683** |
| 4 | `n100-s2` `fixedfw_42` | AAALGSYDPSTVTLRTDYDY | 4 of 4 | 0.92 | 0.63 | MEDIUM | 0.651 |
| 5 | `n100-s2` `fixedfw_18` | AARRSSDSSGYGASYFDY | 4 of 4 | 0.71 | 0.63 | MEDIUM | 0.675 |

Trade-offs worth knowing before ranking the 4 guaranteed slots (this is Leo's decision):
- `fixedfw_48` has the best epitope numbers but the lowest predicted Tm (65.6 °C) and a Boltz-2 ipTM of 0.44 (uninformative, see above).
- `freefw_0` is the closest to the 70% line (0.683 to trastuzumab) and not a fixed-framework design, so my track-risk arithmetic does not cover it.
- `fixedfw_6` is **not** a strong-epitope design (coverage 0.58, in-epitope share 0.65) although Boltz-2 reproduced its pose.
- Sequence IDs are `RUN__EXAMPLE` and example ids repeat across runs, so always carry the run name. Full heavy sequences are in the results CSVs.

## 8. Decisions waiting for Leo

1. Commit the 18 uncommitted files, and push. The results CSVs hold the only copy of the design sequences outside the VM.
2. Merge `origin/nick` (PR #17), or wait for it to land on main.
3. Ranking of the 4 guaranteed slots (section 7).
4. Whether the round-2 queue may use `FORCE=1`: other agents keep chaining `her2-novel` runs, so it may starve.
5. Whether the Protenix v2 refold that EXP-011 recommends is wanted after all (currently excluded by the RFD4-only rule).

## 9. What could be wrong

- **Track assignment is unconfirmed.** The installed `abguard` is not confirmed as the official release, and the local checker measures identity to trastuzumab only. The scaffold is caplacizumab, which the official database probably holds. A design the local checker passes could still be called optimisation.
- **Epitope numbers are geometric.** They measure the designed pose, not a prediction that the pose is real. There is **no binding predictor** for these designs, and no thermostability baseline.
- **The thresholds in `score_designs.py` are mine** (hotspot count, shortlist rules, track-risk arithmetic), except the epitope bar, which is the team's.
- **Refinement is unproven here:** the multi-parent spec has never run, and children may drift off the epitope. Score them like any other design.
- 31 valid designs were excluded as "not eligible" for RFD4 `hotspot_obedience` below 1 or a clash fraction above 0 (four of them only about 0.1% of atoms). That gate was my pre-set choice and is strict.

## 10. Gotchas that cost time

- **ssh and background jobs:** `a && b & c` backgrounds the **whole chain** and the ssh call then hangs. Use `;` and `setsid nohup cmd > log 2>&1 < /dev/null & disown`.
- **`pkill -f pattern` matches the ssh shell's own command line** and kills your own session. Kill by PID.
- **`grep $'\r'` inside `$(...)` in double quotes silently becomes an empty pattern.** Count CRs with `tr -cd '\r' | wc -c`.
- **Windows git has `core.autocrlf=true`:** working files can turn CRLF. Check before `scp`; committed blobs are LF.
- **`python -I` ignores `PYTHON*` variables:** pass `-X utf8` for UTF-8 output on Windows.
- **`refold-prepare` needs a batch to share its conditioned residues** (different CDR lengths do not) and a full design-output tree, including `result.json`. `select_for_refold.py` and `refold_each.py` handle it, but they are unused now.
- **The ledger prints one line per design,** so `validate_and_log.sh` output is long.
- **NanoMelt** needs `~/nick-work/nmenv/bin` on `PATH` for `hmmscan`; run it CPU only (`CUDA_VISIBLE_DEVICES=`).

## 11. Deadline

The VM shuts down at the end of Saturday 14 November. Generated CIFs live only in `runs/` on the VM (git-ignored).
Before then, copy out the CIFs of any design that is shortlisted, and get the results CSVs committed and pushed.
