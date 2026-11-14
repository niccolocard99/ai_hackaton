# AGENTS.md

Starting context for any agent working in this repo. Read this first, then `LEO/README.md` for the full notes.

## What this project is

Team entry for the Oslo AI × Antibodies hackathon (Saturday 14 November 2026). The goal is to design HER2 antibodies that beat trastuzumab, using Apheris Foundry (RFD4-Proteina models) and other tools. Each design is scored against a measured trastuzumab baseline on three properties: binding (SPR), epitope preservation (competitive binding against trastuzumab), and thermostability (direct measurement). Adaptyv expresses designs cell-free and runs the assays; results arrive after the event.

## The two tracks are assigned, not chosen

`abguard` decides each design's track from its sequence. The team does not pick.

- Identity to nearest known antibody **under 70%** → **de novo** track.
- Everything else → **optimisation** track.
- A design is **never moved between tracks**. If it misses its own track's criteria it is rejected, not reclassified.
- A team may submit **any mix** of the two. The tracks are scored independently, with a separate prize pot each.

So "which track are we doing?" is the wrong question. The right one is "what is each design's identity to known antibodies, and does it then pass that track's extra check?"

## Hard constraints (all enforced by `abguard`)

`immunum` first strips tags and linkers to extract the antibody sequence. Then:

**Shared checks, every design**
- Numbering confidence at least 0.1 (is it an antibody?)
- Format is VHH, scFv or Fab
- Total length at most 250 aa (single scFv or VHH, or VH and VL combined)
- No unpaired cysteines

**Identity comparison is type-aware**
- A VHH is compared against known VHHs and heavy chains.
- An scFv or Fab is split into paired chains and compared as paired antibodies.

**De novo track, extra check**
- At least 2 CDRH3 mutations away from trastuzumab

**Optimisation track, extra checks**
- CDRH3 identity to trastuzumab strictly under 70%
- **Strictly more than 2** point mutations from any provided data point

> **Off-by-one warning.** Page 2 of the task description says "at least 2 mutations", but the flowchart and the hard-filter list both say *more than* 2. Use **3 as the minimum** for optimisation designs so nothing is lost to the stricter reading.
>
> **The CDRH3 identity rule binds before either mutation count.** Trastuzumab's CDRH3 is `SRWGGDGFYAMDY` (13 aa), so at equal length a design needs **at least 4 substitutions** to fall under 70% identity (9/13 = 69.2%; 3 substitutions = 76.9%). Indels were not tested. The 3-mutation floor is checked against every CDRH3 in the reference set, not only trastuzumab. Verified 14 Nov with `NICCOLO/scripts/check_optimisation.py`: a 3-substitution test variant scored 77% identity and was also only 2 mutations from a CDRH3 in the FACS table, and failed on both.

**Novelty set**
- A design's CDRH3 must differ from all **525,063** CDRH3s in the HER2 tables (`apheris/hackathon-data/her2-trastuzumab/`). The 524,346 figure elsewhere is the FACS variant count, a different number for a different purpose.

**RFD4-Proteina input rules**
- PDB files are rejected. Use annotated CIF/BCIF, compressed CIF, or a structure directory.
- Pass the design spec file as `input`, not its parent folder. The parent folder silently ignores `contig` and `conditions`.
- Start at `batchSize 1` for the longest input to avoid running out of GPU memory.

**Deadline**
- The VM shuts down at the end of Saturday 14 November. Copy anything worth keeping before then.

## BASE VALIDATION (mandatory, always performed)

Every design must go through base validation **before** it is ranked, refolded, written into a submission list, or reported as a candidate. No exceptions, and re-run it whenever a sequence changes. A design that has not passed is not a candidate.

1. Put the designs in a table with columns `sequence_id, kind, heavy_sequence, light_sequence, design_description` (`kind` is `vhh`, `scfv` or `paired`; a VHH or scFv goes whole in `heavy_sequence`).
2. Run **both** track scripts from `NICCOLO/scripts/` on that table. They wrap `abguard` (numbering, format, length, cysteines, novelty) and add each track's extra rules:

```bash
NICCOLO/scripts/check_denovo.py       designs.csv -o reports/denovo.csv
NICCOLO/scripts/check_optimisation.py designs.csv -o reports/optimisation.csv \
    --provided apheris/hackathon-data/her2-trastuzumab/trastuzumab_her2_spr_all_768.csv
```

3. A design is eligible only if it shows `PASS` in the script for its **assigned** track (`assigned=` in the output). `FAIL` in the other track's script is expected, because a design is never reclassified.
4. Keep the report next to the run notes, and record the verdict with the design.

What the scripts check:
- **Shared (abguard):** numbering confidence at least 0.1, format, length at most 250, no unpaired cysteines, unique ids.
- **De novo:** identity under 70% (abguard `de_novo` class), CDRH3 at least 2 mutations from trastuzumab and every known CDRH3.
- **Optimisation:** identity 70% or more (abguard `novel` class), CDRH3 identity to trastuzumab under 70%, at least 3 mutations from every provided data point (CDRH3 vs the reference, and whole VH+VL vs `--provided`).

Caveats:
- `abguard` and `immunum` are on `PATH` since 14 Nov (`~/.local/bin`, venv `~/.venvs/abguard`, abguard 0.1.0 installed from the wheel in the uv cache; `abguard --help` lists the `--reference-db` novelty gate). The scripts still fall back to the newest uv-cache build with the novelty gate; set `ABGUARD_PY` to override. The older cached build has no novelty gate. This is the build found on the VM; it has not been confirmed to be the organisers' official release.
- abguard measures the 70% identity against trastuzumab only, not all known antibodies, so a design near some other known antibody can still be labelled de novo. Treat the script verdict as necessary, not sufficient, until the organisers' official `abguard` is run.
- The bundled reference is abguard's own CDRH3 index. For the full check, sync `hackathon-data/` and pass `--reference-db` built from it.
- Passing base validation says nothing about binding, epitope or thermostability. It is only the gate.

## GENERATION LEDGER (mandatory, check before every run)

Only 4 designs get tested and the H100 is shared, so the same generation must never be paid for twice. Every run that produces sequences (RFD4-Proteina, a Nebius or JapanFold job, ProteinMPNN, or a scripted or hand-made set of mutants) is logged in the CSV ledger in `LEO/tasks/generation-ledger/`: `runs.csv` (one row per run) and `variants.csv` (one row per sequence per run). `ledger.py` there needs nothing installed. The CSVs are append-only: **never edit or delete a row by hand**. A changed status or verdict is a new row, and the latest wins.

Two runs are the same when method, model image, the **contents** of the input files, seed and settings all match. Two sequences are the same when heavy and light chain match, ignoring case and whitespace.

1. **Before launching, `check-run`.** Exit code 3 means that exact run is already in the ledger: **do not launch it**, reuse the `output` it prints. Change the seed or a setting only when you want a different result. To repeat on purpose (outputs lost, nondeterminism check), add `--setting rerun=2` and say why in `--notes`. A run that only ever `failed` may be retried once its cause is fixed.
2. **Register it with `add-run`** (same flags) before the job starts, so a teammate checking at the same moment sees it. When it ends, `set-status --status done|failed|empty`.
3. **Log what came out** with `add designs.csv --run-id <id>`, using the BASE VALIDATION table. It also lists any sequence an earlier run already produced.
4. **Before refolding, ranking or validating any table of designs, `check designs.csv --new-only fresh.csv`** and carry on with `fresh.csv`. A sequence already in the ledger is not regenerated, refolded or re-validated: reuse the earlier result (the ledger shows its last validation verdict).
5. **After BASE VALIDATION,** add `base_validation` (`PASS` or `FAIL`, for the design's assigned track), `track_assigned` and `cdrh3` columns to the table and run `add` again.

```bash
cd ~/hackathon-team-3/apheris/rfd4             # the folder the model command runs from
L=../../LEO/tasks/generation-ledger/ledger.py
RUN=(--method design --model quay.io/apheris/oslo-hackathon:apheris-rfd4-proteina-v0.3.0-her2-novel-2026-11-09
     --input Data/12-her2/her2_cdrh3.yaml --seed 1 --setting num-replicates=100)

python3 $L check-run "${RUN[@]}"                    # exit 3 = already run, stop
python3 $L add-run   "${RUN[@]}" --run-id leo-20261114-her2novel-s1 --output runs/leo-her2novel-s1
#   ... docker run ... module-run design ...  (the same flags as above)
python3 $L set-status --run-id leo-20261114-her2novel-s1 --status done
python3 $L add   designs.csv --run-id leo-20261114-her2novel-s1
python3 $L check other_designs.csv --new-only fresh.csv   # before refolding or ranking any table
```

How to fill the flags: `--setting name=value` for **every** flag that changes what is generated (`num-replicates`, `sequence-tsr-k`, `recycle-strategy`, and so on), spelled as in the command without the dashes, and never the output path. Pass the spec file as `--input`; the structure file its `input:` line points at is hashed automatically. Give `--model` the full image tag. For ProteinMPNN or a mutant script, use `--method` for what it is and `--model` for tool and version.

Caveats:
- The same seed on a GPU may not reproduce bit for bit. The ledger treats it as the same run anyway, because paying for it twice is the cost it exists to avoid.
- Only runs logged here are known. Nothing before the ledger exists in it except `unknown-20261114-her2-smoke`, the HER2 smoke test, as a run row only; its sequences are on the VM and are not backfilled.
- Log the shortlist, not the screen. The pre-filter in `NICCOLO/eda/` scores millions of enumerated CDRH3s, and the ledger is committed: log that as one run (`--method prefilter`, with its settings) and `add` only the designs that leave it.
- `check` also notes when a new sequence shares its CDRH3 with a logged one. That is information, not a duplicate: the framework may differ on purpose.
- The ledger is committed. Never put keys or credentials in `--setting` or `--notes`.
- Tests: `python3 -m unittest test_ledger` in that folder.

## Only 4 designs are guaranteed to be tested

This governs every ranking decision.

- The team submits a **ranked list of up to 20** designs. Only the **top N** goes to the lab, so rank order is the real output.
- **4 designs are tested whatever happens.**
- Pitch placement adds: 1st +12 (16 total), 2nd +8 (12 total), 3rd +4 (8 total). Each of 6 judges nominates a favourite for +2.
- Ranking after testing counts **S-tier** designs (better on all three properties) first, then A-tier (two), then B-tier (one). Ten designs that each improve one property lose to two that improve all three.

Assume 4 tested. Optimise for a handful of designs that plausibly win on all three properties, not for volume or maximum novelty.

## Pitch day is 40% judgment, 20% reproducible repo

| Weight | Criterion |
|---|---|
| 40% | Judgment — clear reasoning behind the molecules chosen |
| 20% | Execution — a reproducible pipeline in a GitHub repo |
| 20% | User fit — what the design is for, and who needs it |
| 20% | Path to value — a credible route to real use |

Two consequences for agents: **record the reasoning** behind every design that makes the shortlist, and **keep the pipeline runnable end to end** from this repo. Both are directly scored.

## Where things are

| Path | Contents |
|---|---|
| `LEO/README.md` | Full notes: data, docs, setup, rules, submission, tasks, progress log |
| `LEO/.venv/` | Python venv with `pypdf` (git-ignored) |
| `apheris/hackathon-data/` | Synced from S3 (4.7 GB). Never commit it. |
| `apheris/rfd4/` | Synced Foundry docs. Start with `18-vm.md`, then `Tutorial.html`. Run all model commands from here. |
| `apheris/` | **Everything Apheris on this VM.** See "Apheris toolkit" below |
| `NICCOLO/scripts/` | `check_denovo.py`, `check_optimisation.py`: per-track design checks on top of `abguard` |
| `JAHN/pipeline/configs/stage1_frozen.yaml` | **Stage 1 of the pipeline, frozen.** Generate, filter, score and rank in one runnable config. Copy it to experiment; do not edit it |
| `NICCOLO/PIPELINE.md` | **Pipeline status**: what runs, what each stage was trained and tested on, what is rejected and why, and the gaps |
| `NICCOLO/eda/AGENT_GUIDE.md` | **Read before ranking any optimisation design.** The CDRH3 pre-filter: which classifier and ranker to use, how to call them, and the measured reasons not to screen on ipTM |

## Apheris toolkit (`apheris/`)

All organiser-supplied tooling lives in `apheris/`. Paths below are relative to it (`~/hackathon-team-3/apheris/`). Full how-to: `apheris/README_APHERIS.md`.

| Path | What it is |
|---|---|
| `README_APHERIS.md` | Tools, data, model commands, agent workflow. Read it before running anything |
| `rfd4/` | Foundry docs, `Data/` input specs, `Settings/`, `Scripts/run-her2.sh`. **Run every model command from here.** `18-vm.md` first, then `Tutorial.html`; `index.md` lists every page |
| `rfd4/runs/` | Model outputs, one folder per run (`runs/<descriptive-name>/`). Git-ignored. `runs/her2/` is the original smoke test |
| `hackathon-data/` | 4.7 GB from `s3://apheris-oslo-hackathon/data/`. `README.md` explains every file. Git-ignored, read-only |
| `foundry-assistant/` | Offline docs Q&A: `make ask Q="..."` (needs llama.cpp at `/opt/llama.cpp`, models in `models/`) |
| `apheris-gui/` | Streamlit dashboard for the Docker commands (the VM has no Foundry Hub). `./run.sh` serves `localhost:8501`; from a laptop `ssh -L 8501:localhost:8501 <vm>`. Runs one heavy job at a time |
| `challenge/` | The two challenge documents (task description, deck) as `.md`. The PDFs were removed from the repo; use the `.md` |

**How to run models.** No Foundry Hub and no `apheris-foundry` CLI here. `foundry-toolkit` (lists and describes Foundry modules) is installed at `/usr/local/bin`; `foundry-toolkit module list` needs `FOUNDRY_REGISTRY`, `FOUNDRY_USERNAME` and `FOUNDRY_PASSWORD` in the environment, which `~/.bashrc` exports for interactive shells but a one-shot `ssh host cmd` does not. Its `login` is browser-based and not needed. Docker GPU passthrough is verified (`nvidia-smi` runs inside the images). Every model is a Docker image `quay.io/apheris/oslo-hackathon:*`, already pulled, no login needed:

```bash
cd ~/hackathon-team-3/apheris/rfd4
Q=quay.io/apheris/oslo-hackathon
RUN="docker run --rm --gpus all --shm-size 16g --user $(id -u):$(id -g) -e HOME=/tmp -v /etc/passwd:/etc/passwd:ro -v $PWD:/work -w /work"
$RUN $Q:apheris-rfd4-proteina-v0.3.0-her2-novel-2026-11-09 module-run design \
  --input Data/12-her2/her2_cdrh3.yaml --output runs/<name>/design --num-replicates 100 --seed 1
```

Tasks per image: `design`, `refold-prepare`, `refold-score`, `metrics`, `finetune`. Folding images for refolding: Boltz-2 (`apheris-boltz-v1.3.0-boltz2`), OpenFold3, Protenix, Mol-1, ESMFold2. Quick end-to-end check: `bash Scripts/run-her2.sh all` (`N=`, `TOPN=`, `OUT=`, `IMG=` override defaults; set `OUT` so you do not overwrite `runs/her2`).

**Constraints to respect**
- **One H100 80 GB, shared by the whole team.** Run one heavy job at a time, in the background, and poll the log. Check `nvidia-smi` and `docker ps` before starting.
- `docker run` mounts the current directory as `/work`. Run from `apheris/rfd4`, and write outputs under `runs/`. Never write into `hackathon-data/` or the docs.
- Pass the spec file (`.yaml`/`.json`) as `--input`, not its folder. CIF/BCIF only, no PDB.
- Rank on `min_ipae` (lower is better), `iptm`, `binder_plddt`. Compare scores only within one folding model.
- Tools with network need keys from the environment. Never write `AMASS_API_KEY`, the Quay robot token (in `~/.bashrc`) or AWS credentials into any file in this repo.
- Untrusted downloads and archives: unzip into their own new empty directory and read with `python3 -I`. The FACS table is a `.zip`; `zenodo-10831512-trastuzumab-ngs/data.zip` is 4 GB unzipped, so extract outside the repo or into a git-ignored folder.

**`abguard` (the validity gate) and the track checkers**
- `abguard` and `immunum` are on `PATH` (see the caveats under Base validation). The uv cache (`~/.cache/uv/archive-v0/*/bin/abguard`) holds two builds; only the newer has the novelty gate. The scripts find the right build themselves, or set `ABGUARD_PY`. Its reference database is the 525,063 HER2 CDRH3s plus prior art.
- `NICCOLO/scripts/check_denovo.py` and `NICCOLO/scripts/check_optimisation.py` wrap it and add each track's extra rules (see their docstrings). Usage: `NICCOLO/scripts/check_<track>.py designs.csv [-o report.csv] [--strict-exit]`; optimisation also takes `--provided trastuzumab_her2_spr_all_768.csv`.
- Local results are predictions. The real track assignment is `abguard`'s, at submission time.

**MCP tools.** The RCSB PDB MCP (`mcp__rcsb-mcp__*`) is registered for Claude Code: structures, sequences, UniProt (HER2, trastuzumab = 1N8Z).

**Moving or re-creating this setup.** It was created by syncing from S3, not by a script. To rebuild from scratch:

```bash
cd ~/hackathon-team-3/apheris
aws s3 sync s3://apheris-oslo-hackathon/data/ hackathon-data/
aws s3 sync s3://apheris-oslo-hackathon/rfd4-proteina/docs/ rfd4/
```

## Models and weights

Foundry models run as Docker images on the VM. Pick the RFD4-Proteina weight by image:

| Image | Use |
|---|---|
| `release-2026-10-19` | General proteins and miniproteins |
| `antibody-2026-10-30` | Design or redesign CDRs on a fixed framework; the only option for CDRs other than CDRH3 |
| `her2-2026-11-07` | Redesign trastuzumab's CDRH3. `--sequence-tsr-k 2` (balanced) or `1` (exploratory) |
| `her2-novel-2026-11-09` | Mostly new CDRH3 loops; best measured novelty pass rate (86%) and lowest interface error (6.8 Å) |

The HER2 images only redesign CDRH3. Note the tradeoff: `her2-2026-11-07` conservative hits 21% known binders but **0%** pass the novelty filter, so its best-scoring designs are exactly the ones that get rejected. Filter for novelty **before** ranking, never after.

Check designs with `refold-prepare` → fold with a Predict Structure model → `refold-score`. Rank on `min_ipae`, `iptm` and `binder_plddt`; the pass column stays empty for the HER2 exercise. For binders prefer a folding model that **takes a template** (OpenFold3, Protenix, Boltz-2, Mol-1), because ESMFold2 and Boltz-1 have no template support and ESMFold2 folds the target unguided.

## Developability: two of three properties have no built-in proxy

The RFD4-Proteina metrics only speak to binding. Thermostability and epitope preservation do not appear in them, and a design that fails to express cannot be measured at all. There is no hard filter for expression, but Adaptyv reports expression success improved when participants scored candidates with **Aggrescan 3D**, **SoluProt**, **NetSolP**, and **ProteinMPNN with soluble weights**. No thermostability predictor is provided; choosing one is an open decision.

Checked on the VM on 14 Nov: Aggrescan 3D, SoluProt, NetSolP and ProteinMPNN are **not installed** (not on `PATH`; `$HOME` searched to depth 4; a `proteinmpnn` skill folder exists in `~/.claude/skills` but no executable was found and its contents were not inspected), and none appears in the Nebius or JapanFold model lists. Two tools that sound relevant are **not** thermostability predictors: Nebius `stab-ddg` predicts *binding* ddG for point mutants of a complex (`predict-binding-ddg`, an affinity proxy), and JapanFold `opendde` / `opendde-abag` dock protein and antibody-antigen complexes.

## Tools and accounts

- **AMASS:** literature, trial, patent and drug search. One person per team claims the $100 credits. The API key comes from the `AMASS_API_KEY` environment variable. Never write a key into the repo. Budget is roughly 2,000 searches; cache anything re-run. **Not claimed yet:** `AMASS_API_KEY` is unset on the VM (checked 14 Nov).
- **Nebius Scientific AI (set up 14 Nov):** the participant kit is extracted to `~/nebius-kit/participant-03/` (outside the repo) and `install.sh` put the 51-skill bundle plus the `nebius-scientific-ai` client in `~/.claude/skills`. `NEBIUS_CREDS` (exported in `~/.bashrc`) points at the credentials file. Use `python3 ~/.claude/skills/nebius-scientific-ai/forge.py --creds "$NEBIUS_CREDS" models|schema|run|submit|status|wait|fetch`. Tested: the API key authenticates and lists models; one schema (`stab-ddg`) was read. **No Nebius job has been run yet.** Models of interest: `antifold`, `abodybuilder3`, `nanobodybuilder2`, `igbert`, `esm-if1`, `saprot-inverse-folding`, `stab-ddg` (binding ddG, see Developability), `boltzgen`, `rfdiffusion`, `rfd4-proteina`, `alphafold3`, `chai1`.
- **JapanFold (set up 14 Nov):** API at `https://api.japanfold.aiand.com`. Key in `~/.config/japanfold/api_key` (mode 600), exported as `JAPANFOLD_API_KEY` from `~/.bashrc`. Skill at `~/.claude/skills/japanfold/`, `japanfold` CLI 1.1.1 (venv `~/.venvs/japanfold`), and a `japanfold` MCP server registered at user scope (shows Connected). Balance was $1,100.00 before the first job ($1,000 hackathon credit + $100 signup). Models (`GET /v1/models`, no key needed): predict `rf3` (RoseTTAFold3), `openfold3`, `boltz2`, `esmfold2`, `esmfold2-fast`, `mol1`, `openbind`; dock `opendde`, `opendde-abag`; score an existing binder `af2ig`; design `boltzgen`, `rfd3`, `pxdesign`, `rfd4`, with protocols such as `nanobody-anything`, `antibody-anything` and `rfd4-binder`. Verified: an RF3 fold of the trastuzumab Fv (227 residues) succeeded in about 50 s for about $0.0017. Typical p90 cost is roughly $0.001 to $0.02 per target. Never print the key; test it with `[ -n "$JAPANFOLD_API_KEY" ]`.
- **Rebuild commands** for abguard, Nebius and JapanFold are in `README.md` section 4.6.
- **Claude Code:** pre-installed on the VM.

## Rules for agents

- **Read the `.md`, never the PDF.** The challenge documents live as Markdown in `apheris/challenge/`; the PDFs are removed. Always open the `.md`. If you add a new PDF, convert it to `.md` next to it first and use that.
- **The Nebius participant bundle is already extracted and installed** (`~/nebius-kit/participant-03/`, see "Tools and accounts"). Do not unzip it again into the repo. JapanFold is configured separately and is not in that bundle.

- Always run **BASE VALIDATION** (above) on every design before ranking or submitting it.
- Never print or echo a key; refer to `$JAPANFOLD_API_KEY` and `$NEBIUS_CREDS` by name. Never commit secrets or API keys, nor `apheris/hackathon-data/`, `apheris/rfd4/runs/` and `LEO/.venv/`.
- Do not commit or push unless asked. End commit messages with a `Co-Authored-By:` line naming the model you are actually running as.
- Record every design run: input spec, weight image, settings, seed, and resulting scores. Put run notes in the owner's folder, not in shared folders.
- **Check the GENERATION LEDGER (above) before every generation run, and before refolding, ranking or validating a table of designs.** Log each run and its variants in it. A run or sequence it already holds is not regenerated: reuse the earlier output.
- Keep each design's provenance and its **predicted** track (from its identity numbers) with the design. The real track assignment only comes from `abguard`, so treat local labels as predictions and verify them.
- When something is blocked or unclear, record it in the progress log in `LEO/README.md` rather than guessing.
- Prefer verifying a rule against `apheris/challenge/*.md` over trusting a summary, including this one. The documents contradict themselves in at least one place (the mutation count).

## Open items

- Participant field guide did not load when fetched. Open it in a browser and add its contents to `LEO/README.md`.
- The developability write-up linked from the task description ("check the whole write-up here") is not captured. Find it on the hackathon site.
- The JapanFold key was shared in a chat on 14 Nov. Rotate it at `japanfold.aiand.com/account/token`, overwrite `~/.config/japanfold/api_key`, and re-run `claude mcp add` for the `japanfold` server (the MCP entry stores the old key in `~/.claude.json`).
- Confirm with the organisers that the installed `abguard` 0.1.0 is the official release.
- AMASS credits are not claimed.
- AMASS onboarding walkthrough link returned 404.
- No thermostability predictor chosen, and no developability scorer installed (Aggrescan 3D, SoluProt, NetSolP, ProteinMPNN).
