# Leo: Hackathon Info

Owner: Leo

## 1. Overview

- **Data and docs bucket:** `s3://apheris-oslo-hackathon/`
- **Main topic:** nanobody / antibody design with RFD4-Proteina on Apheris Foundry
- **Key constraint:** designed CDRH3 sequences must differ from all 525,063 CDRH3s in the HER2 tables (novelty filter)
- **The track is assigned per design, not chosen by the team.** `abguard` decides: identity to the nearest known antibody under 70% puts a design in the de novo track, everything else goes to optimisation. A design is never moved between tracks; if it misses its own track's criteria it is rejected, not reclassified. Teams may submit any mix of the two tracks.

### Numbers that look inconsistent but are not

- **524,346** — CDRH3 variants in the FACS screen (Chinery et al.), the headline dataset.
- **525,063** — the total across the HER2 tables, which is the set the novelty filter compares against.
- Use 525,063 as the novelty set. Use 524,346 when describing the FACS data itself.

## 2. Data

Download with:

```bash
aws s3 sync s3://apheris-oslo-hackathon/data/ hackathon-data/
```

| Folder | Size | Contents |
|---|---|---|
| `her2-trastuzumab/` | - | FACS binding classes (524k CDRH3s) and SPR affinities (768 variants). Source for the novelty filter. |
| `zenodo-18390239-ab-ag-complexes/` | - | Scores for 1.8M predicted nanobody-antigen complexes |
| `sigma2-ab-ag-champloo/` | - | Per-sample confidence (ipTM, pLDDT, PAE) for the same predictions |
| `sigma2-ab-ag-champloo/spr-supplementary/` | - | Measured SPR for a few of those pairs |
| `zenodo-10831512-trastuzumab-ngs/` | - | Ready-made train/val/test splits of the FACS data |

Total: about 5 GB. Each folder's `README.md` explains its files.

**Status: synced** to `~/hackathon-team-3/apheris/hackathon-data/` (4.7 GB, git-ignored). Never commit it.

## 3. Docs and Instructions

Download with:

```bash
aws s3 sync s3://apheris-oslo-hackathon/rfd4-proteina/docs/ ~/hackathon-team-3/apheris/rfd4
```

Then open `~/hackathon-team-3/apheris/rfd4/Tutorial.html`.

**Status: synced** to `~/hackathon-team-3/apheris/rfd4/`.

| File | What it covers | Where to start |
|---|---|---|
| `foundry.md` | Running RFD4-Proteina through Apheris Foundry (UI, CLI, agents), from design to evaluation | - |
| `18-vm.md` | "Run on your hackathon VM": the VM has no Foundry Hub, but every Foundry model is already on it as a Docker image | **Start here** |

## 4. Setup and Connections

### 4.0 Tool Checklist

| Tool | Status | Link / how to connect | Notes |
|---|---|---|---|
| Apheris Foundry (RFD4-Proteina models) | Pre-installed on VM | `18-vm.md` in the docs | Models run as Docker images on the VM |
| Claude Code | Pre-installed on VM | - | Credits provided by Anthropic |
| abguard (validity CLI) | To be provided by organisers | Check the Slack folder | Runs the hard checks and novelty filter |
| Nebius Tools / Scientific AI | Available | `sciai.nebius.cloud/models` (Workbench > Apps) | No-code; can be connected to Claude as a skill. Setup files in the Slack folder |
| JapanFold | Available | Link not in the deck. Check the Slack folder. | RoseTTAFold3, OpenDDE, AF2-IG. Can be connected as a skill to Claude Code |
| AMASS | Available | See 4.1 | Literature, trial, patent and drug search with source links |
| Adaptyv antibody database | Reference for de novo check | Named in the challenge deck | Used for the under-70%-identity check |

### 4.1 AMASS

- **Event page:** https://amass.tech/oslo-ai-antibodies-hackathon
- **Onboarding walkthrough:** https://amass.tech/onboarding/oslo-ai-antibodies-hackathon (returned 404 when checked; try again or ask the organisers)
- **Connect to Claude (guide):** https://amass.tech/how-to-do-scientific-work-in-claude
  - Add Amass to Claude once and sign in once. The guide lists four steps, but only three are shown on the page.
  - Ask in plain language about a condition, target, drug, trial, mechanism, or company. Sponsor names and NCT numbers also work.
  - Open the linked records behind each answer and check the cited passages.
  - Amass queries six "Cores" during the conversation: 43M+ papers, 1.2M+ trials, 170M+ patents, 7K+ regulatory authorizations, 22K+ molecules, and 43K+ genes.
  - Filters: date, journal quality, retraction status, phase, and sponsor.
  - Requires Claude Pro, Max, Team, or Enterprise.
  - Data is not used for training. Enterprise hosting is EU-based. Certifications: ISO 27001 and SOC 2; NIS2 is in progress.
  - Install guide: https://amass.tech/mcp. Free trial: https://amass.tech/try.
- **Use the Amass app (guide):** https://amass.tech/how-to-use-amass (app: https://app.amass.tech; demo: https://amass.tech/demo)
  - Ask a question, then read the cited sources and filter by data layer.
  - Timelines view: trial landscape, approval timeline, and literature timeline.
  - Steer the research with a suggested next step or your own follow-up.
  - Attach your own PDF, Word, PowerPoint or text file, and ask questions about it.
  - Build a report: describe it, edit the plan, then read the cited sections. Tables are under Exhibits and export as CSV or Markdown. Reports export to PDF.
  - Organise work with pinned tasks, search, saved reports, and projects. Share a read-only task link that needs no account.
  - The guide's 13 videos are not in the text, so their steps are not captured here.
- **Help:** Amass Discord at https://discord.gg/sEGaBHMhWa (Amass joins remotely)
- **Platform sign-up:** https://platform.amass.tech/api/auth/sign-up (use the email registered for the hackathon; a free account is enough to start)
- **Claim team credits:** https://platform.amass.tech/event/ai-x-antibodies-hackathon-oslo. **One person per team claims the $100.**
- **Connect to Claude (no code):** add the Amass connector in Claude and sign in once. Claude needs a paid plan. The organisers provide Claude Max credits through Anthropic. Setup guide: https://amass.tech/how-to-do-scientific-work-in-claude
- **Coding path (API skill):**
  - Install: `npx skills add amass-technologies/public-skills --skill amass-api`
  - Skill file: https://github.com/amass-technologies/public-skills/blob/main/skills/amass-api/SKILL.md
  - API key: create at https://platform.amass.tech, then `export AMASS_API_KEY=...` (never commit the key)
- **Budget and costs:**
  - About $100 per team, roughly 2,000 searches or 10,000 record fetches.
  - Search: $0.05 per 20 results. Get-by-id and lookup: $0.01 each.
  - Cache results you re-run.
- **Amass track:** the most creative use of the data wins extra credits. Show which records your tool read and what it did with them.
- **Databases ("cores"):** six, such as BioMedCore (papers) and TrialCore (clinical trials).

### 4.2 Logistics

- **Date:** Saturday 14 November 2026, 09:00 to 20:00. Doors open at 08:30.
- **Venue:** BI Startup Village, BI Norwegian Business School, Nydalsveien 37, 0484 Oslo.
- **Rules:** https://oslo-ai-antibodies-hackathon.vercel.app/
- **Machines shut down** at the end of Saturday 14 November. Copy your files before then.

### 4.3 Still to Get

- **Participant field guide** (https://oslo-ai-antibodies-hackathon.vercel.app/participant_guide): the page came back empty when fetched, so its contents are not captured here. Open it in a browser.
- **JapanFold** link and login.
- **abguard** CLI and the **Nebius** setup files, both in the Slack folder.

### 4.4 Local Tools on the VM

**`~/hackathon-team-3/apheris/foundry-assistant`:** an offline terminal assistant that answers Apheris Foundry questions from the docs. It runs two small local models (llama.cpp), with no network and no tokens. It doesn't need the Foundry Hub, so it can be used before anything else is set up.

- **Status:** both model files are present in `models/` (about 0.5 GB each). `llama-server` is at `/opt/llama.cpp`. Not yet started.
- **Run it:**

```bash
cd ~/hackathon-team-3/apheris/foundry-assistant
make install
make models  LLAMA_BIN=/opt/llama.cpp
make up      LLAMA_BIN=/opt/llama.cpp     # add NGL=0 for CPU
make ask Q="how do I install the foundry toolkit and log in?"
make tui                                  # interactive: /read <id>, /quit
make down
```

- **Knowledge base:** `kb/foundry-v1/kb/`, covering getting started, the CLI (commands, exit codes, troubleshooting), and the workflows (RFD4-Proteina design and finetune, affinity, AQ, ADMET, and structure preprocessing).
- **Useful pages:**
  - `getting-started-gs-install-login.md` and `getting-started-gs-run-a-workflow.md`: install, log in, and run a workflow.
  - `cli-troubleshooting.md`: auth errors (exit codes 3 and 8), SSH login port forwarding.
  - `workflows-rfd4-proteina-design.md`: the design task parameters.
- **Gotchas from the RFD4-Proteina design page:**
  - **PDB inputs are rejected.** Use annotated CIF/BCIF, compressed CIF, or a structure directory.
  - Pass the design spec file itself as `input`, not its parent folder. The parent folder silently ignores the spec's `contig` and `conditions`.
  - Use `batchSize 1` for the longest input to avoid running out of GPU memory.
  - A custom checkpoint must match the image's architecture and paired autoencoder.

## 5. Reference Documents

### 5a. Apheris Challenge Deck (`docs/2026-11_Oslo_hackathon_Apheris_challenge.pdf`, 9 pages)

- **Task:** build better HER2 binders with Apheris Foundry
- **Tracks** (the deck's prize figures are marked draft; see section 6 for the task description's numbers, which are newer):
  - *Optimisation:* improve trastuzumab. Same format as trastuzumab. Requires CDRH3 identity to trastuzumab strictly under 70% and strictly more than 2 point mutations from any provided sequence.
  - *De novo:* design a new HER2 binder. VHH, scFv or Fab from scratch. Requires under 70% identity to any known antibody, checked against Adaptyv's antibody database.
- **Entry checks:** run automatically at submission. Failing designs are not sent to the lab. No unpaired cysteines; up to about 250 aa in scFv format.
- **Selection:** each team ranks its own designs. Adaptyv makes and tests 100 designs: 64 guaranteed (4 per team × 16 teams) and 36 extra.
- **Scoring:** binding to HER2 by SPR, epitope (still competes with trastuzumab), and thermostability.
- **Winning:** the team with the most S-tier designs (better on all three) wins. Ties go to A-tier, then B-tier. Lab results arrive after the event.
- **Tools mentioned:** RFD4-Proteina (Foundry), Protenix, Boltz, OpenFold3, Mol-1, ESMFold2, and the abguard CLI.

### 5b. Task Description (`docs/ai-antibodies-hackathon-apheris.pdf`, 9 pages)

- **Task:** improve an existing antibody (trastuzumab) with Apheris Foundry, Nebius and JapanFold, using multi-modal data.
- **Two modes** (assigned by `abguard`, not chosen):
  - *Optimisation:* point mutations to trastuzumab to improve affinity and thermostability.
  - *De novo:* a new antibody that keeps the original epitope.
- **Judging (per design, against the trastuzumab baseline):**
  1. Affinity: comparable or better binding to HER2 by SPR.
  2. Epitope preservation: judged by competitive binding assay.
  3. Thermostability.
  - Tiers: S = 3 of 3, A = 2 of 3, B = 1 of 3, unsuccessful = 0 of 3.
- **Formats:** VHH, Fab or scFv.
- **Hard checks, exactly as `abguard` applies them:**
  1. `immunum` extracts the antibody sequence, stripping tags and linkers.
  2. *Shared checks, every design:* numbering confidence at least 0.1 (is it an antibody?); format is VHH, scFv or Fab; total length at most 250 aa (a single scFv or VHH, or VH and VL combined); no unpaired cysteines.
  3. *Type-aware identity comparison:* a VHH is compared against known VHHs and heavy chains; an scFv or Fab is split into paired chains and compared as paired antibodies.
  4. *Track assignment:* identity to nearest known antibody under 70% → de novo; otherwise → optimisation.
  5. *De novo only:* at least 2 CDRH3 mutations away from trastuzumab.
  6. *Optimisation only:* CDRH3 identity to trastuzumab under 70%, **and strictly more than 2 point mutations** from any provided data point.
- **Watch the mutation count.** Page 2 of the PDF says "at least 2 mutations", but the flowchart on page 3 and the hard-filter list on page 6 both say *more than* 2 ("strictly more than 2 point mutations"). Treat 3 as the minimum for the optimisation track, so a design with exactly 2 mutations is not lost to the stricter reading.
- **Supporting data:**
  - General VHH-antigen data: 454 true-negative pairs (21 with AF3 ipTM at or above 0.7), 2 non-cognate labels, 15 cognate affinities (sub-micromolar), and 1.8M predicted complexes (Smorodina et al.).
  - Trastuzumab: 524,346 CDRH3 variants vs HER2 by FACS (Chinery et al.), and about 700 optimisation designs, 190 with SPR affinity.
- **Compute:** one H100 GPU, 16 vCPUs, 200 GB RAM and 500 GB disk per team (Ubuntu 24.04). **Machines shut down at the end of Saturday 14 November. Copy anything you want to keep before then.**
- **Platforms:** Apheris Foundry, Nebius Tools (no-code, also usable as a Claude skill), JapanFold (RoseTTAFold3, OpenDDE, AF2-IG), and AMASS (literature and patent search, $100 credits per team, claimed by one person).

### 5c. RFD4-Proteina Workflow (from the task description)

- **Tasks in each Docker image:** `design`, `refold-prepare`, `refold-score`, `finetune`.
- **Weights:**

| Image | Use it for |
|---|---|
| `release-2026-10-19` | General proteins and miniproteins |
| `antibody-2026-10-30` | Design or redesign CDRs on a fixed framework. Use this to change CDRs other than CDRH3. |
| `her2-2026-11-07` | Redesign trastuzumab's CDRH3. Set `--sequence-tsr-k` to 2 (balanced) or 1 (exploratory). |
| `her2-novel-2026-11-09` | Mostly new CDRH3 loops with the most confident predicted fit to HER2 |

- **Novelty results on the CDRH3 task:**

| Weight and setting | Known binders | Pass novelty filter | Interface error (Å) |
|---|---|---|---|
| `antibody-2026-10-30` (not fine-tuned) | 0% | 100% | 8.2 |
| `her2-2026-11-07`, conservative | 21% | 0% | 9.4 |
| `her2-2026-11-07`, balanced | 21% | 6% | 9.1 |
| `her2-2026-11-07`, exploratory | 0.3% | 74% | 9.1 |
| `her2-novel-2026-11-09` | 4% | 86% | 6.8 |

- **Checking designs:** run `refold-prepare`, fold the inputs with any Predict Structure model, then run `refold-score`. For the HER2 exercise, rank on `min_ipae`, `iptm` and `binder_plddt`. Filter for novelty before ranking, because known binders fail the filter.
- **Which folding model to refold with.** For binders, prefer a model that takes a template: it gives a fairer check, because ESMFold2 folds the target unguided.

| Model | Weights | Template support |
|---|---|---|
| OpenFold3 | OpenBind, preview 2 | Yes |
| Protenix | v2, v1 | Yes |
| Boltz-2 | Boltz-2 | Yes |
| Boltz-1 | Boltz-1 | No |
| Mol-1 | base | Yes |
| ESMFold2 | release | No |

- **`refold-score` pass criteria** are interface PAE, binder pLDDT and RMSD to the design. For the HER2 exercise the pass column stays empty, so rank on the three metrics above instead of relying on a pass flag.

## 6. Submission, Pitch and Prizes

This is the part that decides how much of our work actually gets tested. It was missing from the earlier notes.

### 6.1 Pitch day judging (decides how many designs get tested)

Pitches are judged on four criteria:

| Weight | Criterion |
|---|---|
| 40% | **Judgment** — clear reasoning behind the molecules you chose |
| 20% | **Execution** — a reproducible pipeline in a GitHub repo |
| 20% | **User fit** — what your design is for, and who needs it |
| 20% | **Path to value** — a credible route to real use, for the molecule or the method |

Judgment is the largest single weight, so the reasoning for why we picked each design matters more than the volume of designs. Execution is scored on a reproducible GitHub repo, which is what this repo is for: keep the pipeline runnable end to end.

### 6.2 How many designs get tested

- Each team submits a **ranked list of up to 20 designs** (16 assuming a win, plus 4 spare in case of duplication between teams). Only the **top N** of that list goes to the lab, so ranking order is the decision that matters.
- **Every team gets 4 designs tested, whatever happens.**
- Pitch placement adds more: 1st **+12** (16 total), 2nd **+8** (12 total), 3rd **+4** (8 total).
- Each of the **6 judges** nominates a favourite team for **+2** designs.
- Teams may submit **any proportion** of de novo and optimisation designs. The two tracks are evaluated independently.

Practical consequence: assume 4 tested unless we place. Our top 4 in rank order should be our highest-conviction designs, not our most novel or most numerous.

### 6.3 Experimental evaluation

Adaptyv does **cell-free expression** of each design, then measures:

| Property | How it is measured |
|---|---|
| Binding to HER2 | SPR |
| Epitope preservation | Competitive binding against trastuzumab |
| Thermostability | Direct measurement |

Trastuzumab itself is measured too, and is the baseline every design is compared against.

### 6.4 Post-evaluation ranking

- **S-tier:** improves on all three (affinity, epitope preserved, thermostability)
- **A-tier:** improves on two of three
- **B-tier:** improves on one of three

Teams are ranked by most S-tier designs, then most A-tier, and so on. Ten designs that each improve one property lose to two designs that improve all three. So a few designs that win on all three beat a large pile of partial wins.

### 6.5 Prizes

- Total fund **€3000**, split equally between the de novo and optimisation tracks (€1500 each).
- Per track: **€800** (1st), **€500** (2nd), **€200** (3rd).
- **Best Pitch:** an extra $500 on the day.
- **AMASS:** picks the team that best uses their platform and awards one month free.
- Note: the challenge deck quotes $720/$360/$270 and marks the figures as draft. The task description's euro figures above are the newer ones.

## 7. Developability and Expression

Two of the three judged properties (thermostability and epitope preservation) have no proxy in the RFD4-Proteina scoring metrics, and a design that does not express cannot be measured at all. There is **no hard filter** for expression, but Adaptyv reports that expression success improved once participants scored candidates with:

- **Aggrescan 3D** — aggregation propensity from structure
- **SoluProt** — solubility prediction
- **NetSolP** — solubility prediction
- **ProteinMPNN with soluble weights** — sequence design biased toward solubility

The task description links a fuller write-up ("Please check the whole write-up here"), which is not captured here. Find that link on the hackathon site.

Open gap: nothing in the provided tooling predicts **thermostability** directly. Worth deciding early whether to add a thermostability predictor, since it is one of the three S-tier conditions.

## 8. Tasks and Ideas

Goal: do interesting things with the data.

- [ ] Data analysis on the HER2 binding classes and SPR affinities
- [ ] Binding prediction
- [ ] Bring in other models that seem promising
- [ ] Design CDRH3s with RFD4-Proteina and pass the novelty filter against the 525,063 HER2 CDRH3s

## 9. Progress Log

| Date | Item | Status |
|---|---|---|
| 2026-11-14 | Hackathon info added to this folder | Done |
| 2026-11-14 | Two hackathon PDFs added to `docs/` | Done |
| 2026-11-14 | Setup and connections (AMASS, tool checklist) organised | Done |
| 2026-11-14 | Participant guide: fetch returned empty | Open |
| 2026-11-14 | Step 1 setup: connect Foundry, AMASS, JapanFold, Nebius, abguard | In progress |
| 2026-11-14 | Local foundry-assistant reviewed and documented (not yet started) | Done |
| 2026-11-14 | Full PDF re-read; corrected track assignment, mutation count, hard checks | Done |
| 2026-11-14 | Added submission/pitch/prize rules and developability tools (were missing) | Done |
| 2026-11-14 | Sync `hackathon-data/` and `~/hackathon-team-3/apheris/rfd4` from S3 | **Done** |
| 2026-11-14 | Pick a thermostability predictor | Open |
| 2026-11-14 | Generation ledger added (`LEO/tasks/generation-ledger/`); rule added to `AGENTS.md` | Done |
| 2026-11-14 | RFD4 via JapanFold + base validation (`LEO/tasks/rfd4-japanfold/`): 9 of 18 candidates eligible on optimisation, 0 on de novo | Done |
| 2026-11-14 | Design loop with NanoMelt thermostability scoring (`LEO/tasks/design-loop/`) | In progress |
| 2026-11-14 | Fixed-framework VHH batch `vhh-fixedfw-n100-s2` (100 designs) on the H100; only 4 smoke designs validated so far | In progress |
| 2026-11-14 | Results of the live JapanFold loop and the 100-design VHH batch not yet scored or logged | Open |
