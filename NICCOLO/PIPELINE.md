# Design pipeline: current status

Status of the HER2 antibody pipeline as built in `NICCOLO/`, as of **14 November 2026**. It records what
runs, what it was trained on, what it was tested on, and what is still missing. Written to be read before
anyone trusts a number from it.

Scope note: everything here is the **optimisation track** — 10-residue CDRH3 variants on trastuzumab's
VH/VL framework. No part of it is validated for de novo designs (<70% identity to known antibodies).

---

## 1. Status at a glance

| Stage | What it does | Status | Validated on |
|---|---|---|---|
| 0. Hard rules | novelty vs 525,063 CDRH3s; ≥4 substitutions | **ready** | rule-derived, verified against `check_optimisation.py` |
| 1. Binder classifier | P(measurable KD) | **ready** | 620 measured SPR outcomes, grouped CV |
| 2. Affinity ranker | predicted log10 KD | **ready, with a hard limit** | 179 measured KDs, grouped CV + bootstrap |
| 3. BASE VALIDATION | `abguard` + per-track checks | **ready** (not mine) | `NICCOLO/scripts/` |
| 4. Thermostability | NanoMelt Tm prediction | **tool works, no actionable lever** | NanoMelt's 639 measured nanobody Tm |
| 5. Epitope preservation | — | **not built** | nothing |
| 6. Refold for provenance | Boltz-2 / Protenix structure | **ready, not a filter** | 240 folded structures |

**One-line summary:** stages 0–3 screen sequences fast and are honestly validated; stage 4 says CDRH3
cannot move thermostability; stage 5 does not exist; stage 6 earns its place only as provenance.

---

## 2. The pipeline, stage by stage

### Stage 0 — hard rules (free)

- **Novelty:** a design's CDRH3 must differ from all **525,063** CDRH3s in the HER2 tables.
- **Mutation floor: ≥4 substitutions.** The 70% CDRH3 identity rule binds before the ">2 mutations" rule:
  trastuzumab's CDRH3 is `SRWGGDGFYAMDY` (13 aa), so 3 substitutions give 76.9% identity and fail, while 4
  give 69.2%. Indels untested.
- Cheapest stage, so it runs first.

### Stage 1 — binder classifier (`eda/prefilter.py`)

Logistic regression, L2, `C=0.3`, on one-hot CDRH3 (10 positions × 20 residues = 200 features).

| | |
|---|---|
| Trained on | 620 SPR outcomes: 179 `quantified_binding` vs 441 `no_reportable_kd` |
| Grouped-CV AUROC | **0.927** overall, **0.913** inside the ≥4-substitution band |
| Distance extrapolation | **0.904** training on ≤4 substitutions and testing beyond |
| Baselines beaten | mutation count 0.800, FACS score 0.696 |
| Throughput | 1.5M sequences/s, one laptop process, no GPU |

This is the **robust** component. It degrades gracefully as test sequences get more distant.

### Stage 2 — affinity ranker (`eda/prefilter.py`)

Ridge regression, `alpha=3.0`, on the same 200 one-hot features plus mutation count (201).

| | |
|---|---|
| Trained on | the 179 variants with a measured KD; target log10 KD in nM |
| Grouped-CV Spearman | **+0.687** overall, **+0.629** in band |
| Selection quality | top 20 has median KD **6.5 nM** in band (1.4 nM unrestricted) against **86 nM** overall |

**Its hard limit: it interpolates.** With training-set size fixed, varying only train/test similarity:

| Separation | Classifier AUROC | Ranker Spearman |
|---|---|---|
| random CV | 0.933 | +0.730 |
| Hamming ≤2 | 0.930 | +0.671 |
| Hamming ≤4 | 0.902 | +0.433 |
| Hamming ≤5 | 0.845 | **−0.057** |

Beyond ~4 substitutions' separation from its training set the ranker has **no demonstrated power**. Check a
design's minimum Hamming distance to the 179 measured variants before believing its `pred_log10_kd`; at ≥5,
treat the output as uninformative rather than as weak signal.

Note also that the legal band costs affinity before any model is involved: the tightest binders in this
library are the ones closest to trastuzumab, which the identity rule forbids. Top-20 median goes from
1.4 nM unrestricted to 6.5 nM in band.

### Stage 4 — thermostability (`thermo/`)

- **NanoMelt** is the only usable predictor found. Neither Nebius nor Foundry ships a Tm model: Nebius'
  `stab-ddg` is StaB-ddG (*binding* ddG), and Foundry's AQ workflows are small-molecule tools.
- Installed on the VM, **verified to reproduce its own reference output exactly**, ~45 sequences/s.
- Published performance: 639 measured nanobody Tm, held-out Pearson 0.853, MAE 4.1 °C.
- **CDRH3 cannot move thermostability.** 41 variants of the identical 120-aa VH differing only in CDRH3:
  predicted Tm spread 4.63 °C, **0 of 40 beyond the 4.1 °C MAE**, Spearman(predicted Tm, measured KD) +0.001.
- **Control proving the model is not blind:** single point mutations, framework (n=57, mean −0.76 °C, 3
  beyond MAE) vs CDRH3 (n=20, mean −0.09 °C, 0 beyond), Mann-Whitney **p = 0.002**.
- **Consequence: spend CDRH3 on affinity and epitope only.** Framework mutations can move Tm on paper, but
  no single one clears the error bar and the combinations that do are dominated by charge-adding
  substitutions a VHH-trained model is biased toward. Not submission-ready.

### Stage 6 — refolding (provenance only)

Boltz-2 or Protenix v2 against HER2 domain IV with the 1N8Z template. Keep it for the designs we submit,
because provenance is scored. **Do not use it to screen or rank** — see §4.

---

## 3. Data: what trained on what, and what was tested on what

### Used to fit the pipeline

| Dataset | Size | Role |
|---|---|---|
| `her2-trastuzumab/trastuzumab_her2_spr_all_768.csv` | 768 rows → **732** 10-mer variants | the only training data for stages 1–2 |
| ↳ of those | 179 `quantified_binding`, 441 `no_reportable_kd` | classifier: 620; ranker: 179 |
| ↳ excluded | 8 trastuzumab reference rows, 8 rows of an unrelated antibody (`ARDRIVDDYYYYGMDV`), 20 9-mers, 99 `inconclusive`, 13 `not_reported` | — |

Trastuzumab's own 8 rows (measured KD 0.029–0.178 nM) are held out of training and used as a reference
point. All 732 variants share **one** VH framework and **one** VL, so only the 10 CDRH3 residues vary —
which is why a one-hot encoding spans the design space exactly.

### Used to test, not to fit

| Dataset | Size | What it was used for | Outcome |
|---|---|---|---|
| `tz_her2_facs_3classes_524346.csv` | 524,346 unique CDRH3s (high 172,149 / medium 190,011 / low 162,186) | novelty reference; tested as a feature and as a model | **rejected as a feature** |
| Mehta clonotype splits `0.7_0.15_0.15` | 367,044 / 78,651 / 78,651, 32.83% positive, 13,829 / 251 / 706 clusters | leakage-controlled test of a FACS-trained model | exposed the in-library illusion |
| 37 CDRH3s shared SPR ↔ FACS | 37 | external check of FACS labels against measured KD | FACS classes do not track KD |
| `shanehsazzadeh_et_al/test.csv` | 198 CDRH3s, all binders, with affinity | external affinity check | Spearman +0.307 |
| `iptm_confidence_scores.zip` | 91 × 91 matrices × 3 fold models | published real-vs-shuffled benchmark | AF3 0.75, Boltz-2 0.64, Chai-1 0.57 |
| `topmodel_structure_accuracy.zip` | ~720k scores per model | structure accuracy, real vs shuffled | weaker than ipTM, near chance |
| `spr-supplementary/` | 29-row Table S1, 12-row cross-reactivity | sanity check only | too small to train on |

### Structures folded on the VM (all by us, this pipeline)

| Run | Model | n | Wall time | Per structure |
|---|---|---|---|---|
| Calibration, mixed outcomes | Boltz-2 | 60 (30 with KD, 30 without) | 475 s | **7.9 s** |
| Same 60, second model | Protenix v2 | 60 | 772 s | 12.8 s |
| All measured binders | Protenix v2 | 180 | 2,050 s | 11.3 s |

Input: chain A = VL (107 aa), B = VH (120 aa), C = HER2 489–580 (92 aa), D = HER2 591–607 (17 aa); 336
tokens; 1N8Z template; 1 diffusion sample, seed 42. HER2 is split at its unresolved 581–590 loop.
**All 240 structures share 326 of 336 tokens** — only the CDRH3 varies.

### Not available, and it matters

- **No thermostability labels** anywhere in the hackathon data.
- **No epitope-preservation labels** anywhere in the hackathon data.
- The **454 VHH–antigen true negatives** and **15 cognate affinities** described in the organisers' summary
  are not in the delivered files; only the 29-row and 12-row supplementary tables exist. Asked about, unresolved.

---

## 4. Tested and rejected (with the numbers)

| Candidate | Verdict | Evidence |
|---|---|---|
| **FACS score as a ranker feature** | rejected | ranker +0.671 with it, +0.673 without, +0.012 alone. Dropped: removes an 11 MB dependency and triples throughput |
| **FACS-trained classifier as the screen** | rejected | AUROC 0.995 on its own grouped test set but **0.757** externally, below mutation count's 0.802 |
| **FACS classes to rank affinity** | rejected | among the 37 shared CDRH3s, median KD 7.2 nM `high`, 4.2 nM `medium`, 15.4 nM `low` — not monotonic |
| **ipTM / fold confidence as a gate** | rejected | 0.502 AUROC; Boltz-2 ipTM ≥0.7 passes **100%** of true non-binders |
| **ipTM / fold confidence as a ranker** | rejected | at n=180, ipTM vs log10 KD r = **−0.065**, CI [−0.190, +0.082]. Adding it to the ranker *harms*: +0.683 → +0.666, CI [−0.031, −0.002] |
| **`min_ipae` (interface PAE)** | rejected | 0.497 AUROC despite real dynamic range (spread 3.0) |
| **T0 sequence-only ddG features** | rejected | all four cluster-bootstrap CIs span zero |
| **Cheap features for Tm** | rejected | Pearson +0.40, MAE 7.6 °C against 8.3 °C for predicting the mean |
| **Gradient boosting / random forests** | rejected | lost to linear on both tasks (0.915 vs 0.927; +0.628 vs +0.667) |

---

## 5. How the validation is done, and why to trust it

Random cross-validation is meaningless here: variants differ by 1–2 residues, so neighbours straddle a
random split. Guards:

1. **Grouped CV.** Single-linkage clusters at Hamming ≤2 held out together. Verified: 0 identical sequences
   shared between folds, closest train–test pair Hamming 3.
2. **Permutation controls.** Classifier 0.486 ± 0.023, ranker +0.045 ± 0.028 on shuffled labels.
3. **Nested CV.** Hyperparameters chosen inside each fold: 0.925 / +0.649 against 0.930 / +0.671 fixed, so
   model-selection optimism is ~0.005 AUROC and ~0.02 Spearman.
4. **Cluster bootstrap** for every accept/reject decision, because repeat-to-repeat CV spread measures only
   fold-assignment noise and understates sampling uncertainty by 2–3x.
5. **Band reporting.** Every figure also given inside the ≥4-substitution legal band.

Two traps found the hard way, both of which first looked like wins: standardising the binary one-hot block
handicaps the baseline (0.844 vs 0.904), and tight CV repeat spread made a +0.07 gain look significant when
the bootstrap CI spanned zero.

---

## 6. Coverage of the three judged properties

| Property | Can we predict it? | Can we act on it? |
|---|---|---|
| **Binding affinity** | yes, within ~4 substitutions of measured variants | yes — CDRH3 design, stages 1–2 |
| **Thermostability** | yes, NanoMelt, ±4.1 °C | **no** via CDRH3; framework only, and nothing clears the error bar |
| **Epitope preservation** | **no** | **no** |

Designs are tiered by how many of the three they improve, so this is the pipeline's biggest gap: we can
rank on one property of three.

---

## 6b. Per-position design guidance (from the full 131 legal variants)

From `LEDGER_CROSSCHECK.md`, Mann-Whitney with BH correction over all 10 positions:

- **Spend substitutions on G100 and G101** — the only positions where mutating associates with *tighter*
  binding (q = 0.062, 0.077).
- **Then D102, A106, D108** — inert (q > 0.8), so they are free filler to reach the ≥4 floor.
- **Conserve W99, Y105, F104, M107** — mutating any of them associates with *weaker* binding
  (q < 0.001, < 0.001, 0.028, 0.040).

This is **design guidance, not a better predictor.** The ranker already encodes the same per-position
structure — its implied costs agree with the measured effects at Spearman +0.794 across the 10 positions —
and adding a costly-position counter to it changes nothing out-of-fold (+0.683 → +0.684, CI [−0.005,
+0.011]). Use the table to decide *where to put substitutions when generating* a design; use the ranker to
score the result.

Note that only **1 of 131** legal variants beats trastuzumab's measured KD, and **0 of 131** beat its best
replicate — and all 131 are in the provided data, so none is itself submittable.

## 7. What would actually move the needle

1. **More measured KDs in the ≥5-substitution band.** With 179 labels, differences below ~0.1 Spearman are
   unresolvable, and the ranker is blind exactly where we are allowed to design. This dominates everything
   else; no feature engineering substitutes for it.
2. **Any epitope-preservation signal.** Stage 5 is empty, and it is a third of the score.
3. **A thermostability lever that clears 4.1 °C.** Framework combinations are the only candidate and are
   not submission-ready.

Explicitly *not* worth more effort, on the evidence above: fold-confidence metrics, the 524k FACS table,
and further physicochemical features.

---

## 8. File map

| Path | What it is |
|---|---|
| `eda/AGENT_GUIDE.md` | **how to use the pre-filter** — read before ranking any design |
| `eda/prefilter.py` | the deliverable: `PreFilter().fit()`, `.score()`, `validate()`, `benchmark()` |
| `eda/README.md` | full results and caveats |
| `eda/p1`–`p14_*.py` | the experiments, including `p10_leakage_audit.py` and `p13_bootstrap.py` |
| `eda/seqfeat.py` | T0 feature builder (tested, rejected) |
| `features/README.md` | ddG feature cost tiers, from a different task |
| `thermo/README.md` | thermostability evidence base and the CDRH3 control |
| `LEDGER_CROSSCHECK.md` | cross-check of LEO's ledger against the full SPR-768 table; corrects TH-019 and IDEA-008 |
| `scripts/` | `check_denovo.py`, `check_optimisation.py` (BASE VALIDATION) |

Inputs live in `eda/data/` and are git-ignored; `eda/AGENT_GUIDE.md` §2 has the `scp` commands.
