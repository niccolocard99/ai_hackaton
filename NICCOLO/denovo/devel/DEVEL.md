# Developability + thermostability scoring for the de novo HER2 VHH campaign

P0-worker (g), verified on this VM 14 Nov 2026. Everything below was run, not assumed.

**Output rule for this whole campaign: predicted and measured values never share a
column.** Every column `devel.py` emits carries a `pred_` prefix, and the CSV header
comment names what produced it. A lab number must arrive under a different prefix
(`meas_`).

---

## 1. The four Adaptyv-recommended tools: what actually happened

| Tool | Installed? | Exactly what was run | Outcome / failure mode |
|---|---|---|---|
| **ProteinMPNN (soluble weights)** | **YES — working** | `git clone --depth 1 https://github.com/dauparas/ProteinMPNN.git` into `…/scratchpad/pmpnn_dl/`; ran with the existing torch 2.14.1+cu126 in `/home/hackathon/.venvs/abmodels` (CUDA available) via `PYTHONPATH=<repo> python protein_mpnn_run.py … --use_soluble_model` | Works. All 4 soluble checkpoints (`soluble_model_weights/v_48_002/010/020/030.pt`, 6.7 MB each) ship in the repo — no separate weight download. Verified end to end on trastuzumab VH. No new venv needed. |
| **Aggrescan3D** | **NO** | `uv venv --python 3.12 /home/hackathon/.venvs/devel_a3d`; `uv pip install aggrescan3d` | Hard fail. The PyPI package is **Python 2 only** — its `setup.py` contains a bare `print "data folder not available…"`, so the build dies with `SyntaxError: Missing parentheses in call to 'print'`. No Python 2 on this VM. A3D also needs FoldX (licensed) or CABS-flex for its dynamic mode, neither present. |
| **SoluProt** | **NO** | `curl https://pypi.org/pypi/soluprot/json` → **404** | Not on PyPI. Distributed only as a Loschmidt-Laboratories web server / standalone tarball; not installed, not reachable as a package. |
| **NetSolP** | **NO** | `curl https://pypi.org/pypi/netsolp/json` → **404** | Not on PyPI. DTU web server + a standalone bundle whose weights are an ESM-1b/ProtBert fine-tune (multi-GB); not installed. |

### Premises in the brief that were WRONG

1. **The hosted Nebius `proteinmpnn` is NOT a route to soluble weights.** Its schema
   (`forge.py schema proteinmpnn`) describes itself as the *"ProteinMPNN **v_48_020**
   portable adapter"* and accepts only `input_pdb`, `input_pdb_chains`,
   `num_seq_per_target`, `sampling_temp`, `random_seed`, `omit_AAs`. **There is no
   weights / `use_soluble_model` field.** Hosted = vanilla weights only. The
   *local* install is the only way to get the soluble weights, and it is also the
   easier route here because torch is already present. The `proteinmpnn` skill
   pointed at the hosted model, which would have silently given us vanilla scores.
2. **The caplacizumab VHH cannot be extracted from
   `apheris/rfd4/Data/03-vhh/pmvk_caplacizumab.cif.zst`.** That file's VHH entity is
   **CDR-stripped**: `_entity_poly.pdbx_seq_one_letter_code` for entity 2 reads
   `EVQLVESGGGLVQPGGSLRLSC` **`X`** `WFRQAPGKGRELVA` **`X`** `YPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYC` **`X`** `WGQGTQVTVS`
   — three `X` placeholders where CDR1/2/3 belong, 83 observed residues instead of ~128.
   It is an RFD4 *design scaffold* against PMVK (chain A = phosphomevalonate kinase),
   not a complete caplacizumab structure. The sanity check therefore uses the
   published caplacizumab VHH sequence (in `reference_knowngood.csv`); its
   framework-2 window `WFRQAPGKEREFVA` matches the CIF's framework exactly.
3. **`pmvk_caplacizumab.cif` also decompresses to a different target than the name
   suggests** — chain A is PMVK, so it is a PMVK binder scaffold, not a HER2 or vWF
   complex. Also note chain B of `her2_fv_1n8z.cif` really is VH 1–120 and chain A
   really is VL 1–107 — that premise was correct.
4. Minor: the brief's `C22–C96` pair is **not** IMGT numbering, it is Kabat. IMGT
   calls the same pair C23–C104. `devel.py` finds it at sequence positions 22 and 96
   in both references and the structure confirms it closed at **2.04 Å** — so the
   numbers in the brief are right, the label is not.

### TEMPRO: found, not usable. We have NO real Tm predictor.

Repo located at `github.com/Jerome-Alvarez/TEMPRO` (Alvarez & Dean, *Sci Rep* 14:19074,
2024; reported MAE 4.03 °C — the only published nanobody Tm model). Cloned; it ships
only archives. `user.zip` unpacks to six pre-trained Keras models
(`ESM_650M.keras`, `ESM_3B.keras`, `ESM_15B.keras`, `af2.keras`, `nsp3.keras`,
`pchars_maestro.keras`, 2–10 MB each, dated Apr 2024).

Two independent blockers, both verified:

1. **The models will not load.** `uv venv --python 3.12 /home/hackathon/.venvs/devel_tempro`
   + `uv pip install tensorflow-cpu` gives TF 2.21.0 / **Keras 3.15.1**;
   `keras.saving.load_model("ESM_650M.keras", compile=False)` fails with
   `ValueError: File not found … Please ensure the file is an accessible .keras zip file`
   **even though `os.path.exists` is True** — i.e. Keras 3 does not recognise these
   Keras-2.10 legacy archives. The repo README says so itself: *"it is important to
   install the correct Python version (version 3.9) to properly load the Keras
   models."* `requirements.txt` pins `python==3.9.15, tensorflow==2.10.1,
   keras==2.10.0`. TF 2.10 has no Python 3.12 wheels, and this VM has only Python
   3.12 and no conda/mamba on PATH (`/home/hackathon/.conda` and `.mamba` are empty
   config dirs; the only binary is `/home/hackathon/nick-work/bin/micromamba`, which
   was not exercised inside the time box).
2. **Even if they loaded, the features are missing.** The model's input is a
   mean-pooled **ESM-2 embedding** (650M = 1280-d, 3B = 2560-d, 15B = 5120-d). We
   would need `fair-esm` plus the ESM-2 checkpoint (2.5 GB for 650M, 11 GB for 3B)
   and the embedding pipeline from `embedding_generator.zip`. The author also states
   the 15B embeddings were never uploaded.

**Conclusion: there is no working thermostability predictor on this VM.** Confirming
the brief: Nebius `stab-ddg` predicts *binding* ddG for point mutants (an affinity
proxy, not folding stability) and JapanFold `opendde`/`opendde-abag` are docking
models. Neither is a Tm predictor. `pred_tm_c` is therefore **empty** and
`pred_tm_source` reads `no_tm_predictor_available`.

**Recovering a real Tm number is the single highest-value follow-up.** Cheapest
credible path, in order: (a) `micromamba create -n tempro python=3.9 tensorflow==2.10.1
keras==2.10.0` using `/home/hackathon/nick-work/bin/micromamba`, then re-try the
`.keras` load; (b) generate the 650M embeddings either with `fair-esm` locally or via
the hosted Nebius `esm2-650m` / `esmc-600m` models (both advertise embeddings from a
FASTA in a tar); (c) only then is `pred_tm_c` populated via `--tempro-csv`.

---

## 2. What replaced each tool

| Missing tool | Proxy now implemented in `devel.py` | Honest limitation |
|---|---|---|
| Aggrescan3D | **Structural exposed-apolar-patch metric**: Shrake-Rupley SASA (biotite) → apolar residues (AVILMFWYC) with relative SASA > 0.25 (Tien 2013 max-SASA reference) → single-linkage clustering on CB–CB ≤ 8 Å → largest cluster's residue count and total SASA (Å²), plus total exposed apolar SASA. | Measures apolar **surface contiguity**, not an aggregation free energy. Direction of effect agrees with A3D (bigger patch = worse); the units and the per-residue structure-corrected scale do not. Reported as `pred_hydrophobic_patch_*`, never as an A3D score. |
| Aggrescan3D, no structure | `pred_hydrophobic_patch_seq_proxy` = max mean Kyte-Doolittle over a 5-residue sliding window, with its start position. | Explicitly labelled a proxy in its own `_note` column. Cannot know whether the window is buried — strictly a screen. Use the structural metric whenever a backbone exists. |
| SoluProt / NetSolP | **ProteinMPNN soluble-weight score** (the real tool, run locally) as the expression/solubility prior, plus the charge/GRAVY/patch block: `pred_net_charge_ph74`, `pred_pI`, `pred_gravy`, patch metrics. | ProteinMPNN's soluble score is a *designability-under-a-soluble-prior* log-likelihood, not a calibrated solubility probability. It ranks; it does not predict percent-soluble. |
| TEMPRO | `pred_tm_surrogate_rank` — an ordinal 0–100 rank from the disulfide, the FR2 hallmark count, CDR-H3 length, patch metric and GRAVY. | **Never fitted to any Tm dataset.** It is a judgement-weighted rank, not a temperature; its `_note` column says so. Only valid for ranking designs against each other inside this campaign. `pred_tm_c` stays empty. |

### Everything `devel.py` computes

`devel.py --csv designs.csv --out scores.csv [--structure-dir DIR | --structure P] [--tempro-csv P]`
on a CSV of `sequence_id, kind, heavy_sequence, light_sequence` → 59 columns:

- **Charge/hydropathy**: `pred_net_charge_ph74` (Henderson-Hasselbalch, EMBOSS pKa set,
  termini included), `pred_pI` (bisection), `pred_gravy` / `pred_mean_hydrophobicity_kd`.
- **Patch**: `pred_hydrophobic_patch_n_res`, `_sasa_A2`, `_residues`,
  `pred_exposed_apolar_sasa_A2` (structure) and `pred_hydrophobic_patch_seq_proxy`
  (+`_start`, +`_note`) always.
- **CDR composition**: `pred_cdr_frac_arg`, `pred_cdr_frac_lys`,
  `pred_cdr_frac_hydrophobic`, `pred_cdr_net_charge_ph74`, `pred_cdrH1/H2/H3`,
  `pred_cdrL1/L2/L3`, `pred_cdrH3_len`, `pred_cdr_total_len`, `pred_cdr_parse_ok`.
  CDRs are located by **regex anchors on the conserved framework motifs, not ANARCI**;
  the CDR3 definition (conserved Cys+1 … W of WGxG−1) is IMGT-like, so trastuzumab
  CDR-H3 comes out as the 13-mer `SRWGGDGFYAMDY` rather than the Kabat 11-mer.
  **`pred_cdr_parse_ok == 0` means the anchors missed — do not trust the CDR columns
  on that row.**
- **Liabilities with 1-based positions** (`pred_liab_detail`, chain-tagged, `(CDR)`
  marked): N-glycosylation `N[^P][ST]`, deamidation `NG` and `NS`, isomerisation `DG`,
  acid-labile `DP`, total/odd/extra cysteine, Trp (total and in-CDR, the oxidation
  risk), Met (total and in-CDR), and `pred_liab_total`.
- **VHH-specific**: `pred_vhh_hallmark_n` (0–4) and `pred_vhh_hallmark_detail` reading
  Kabat 37/44/45/47 positionally from the conserved W (indexed as Kabat 36), with
  `pred_fr2_window`; `pred_vhh_like`; `pred_cys_fr1_pos` / `pred_cys_fr3_pos` /
  `pred_canonical_ss_seq_ok` for the C22–C96 (Kabat; IMGT C23–C104) pair, and with a
  structure `pred_canonical_ss_struct_closed`, `pred_struct_n_disulfide`,
  `pred_struct_free_cys` and the measured SG–SG distances in `pred_struct_note`.
- **Tm**: `pred_tm_c` (+`pred_tm_source`), `pred_tm_surrogate_rank` (+`_note`).

### ProteinMPNN soluble score — the exact verified command

```bash
PYTHONPATH=<SCRATCH>/pmpnn_dl/ProteinMPNN \
/home/hackathon/.venvs/abmodels/bin/python \
  <SCRATCH>/pmpnn_dl/ProteinMPNN/protein_mpnn_run.py \
  --pdb_path design.pdb --pdb_path_chains "B" \
  --use_soluble_model --score_only 1 \
  --out_folder out/ --seed 1 --batch_size 1
```
Writes `out/score_only/<name>_pdb.npz`. Drop `--score_only 1` and add
`--num_seq_per_target N --sampling_temp 0.1` to redesign. Run the same backbone with
and without `--use_soluble_model`: the **gap** between vanilla and soluble score is
the informative quantity (a design the soluble prior dislikes more than the vanilla
prior is the expression risk). Report these as `pred_mpnn_soluble_score`,
`pred_mpnn_vanilla_score`, `pred_mpnn_soluble_gap` — they are not yet columns of
`devel.py` because the scorer needs a per-design backbone; wire them in when the
design PDBs exist.

---

## 3. Reference numbers — the known-good anchors

Both pass, as they must: trastuzumab is an approved IgG1 and caplacizumab an approved
VHH. **These are the calibration points; a threshold that fails either is wrong.**

| Metric | trastuzumab Fv (VH 120 + VL 107, from `her2_fv_1n8z.cif` B+A) | caplacizumab VHH (128 aa, published seq) |
|---|---|---|
| `pred_net_charge_ph74` | **+2.82** | **+1.78** |
| `pred_pI` | **8.48** | **8.64** |
| `pred_gravy` | **−0.334** | **−0.381** |
| `pred_cdrH3` / len | `SRWGGDGFYAMDY` / **13** | `AAAGVRAEDGRVRTLPSEYTF` / **21** |
| `pred_cdr_frac_arg` | **0.060** | **0.104** |
| `pred_cdr_frac_hydrophobic` | **0.388** | **0.354** |
| `pred_cdr_net_charge_ph74` | +1.2 | 0.0 |
| `pred_liab_total` | **5** (NG@H55 CDR, NS@H84, DG@H102 CDR, Trp@H99 CDR, Met@H83/H107) | **3** (NS@H84, DG@H105 CDR, Met@H34 CDR) |
| `pred_liab_nglyc` / `cys_odd` / `cys_extra` | 0 / 0 / 0 | 0 / 0 / 0 |
| `pred_vhh_hallmark_n` (FR2 window) | **1** (`WVRQAPGKGLEWVA`; K37=V, K44=G, K45=L, K47=W — all VH) | **4** (`WFRQAPGKEREFVA`; K37=F, K44=E, K45=R, K47=F) |
| `pred_cys_fr1_pos` / `pred_cys_fr3_pos` | **22 / 96** | **22 / 96** |
| `pred_hydrophobic_patch_sasa_A2` (struct) | **449.8 Å² over 5 residues** (`L11;V15;F83;Y105;I106`) | n/a (CIF scaffold is CDR-stripped) |
| `pred_struct_n_disulfide` / free Cys | **2** / 0 — `A23–A88 @ 2.09 Å` (VL), **`B22–B96 @ 2.04 Å`** (VH) | 1 / 0 in the stripped scaffold (`H22–H75 @ 2.05 Å`) |
| `pred_tm_c` | *(empty)* `no_tm_predictor_available` | *(empty)* `no_tm_predictor_available` |
| `pred_tm_surrogate_rank` | 60.7 | 56.5 |
| ProteinMPNN **soluble** score (VH chain B) | **mean 1.2610, global 1.2768** | not run (no usable backbone) |

Neither is flagged as a disaster, and the hallmark count cleanly separates VH (1) from
VHH (4). The `pred_tm_surrogate_rank` values (60.7 vs 56.5) illustrate why that column
is only a rank: caplacizumab is in reality a *very* stable VHH, and the surrogate
scores it slightly below an Fv purely because of its 21-residue CDR-H3. **Do not read
it as "caplacizumab is less stable than trastuzumab."**

---

## 4. Recommended funnel thresholds

Provisional — each is anchored on the two references above, which means most are
calibrated on **one** data point per format. Use them to rank and to catch outliers,
not as physics.

| Gate | Pass | Watch | Fail | Why this number |
|---|---|---|---|---|
| `pred_cdr_parse_ok` | 1 | — | 0 | A row whose CDRs could not be located has meaningless CDR columns. Quarantine, re-run with ANARCI if it matters. |
| `pred_nonstandard_residues` | empty | — | any | `X`/gap means the design or its extraction is incomplete — exactly the trap the caplacizumab CIF sets. |
| `pred_net_charge_ph74` | −2 … +6 | +6 … +9 | <−4 or >+9 | Refs sit at +1.8 and +2.8. Jain 2017's developability survey of clinical mAbs clusters in this band; strongly positive surfaces drive non-specific tissue binding and fast clearance, strongly negative ones hurt expression. |
| `pred_pI` | 6.0 … 9.0 | 5.5–6.0, 9.0–9.5 | <5.5 or >9.5 | Refs 8.48 / 8.64. A pI within ~1 unit of a pH 6–7 formulation is the classic aggregation window; pI > 9.5 is the polyreactivity end. |
| `pred_gravy` | ≤ −0.25 | −0.25 … −0.10 | > −0.10 | Refs −0.334 / −0.381. Both approved molecules are comfortably hydrophilic; a VHH drifting above −0.1 has an unusual surface for the fold. |
| `pred_hydrophobic_patch_sasa_A2` **VHH, with structure** | ≤ 450 | 450 … 700 | > 700 | The two-chain trastuzumab **Fv** — twice the mass and with a VH/VL interface — reaches only 449.8 Å². A single-domain VHH exceeding that is carrying excess exposed apolar surface. Re-calibrate the moment you have 5+ real VHH backbones. |
| same, **Fv/scFv format** | ≤ 600 | 600 … 900 | > 900 | 1.3× / 2× the measured trastuzumab Fv value, to allow for scFv linker exposure. |
| `pred_hydrophobic_patch_seq_proxy` (no structure) | ≤ 2.6 | 2.6 … 3.2 | > 3.2 | Trastuzumab 2.16, caplacizumab 3.02 — note caplacizumab, an approved drug, is at 3.02, which is exactly why this sequence proxy must never be a hard gate. Use it only to triage which designs get folded. |
| `pred_cdr_frac_arg` | ≤ 0.12 | 0.12 … 0.18 | > 0.18 | **Caplacizumab is 0.104**, so any gate at or below 0.10 would reject an approved drug. Above ~0.18 the Arg density is the known non-specific/polyreactive-binder signature — Arg-driven affinity is charge-driven, not shape-driven, and tends not to survive a specificity counter-screen. |
| `pred_cdr_frac_hydrophobic` | ≤ 0.45 | 0.45 … 0.55 | > 0.55 | Refs 0.388 / 0.354. Above ~0.55 the paratope is a hydrophobic sticker — the dominant de novo failure mode for both aggregation and off-target binding. |
| `pred_liab_nglyc` | 0 | — | ≥ 1 | A new NxS/NxT sequon in a CDR changes the molecule the lab makes. Both refs are 0. Free to fix. |
| `pred_liab_cys_odd` + `pred_liab_cys_extra` | 0 | — | ≥ 1 | Both refs have exactly the framework pair per domain and zero free Cys. An odd or extra Cys gives disulfide scrambling and covalent dimers — a direct expression and measurability killer. |
| `pred_canonical_ss_seq_ok` and `pred_canonical_ss_struct_closed` | 1 | — | 0 | The C22–C96 (Kabat; IMGT C23–C104) pair is the single largest contributor to VHH Tm. Confirmed closed at 2.04 Å in trastuzumab VH. Lose it and the thermostability score is gone. **Hard fail.** |
| `pred_vhh_hallmark_n` (VHH designs only) | ≥ 2 | 1 | 0 | Caplacizumab has 4/4, trastuzumab's VH 1/4. The hydrophilic tetrad is what lets a VHH live without a VL; a design at 0 has a VH-like, exposed former-VL interface — the classic VHH aggregation route. Not all VHHs have 4 (many natural ones have 2–3), hence the gate at 2. |
| `pred_cdrH3_len` | 8 … 22 | 22 … 26 | > 26 | Caplacizumab is 21 (IMGT-like), so the pass band must include it. Long CDR-H3s raise entropy cost and lower Tm; past ~26 they also become hard to fold confidently. |
| `pred_liab_trp_cdr` / `pred_liab_met_cdr` | ≤ 1 each | 2 | ≥ 3 | Trastuzumab carries 1 CDR Trp and 1 CDR Met and is a marketed product, so these cannot be zero-tolerance; ≥3 oxidisable CDR residues is a real stability-on-the-shelf liability. |
| `pred_liab_total` | ≤ 5 | 6 … 8 | ≥ 9 | Trastuzumab scores 5, caplacizumab 3. The gate is set *at* the worse approved reference deliberately — a de novo design has no clinical track record to spend. |
| `pred_mpnn_soluble_gap` (once wired) | — | — | — | **No threshold yet — do not invent one.** Rank-order only until we have ≥20 designs to see the distribution. |
| `pred_tm_c` | — | — | — | Empty for every design. **There is no Tm gate.** |
| `pred_tm_surrogate_rank` | rank only | — | — | Never gate on it. Use it to order ties after the hard gates. |

Suggested order: hard structural/chemical fails first (`nonstandard`, `cys_odd/extra`,
`canonical_ss`, `nglyc`, `cdr_parse_ok`), then the VHH-format gates
(`vhh_hallmark_n`, `cdrH3_len`), then the continuous developability band
(charge, pI, GRAVY, patch, CDR composition), then rank the survivors by
ProteinMPNN soluble score and `pred_tm_surrogate_rank`.

---

## 5. What we CANNOT assess

- **Absolute thermostability.** No Tm, no ΔG of unfolding, no T<sub>agg</sub>, no DSF
  prediction. `pred_tm_c` is empty by design. Any Tm number quoted for these designs
  would be fabricated. This is the campaign's biggest open gap, and one of the three
  judged properties.
- **Calibrated aggregation propensity.** No Aggrescan3D score, no CamSol intrinsic
  solubility, no TANGO/WALTZ amyloid call. We have apolar-surface contiguity only.
- **Calibrated solubility / expression probability.** No SoluProt, no NetSolP. We can
  rank with ProteinMPNN's soluble prior; we cannot say "this will express."
- **Viscosity, self-association (AC-SINS/kD), polyreactivity, PK, immunogenicity
  (no T-cell epitope / MHC-II prediction), chemical-stability kinetics** (we flag
  motifs, not rates; a DG in a rigid framework loop is not the same risk as a DG in an
  exposed CDR tip), **conformational heterogeneity, and host-cell-specific expression.**
- **Numbering fidelity.** No ANARCI/IMGT on this VM, so CDR boundaries are regex
  estimates and the hallmark positions are read positionally from a single anchor. A
  heavily non-germline de novo framework can break the anchor; that is what
  `pred_cdr_parse_ok` is for.
- **Thresholds themselves.** Calibrated on two molecules. They will mis-rank.

## 6. Which judged properties this tooling speaks to

| Judged property | Covered? | By what |
|---|---|---|
| **Binding** | **No — nothing here predicts affinity.** | `devel.py` deliberately contains no affinity term. The only partial, *negative* contribution: high `pred_cdr_frac_arg` / `pred_cdr_frac_hydrophobic` flags affinity that is likely charge- or stickiness-driven and so unlikely to survive a specificity screen. Binding must come from the folding/docking/affinity models (Boltz-2, AF2-IG, BoltzGen), not from here. |
| **Epitope preservation** | **No.** | No epitope, interface or contact-map term whatsoever. `devel.py` never sees the target. Epitope checks belong to the complex-structure lane (interface contacts vs. the trastuzumab/1N8Z epitope). |
| **Thermostability** | **Partially, and only as proxies.** | Speaks to it: the C22–C96 disulfide check (hard gate, structure-confirmed), the VHH FR2 hallmark count, CDR-H3 length, GRAVY and the exposed apolar patch — all literature-supported *directions* on VHH Tm. Does **not** speak to it: any actual temperature. No installed predictor produces one. |

**Net: this tooling is a developability and manufacturability filter plus a
thermostability *proxy*. It de-risks "the design cannot be measured at all", which is
its real job — a design that does not express scores zero on all three judged
properties. It does not score binding or epitope preservation, and it does not give a
Tm.**

---

## 7. Files

- `devel.py` — the scorer. Stdlib-only for the sequence block; needs numpy + biotite
  for the structure block (present in `/home/hackathon/.venvs/abmodels`, which also
  carries the torch used for ProteinMPNN). Run it as
  `/home/hackathon/.venvs/abmodels/bin/python -I devel.py …`.
- `reference_knowngood.csv` — trastuzumab Fv + caplacizumab VHH, the calibration set.
  Re-run it after any change to `devel.py`; if either molecule starts failing gates,
  the change is wrong.
- Venvs created (outside the repo): `/home/hackathon/.venvs/devel_tempro`
  (TF 2.21 / Keras 3.15, kept as evidence of the `.keras` load failure) and
  `/home/hackathon/.venvs/devel_a3d` (empty — the aggrescan3d build failed).
  `uv` at `/home/hackathon/.local/bin/uv` is the working package manager on this VM:
  **there is no system `pip` and `python3 -m venv` fails** (`ensurepip` missing), so
  use `uv venv` / `uv pip install --python <venv>/bin/python`.
