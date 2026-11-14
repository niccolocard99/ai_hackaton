# Prediction tools, expression/solubility, CDRH3 loop design rules

From the six `cache/biomedcore__*` files on predictors, solubility, ProteinMPNN, cell-free
expression, CDRH3 loops and de novo loop design. 120 records, **none retracted**. **Abstracts only.**

## THE FREE WIN: a published stabilising mutation on trastuzumab scFv

**pmid 37396255 — trastuzumab scFv `VL-R66G`: Tm +5 °C, antigen-binding affinity unchanged.**

Designed by short-MD RMSF rigidification, with the explicit rule "rigidify flexible regions
**outside the CDRs and outside the VH–VL interface**". That makes it **orthogonal to our CDRH3
rewrite** — it cannot cost us binding or epitope preservation. This is the single most defensible
thermostability action available to us.

**Caveats:** single mutation, single paper (14 citations); Tm measured on *E. coli*-expressed scFv,
not a cell-free product; effect measured on **WT CDRs**, so additivity with a heavily rewritten H3
is untested.

## Thermostability predictors — the honest ceiling

**Calibration first (pmid 41724677, Ginkgo AbDev blinded competition, 113 teams, 80 blinded clinical
antibodies):** best **thermostability Spearman ρ = 0.392**; titer ρ = 0.310. Cross-validation scores
"consistently exceeded" held-out test scores — i.e. the field overfits.

**With 4 test slots and a ρ≈0.4 ceiling, a thermostability predictor should drop outliers and break
ties, never fine-rank.** Spend scoring confidence on expression filters instead, where a failed
design yields no measurement at all.

| Tool | Input | Antibody-validated? | Performance | pmid |
|---|---|---|---|---|
| **Antibody-specific pLM + Rosetta-feature CNN** | sequence; seq+Rosetta | **Yes — scFv**, trained on TS50 | OOD Spearman **ρ 0.52** (Ab-specific LM), **0.40** (CNN), vs **0.15** general pLMs. Transferred: recovered **18/20 positions, 5/20 exact substitutions** on an unseen mAb | **36683173** |
| **TEMPRO** | sequence (ESM + NetSurfP3 + AF2 pLDDT) | **Yes — VHH only**, 567 nanobodies | **MAE 4.03 °C, RMSE 5.66 °C** | 39154093 |
| **DeepSP** | sequence only | **Yes — antibody Fv**, 20,530 seqs | SAP/SCM vs MD r **0.76–0.96, mean 0.87** | 38827232 |
| **DeepAb + DMS pipeline** | sequence → Fv structure | Yes — anti-HEL, 200 variants made | **91%** increased thermal+colloidal stability; **10%** had both affinity (5–21×) and Tm (>2.5 °C) gains | 38899735 |
| TemStaPro | sequence | **No** — CRISPR effectors | trained on organism growth temperature | 38507682 |
| ESMStabP | sequence | **No** | R² 0.95 on its own dataset; **bioRxiv preprint** | 40027616 |
| Rosetta Cartesian_ddG | structure | **No** | explicit verdict: **"high false discovery rate"**; Pythia better precision/speed | 40643334 |

**Recommendation (decisive):**
1. **Hard-code `VL-R66G`** (37396255) — the measured, orthogonal +5 °C.
2. **Rank with 36683173** if designs are scFv/Fv (only scFv-trained, scFv-validated ranker with OOD
   numbers); **TEMPRO** only if VHH. **Caveat:** 36683173's abstract does not state availability —
   verify a runnable implementation exists before building a pipeline on it. And note TEMPRO's
   **5.66 °C RMSE is larger than the +5 °C effect** a real stabilising mutation produced, so it
   cannot resolve the differences we care about.
3. **Add DeepSP** as the orthogonal antibody-specific aggregation/viscosity screen (free web app + code).
4. **Do not use** TemStaPro, ESMStabP, or raw Rosetta ΔΔG as filters.

**Any absolute Tm number we quote is indefensible** — there is no validated paired-VH/VL Tm predictor.

## Expression and solubility

| Tool | Antibody evidence | pmid |
|---|---|---|
| **ProteinMPNN, framework-only redesign with CDRs frozen** | **Strongest result in the corpus:** converted **19 of 26** antibody sequences into functional scFv intrabodies — and **18 of those 19 had failed** with the standard approach | **39975170 / 41481731** |
| **Aggrescan3D / A3D 2.0 / A4D** | Yes, experimentally — designed a human VH with significantly reduced aggregation; A4D has an explicit antibody case study (no numbers in abstract) | 25883144, 30825368, 39324697, 30036481 |
| **DeepSP** | Yes — antibody-specific, 20,530 sequences | 38827232 |
| **Protein-Sol** | No antibody validation, but **the only tool trained on E. coli solubility measured in a cell-free system** | 28575391 |
| **SoluProt** | No antibody validation. **Accuracy 58.5%, AUC 0.62** — barely better than a coin flip | 33416864 |
| NetSolP | No antibody validation; no number in abstract | 35088833 |
| **"soluble ProteinMPNN"** | **Appears in this corpus only inside the SPASE abstract. No standalone paper, no benchmark, no antibody validation** | 42560022 |
| SPASE | Benchmarks *predicted* vs *predicted* properties — no experimental validation | 42560022 |

**Note for the team:** the task description's suggested stack (Aggrescan3D, SoluProt, NetSolP,
ProteinMPNN-soluble) is uneven. **A3D is well supported; SoluProt at AUC 0.62 is near-useless as a
filter with only 4 slots** — more likely to discard a good design than catch a bad one. And
"ProteinMPNN with soluble weights" has **no validation record at all** in 120 records.

**Suggested stack:** AF3/ImmuneBuilder model → **A3D 2.0 or A4D** (structure-based aggregation
surface) → **DeepSP** (antibody-specific SAP/SCM) → **Protein-Sol** (cell-free-trained), then freeze
the rewritten H3 and let **ProteinMPNN re-pick framework/non-paratope positions** (per 39975170).

## Cell-free expression — the determinants are real and mostly not sequence

**Base rate to plan around:** de novo CDR-designed scFvs — of 50 synthesised, **"roughly half
resulted in detectable scFv expression"**, and **4 bound** the target (pmid 28860479).
**With 4 slots, expect 1–2 silent failures unless filtered. This, not Tm ranking, is where filtering
buys the most information.**

1. **Disulfide formation is not automatic.** In dialysis-mode CFPS, "**the concentration of
   supplemented disulfide isomerase, DsbC, was critical** to obtaining soluble antibody fragments";
   with it, anti-EGFR Fv and Fab reached >1 mg/mL (pmid **29462206**). *So tighten our rule: keep
   only the two canonical intradomain disulfides and introduce no additional Cys anywhere, including
   in the rewritten H3 — a free thiol in a redox-tuned, isomerase-limited lysate is a mis-pairing and
   aggregation hazard.*
2. **"Expressed" ≠ "soluble".** A 2026 E. coli-CFPS workflow needed a two-stage refolding step to
   recover "functional scFvs **from insoluble fractions**" (pmid 42697344). *The dominant failure
   mode for our format is aggregation, not transcription — so the aggregation-surface filters matter
   more than a generic solubility classifier.*
3. **Proteolysis by the lysate** — serine proteases degrade fragments; ΔdegP ΔompT extract
   significantly increased scFv/Fab yield and activity (pmid 16233776).
4. **Codon usage and 5′ context** — codon optimisation + rare-tRNA extract gave **15.2×** higher
   soluble yield, up to **23×** with full optimisation (pmid 33195157).
5. **Format and chain stoichiometry** — HC:LC plasmid ratio is "a critical parameter" for IgG
   assembly (pmid 36978695). *A single-chain format (scFv or VHH) eliminates this failure mode.*
6. **Eukaryotic lysates need a signal peptide** — in CHO-lysate CFPS, signal-peptide-induced ER
   translocation "**was found to be the prerequisite**" for antibody chain assembly (pmid 28931913).
   **→ Confirm with Adaptyv whether their lysate is bacterial or eukaryotic: the signal sequence is
   mandatory in one case and harmful in the other.**

## CDRH3 loop design rules — the hardest numbers we have

**Length. Cap the rewritten H3 at ~11 residues.** RosettaAntibody median CDR-H3 backbone RMSD by
length (pmid **19062174**): **1.6 Å (4–6 aa), 1.9 Å (7–9), 2.4 Å (10–11), 3.1 Å (12–14), 6.0 Å
(17–22)**. *Beyond ~14 residues our own structure model is wrong by >3 Å, and every score computed on
it — epitope preservation, A3D surface, MD, Tm proxy — degrades with it.*

**Grafting onto a non-native scaffold is the harder regime:** ab initio H3 prediction was 91% <2.0 Å
in the native crystal environment but only **70% <2.0 Å when H3 was swapped onto a homologous
scaffold** (pmid 23255066). Rewriting H3 on trastuzumab's framework is the 70% case.

**Conformation — what to preserve:**
- **Preserve the C-terminal kink.** It "disrupts the β-strand pairing at the base of the loop" and is
  "a critical driver of the observed structural diversity in CDR H3" (pmid 25579815). Constraining H3
  to a kinked conformation "allows near-native conformations to be sampled more frequently" (24519881).
- **Keep the conserved tryptophan at the H3 base.** It "extend[s] the kinked base of these β-bulges
  through a stabilizing intramolecular hydrogen bond and CH/π interaction"; the β-bulge is a
  **molecular hinge** between the rigid stem and the flexible apex (pmid 38916527).
  **→ Concentrate the <70%-identity rewriting in the apex; leave the stem, bulge and Trp intact.**
- **Do not rely on the Arg/Lys94 + Asp101 anchor rule** — "all four combinations lead to a majority
  of conformations that are bulged" (pmid 21035459).
- **Hold L1/L2/L3/H1/H2 at parental lengths and canonical clusters** so H3 is our only source of
  structural uncertainty (~85% of non-H3 sequences are assignable, 21035459).
- **Trastuzumab's VL is an asset** — VL subtype and CDR canonical forms "influence the structural
  diversity of the CDR-H3 loop"; the rigid minor loops "fine-tune the structure of CDR-H3"
  (pmid 37564639). Holding VL fixed sterically constrains our rewritten H3 and makes it predictable.
- **Score over an ensemble, not one pose** — CDR loops interconvert between canonical clusters on the
  **µs–ms timescale, "independent of length and sequence composition"** (pmid 32264741).

**What ML loop design actually achieves:**
- **RFantibody** (RFdiffusion for antibodies, pmid **41193805**): de novo VHHs/scFvs to specified
  epitopes, cryo-EM-confirmed, atomically accurate CDRs — but **"initial computational designs
  exhibit modest affinity (tens to hundreds of nanomolar Kd)"**, needing display screening +
  OrthoRep maturation to reach single-digit nM. *Sobering for a one-shot, 4-design submission.*
- **RFantibody's clash problem, quantified by a competitor** (pmid 42635242): it "frequently
  produce[s] severe steric clashes requiring extensive post-hoc filtering"; TiDE-Ab reached epitope
  recall 0.935 vs RFantibody 0.878 **with over 95% fewer steric clashes**. *Budget for aggressive
  clash filtering.*
- **Use an antibody-specific inverse-folding model, not vanilla ProteinMPNN.** Benchmark (pmid
  40465588): **AntiFold best for Fab, LM-Design most adaptable incl. VHH**, while "models trained on
  general protein datasets (ProteinMPNN, ESM-IF) **struggle with antibody-specific nuances**".
  **Hard warning for our exact task** (pmid 40576473): on CDR3 of nanobody huJ3, AntiBMPNN reached
  **1.7 nM**, AbMPNN 51.2 nM, and **AntiFold and ProteinMPNN gave no detectable activity**.
- **Diversify the input backbone:** MD-cluster-seeded RFdiffusion gave stable designs at **60%** vs
  **30%** for static single-structure design (pmid 41473986).
- **Fragment-based combinatorial loop design** (pmid 36367941): 6 single-domain antibodies, **all
  stable and binding at nanomolar affinity without any in vitro maturation**; runs on a laptop.
  *With a 4-design budget, this precision-first posture fits better than generate-and-screen.*
- **Humanise/fix the framework before designing loops** (pmid 41096532).
- **Non-paratope CDR positions are free design space** — redesigning the loops' "backside" proximate
  to the framework "yield[s] an antibody with binding better than the CDR graft" (pmid 23816785).
  *A legitimate route to <70% H3 identity without touching paratope contacts.*
- **HER2 domain IV de novo precedent exists** (pmid **41301505**): RFdiffusion + ProteinMPNN + AF2
  minibinders against **HER2 domain IV** — the trastuzumab epitope region — ITC-validated, minibinder
  `0_703_6` at nanomolar affinity, ~3× higher affinity and ~3× smaller than existing drugs.
- **Aggregation resistance and structural stability anti-correlate** (r = −0.82, p < 0.01) in a loop
  redesign study that reached 2.3× better aggregation resistance and >98% monomer vs 68% WT
  (pmid 41563298).

## Gaps

1. **No validated Fv/scFv Tm predictor.** TEMPRO is VHH-only; 36683173 predicts TS50 classes, not Tm;
   blinded ceiling on real antibodies is ρ ≈ 0.39.
2. **No predictor of cell-free expression failure.** Protein-Sol is the only cell-free-trained tool
   and its abstract gives no accuracy figure.
3. **The determinants of cell-free failure are reaction/extract/DNA-level** (DsbC, proteases, codon
   usage, RBS, HC:LC ratio, signal peptide), **not properties of our amino-acid sequence.**
4. **Unpaired cysteine is unquantified** — no record measures the effect of a free Cys on antibody
   expression, solubility or Tm. Support for the rule is indirect.
5. **No per-design hit rate for RFantibody.** The only figure ("10–20%") is secondhand in a no-pmid
   conference abstract.
6. **Epitope-preservation scoring is unvalidated for a rewritten H3.** AF3 is the best antibody-complex
   predictor tested yet produced medium/high-quality interfaces for only **46%** of 50 complexes
   (pmid 42719814).
7. **No trastuzumab-specific CDRH3 length or identity-tolerance data.** The only quantitative
   length-tolerance number is 2F5's ±4 residues (22829767) — different antibody, n = 1.
