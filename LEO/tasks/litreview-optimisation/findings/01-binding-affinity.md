# Binding affinity + CDRH3 affinity maturation

From `cache/biomedcore__trastuzumab-affinity-maturation-cdr-mutations.json`,
`...anti-her2-antibody-affinity-improvement-spr-kd-variant.json`,
`...trastuzumab-cdrh3-engineering-her2-binding.json`,
`...antibody-cdr3-randomization-library-affinity-her2.json`.
**Abstracts only — not yet verified against full text.**

Scope arithmetic *(inference, not abstract-stated)*: trastuzumab CDRH3 is ~11–13 residues
(Kabat H95–H102), so **<70% identity means ≥4 simultaneous substitutions**. Any result built on
1–3 point mutations is mechanism-only, never copyable.

## Headline: the constraint matches what actually worked

The literature's best affinity gains came from **many simultaneous CDR substitutions**, not from
single points. Our <70% identity rule is **aligned with** published practice, not opposed to it.

| Evidence | Result | pmid |
|---|---|---|
| **anti-HER2 scFv C6.5**: 20-aa VH CDR3 split into **4 sub-libraries**, randomised separately, winners recombined | VL-CDR3 alone 16-fold (KD 1.0 nM); best VH-CDR3 library 9-fold (KD 0.11 nM); recombining across sub-libraries → **1230-fold total, KD 13 pM** | **8918938** |
| Permissive-position libraries (Ala-scan defined) in a VHH | Greatest gains (**>5-fold**) came from variants carrying **4–6 CDR mutations**. Gains were "generally **not** due to direct interactions involving the mutated residues but indirect effects" enhancing contacts of unmutated residues | **28928732** |
| Multi-CDR perturbation + recombination, anti-ErbB2 scFv HuA21 | **158-fold, Kd 25.5 pM** | 26046845 |
| anti-HER2 mAb 1E11, CDR-L3/H3 randomised | Winner carried **4 aa substitutions** → 10-fold gain (but in **L3**; selection converged away from H3) | 26225765 |
| Unbiased all-six-CDR recombination, anti-ARG2 | Crystal structures show **large reorientation of the paratope**, increased contact surface and shape complementarity | 32616569 |

**Trastuzumab's paratope is demonstrably rewritable.** bH1/bH1-44 — radically re-engineered Herceptin
binding sites — **kept Herceptin-level HER2 affinity and kinetics** despite much greater conformational
plasticity, switching from **enthalpy-driven to entropy-driven** binding without the expected entropic
penalty (**21526167**). Combined with the VHH-at-the-same-epitope result in the epitope findings, this
is the core scientific license for our approach.

## Trastuzumab's own numbers (the baseline to beat)

- **Parent KD 0.058 nM**; computationally designed **triple mutant 32-fold** better, by SPR.
  >90% of intermediate single/double mutants beat parent; only 20–30 mutants made (**28750054**, ADAPT).
  **The winning positions are not in the abstract — this is the single most important fact to pull
  from full text**, both as a position prior and to check the ≥3-mutations-from-published rule.
- **193 single-mutant trastuzumab affinity values** exist, used to train an epistasis-aware model;
  test R²pred = 0.71 on combination mutants. Key mechanism: affinity variation traced to
  **mutation-induced secondary-structure change in CDR-H3**, not direct contact loss (**34392800**).
- A structure-guided **SPR-ranked Herceptin mutant panel** (2 strong / 2 moderate / 3 weak) exists —
  a source of permissive positions *and* of published variants to stay ≥3 mutations from (**31891584**).

## Design rules

**Where to spend substitutions**
- **Split CDRH3 into 3–4 contiguous sub-regions, diversify each, recombine winners.** Gains were
  multiplicative across sub-regions (8918938: 9-fold × 9-fold).
- **Concentrate change at the apex and C-terminal half.** C-terminal CDR3 residues were
  contact-enriched (11335716); the apex controlled specificity in a long H3 (25484051).
- **Treat the H3 base as structural, not tunable.** Residue **H95 sets a ~100-fold "ceiling
  affinity"** and acts as a class switch: mutations elsewhere had **reciprocal** effects between
  Tyr95-type and Gly95-type antibodies (**20961619**).
- **Pick positions by contact/permissiveness, never by SHM hotspot propensity** — randomising 9
  AGY/RGYW hotspot codons gave only 30-fold vs 450–6000-fold from non-targeted approaches (15777944).
- **Changes outside H3 are free against the identity budget** — L3, H2, VH/VL family, framework are
  independent levers (26225765, 26406987, 29593727, 40653772).

**Length**
- Functional CDR3 length spans **7–24 aa** in a VH scaffold (18687338); 14 aa was chosen deliberately
  to give a convex topology (35351202); indel-inclusive recombination beat point mutagenesis
  (20645027: 1000-fold in a single round).
- **But length change rewires epistasis.** The H95 ceiling effect held for a 9-residue H3 and
  **vanished at 11 residues** (20961619). Change length and any inherited "good mutation" list may
  invert sign.

**Residue chemistry**
- **Arginine: place in hydrophobic CDR microenvironments, and don't let affinity depend on it.**
  High-specificity variants put Arg in the most hydrophobic CDR regions; low-specificity ones in
  hydrophilic regions and lost affinity on Arg→Ala (**28778924**).
- **Pair an introduced Arg with an Asp** — Arg–Asp salt bridges gave highly localised rigidity
  increases during real maturation (26132144).
- **Conservative charge geometry beats bulk aromatics.** Asp→Glu at a charged contact gave
  1.5 nM → 0.2 nM (40267400). Adding two interface tryptophans made affinity **6-fold worse**
  (326 nM vs 53 nM parent, **42909185**).

**Hard pre-synthesis filters**
- **No new N-glycosylation sequons (N-X-S/T, X≠P).** Three otherwise-good variants were lost to a
  CDR-H2 glycan that phage/bacterial screening could not see (**24492299**).
- No free cysteines (also a track rule).
- **Avoid engineered H3 disulfide staples** — they form quantitatively but over-constrain the loop
  (11335716).

**Mechanism to aim at**
- Target **conformational pre-organisation of H3 in the bound state**, not raw contact count
  (25524709: as few as two mutations stabilise the loop, but **none is detectable alone**).
- Don't over-weight rigidity: repertoire-scale analysis found only a **slight** H3 rigidification on
  maturation (29545810), and VL/CDR-L2 become *more* flexible (26132144).
- **Design and score combinations, not singles.** Single-site scans missed the key mutations entirely
  (25524709); trastuzumab combination affinity needs an epistasis-aware model (34392800). Note the
  tension: 28750054 found contributions "roughly additive" — don't assume additivity past 3 mutations.

## Validation warnings (read before trusting any ranking)

- **ELISA signal is not affinity.** Both failed designs in **42909185** gave *higher* indirect-ELISA
  signal than parent while being 2–6-fold **worse** by Biacore.
- **Fixed-backbone ΔΔG ranking (FoldX-style) on a modelled complex failed** in that same benchmark —
  it cannot see the repacking cost.
- **Antibody language models lacked framework-mutagenesis insight** vs structure-based methods (40653772).
- **SPR format matters:** capture the antibody via the light chain / Fab region and flow soluble HER2
  over it, to suppress avidity artefacts. Critical here because trastuzumab's HER2 complex is
  sterically crowded and effectively **monovalent**, unlike pertuzumab's (37569408). *(Format advice
  is from a conference abstract with no pmid — unverified.)*

## Read in full, in order

| pmid | Why |
|---|---|
| **28750054** | Trastuzumab's own parent KD and 32-fold designed variants — **need the mutated positions** |
| **8918938** | The 4-sub-library CDRH3 template that gave 1230-fold on HER2 |
| **28928732** | 4–6 mutation variants; indirect-effect mechanism — the license for heavy rewriting |
| **34392800** | Epistasis-aware scoring keyed on H3 secondary structure; 193-point training set |
| **20961619** | H95 ceiling classes and length-dependent epistasis |
| **21526167** | bH1/bH1-44 — how much Herceptin paratope rewriting is tolerated, and the thermodynamic switch |
| **32616569** | Paratope reorientation under extensive sequence change |
| **42909185** | The negative result: trust neither FoldX ranking nor ELISA pre-screen |

Runners-up: 28778924 (Arg placement), 35351202 (per-position codon/length recipe),
24492299 (glycosylation trap), 31891584 (SPR-ranked Herceptin panel).

## Gaps

- **No trastuzumab CDRH3 sequence or per-position contact map in any abstract.**
- **The winning triple-mutant positions (28750054) are not in the abstract** — biggest missing fact.
- **No abstract reports a trastuzumab variant with ≥4 CDRH3 substitutions and a measured KD.** Every
  trastuzumab-specific gain here is 1–3 mutations; all 4+ evidence comes from *other* antibodies.
  **The <70%-identity regime is literally unmapped for trastuzumab in this corpus.**
- No quantitative CDRH3 length-vs-affinity tolerance curve, and none HER2-specific.
- No published hit rate for 4+ simultaneous H3 substitutions.
- No k_on/k_off decomposition for anti-HER2 variants — so no guidance on engineering electrostatic
  steering vs interface burial against domain IV.
- **Epitope mismatch risk:** several high-affinity anti-HER2 results are at *other* epitopes
  (domain II 40267400, domain II/III 28531814, p95HER2 42909185). Transferability to the domain-IV
  paratope is not addressed by any abstract.

No retracted records. Three load-bearing items are conference abstracts with no pmid (the SPR
format advice, a JCO 2006 ADCC abstract, SABCS P2-19-01).
