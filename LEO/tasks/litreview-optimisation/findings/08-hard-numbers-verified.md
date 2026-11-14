# Hard numbers from full text — and a correction to my own advice

Mined from the 195-record full-text corpus. **This supersedes parts of
[findings 06](06-structural-paratope-measured.md).** Read §2 before designing anything.

## 1. CDRH3 sequence independently confirmed — twice

My structural derivation from 1N8Z was `WGGDGFYAMDY` (11 aa). Two independent sources agree:

- **Patent US-8840896-B2** (Genentech, "Anti-HER2 antibody variants"), **SEQ ID NO: 24**, for
  huMAb4D5-8 = trastuzumab, explicitly "according to both a standard sequence definition (Kabat 1987)
  and a structural definition (Chothia & Lesk 1987)".
- **pmid 38717571** (NMR assignment of the trastuzumab scFab) gives the complete chain sequence.

Full CDR set (Kabat/Chothia definition, per US-8840896-B2):

| CDR | Sequence | Len | SEQ ID |
|---|---|---|---|
| **H3** | **`WGGDGFYAMDY`** | 11 | 24 |
| H2 | `RIYPTNGYTRYADSVKG` | 17 | 23 |
| H1 | `GFNIKDTYIH` | 10 | 22 |
| L1 | `RASQDVNTAVA` | 11 | 19 |
| L2 | `SASFLYS` | 7 | 20 |
| L3 | `QQHYTTPPT` | 9 | 21 |

**Numbering crosswalk (my arithmetic from the patent's named Kabat positions):**

| Kabat | W95 | G96 | G97 | **D98** | **G99** | F100 | Y100a | A100b | M100c | D101 | Y102 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **1N8Z auth** | **W99** | G100 | G101 | **D102** | **G103** | **F104** | **Y105** | A106 | M107 | D108 | Y109 |

Do not use `US-6800738-B1` for CDR text — its parsed sequences are **OCR-corrupted** (`SASFLES`,
`FYAMDV`).

## 2. ⚠ CORRECTION: a measured alanine scan overturns two of my calls

**A real alanine scan exists** (Kelley & O'Connell, *Biochemistry* 32:6828-35 (1993), reported in
US-8840896-B2), with K_D fold-effects — this is the **energetic** budget I said was missing. Combined
with my geometric analysis:

| Kabat | 1N8Z | aa | Ala-scan K_D effect | My geometric call | **Verdict** |
|---|---|---|---|---|---|
| W95 | W99 | W | **50–5000×** | contact, untouchable | **UNTOUCHABLE — confirmed hotspot** |
| G96 | G100 | G | not in scan | free | free |
| G97 | G101 | G | not in scan | free | free |
| **D98** | **D102** | **D** | **only 1.5–2×** | "avoid, contacts at 5.5 Å" | **❌ I WAS WRONG — see below** |
| **G99** | **G103** | **G** | not in scan | untouchable (shortest contact) | **UNTOUCHABLE — and now for a second reason** |
| F100 | F104 | F | **only 1.5–2×** | lowest-risk contact | **confirmed lowest-risk contact** |
| **Y100a** | **Y105** | **Y** | **50–5000×** | "high risk" | **❌ UPGRADE TO UNTOUCHABLE — hotspot** |
| A100b | A106 | A | not in scan | free | free |
| M100c | M107 | M | not in scan | free | free |
| D101 | D108 | D | not in scan | free | free |
| Y102 | Y109 | Y | not in scan | free | free |

**The four critical paratope residues by alanine scan (50–5000× K_D effects): heavy R50, W95, Y100a;
light H91.** Two are in CDRH3 (W95/W99 and Y100a/Y105). Light-chain H91 alone costs **200-fold** when
mutated to Ala, and has "only 22% exposed surface area" — another reason the light chain must stay fixed.

### Correction A — `D102` is not a position to avoid. It is the best position to mutate.

I told you to avoid `D102` because it picks up contacts at a 5.5 Å cutoff. The measured data says the
opposite, and emphatically:

- **Ala-scan effect only 1.5–2×** — it is not a hotspot.
- US-8840896-B2: D98 was **"the most poorly conserved residue"**, mutated in **98% of selected phage
  clones**, most frequent substitution **Trp (23%)**.
- **`D98W` gives a 3-fold improvement in binding affinity over wild-type hu4D5-8 Fab, measured by
  SPR** (BIAcore, 37 °C).
- `D98F` also gains affinity (pmid 31891584), by removing electrostatic repulsion with HER2 D586/E598.
- It sits at **"the tip of variable loop 3"** — the apex, exactly where the loop-design literature says
  to concentrate change.

**So `D102` belongs in the mutable set, and `D102W` or `D102F` is a known affinity-improving move.**
**Caveat:** `D98W` is published, so a design must stay ≥3 mutations from that single mutant — fine if
`D102W` is combined with 3+ other substitutions, but a design whose only change is `D102W` is invalid.

### Correction B — `G103` is untouchable for two independent reasons

My geometric reason: it makes the interface's shortest contact (2.42 Å to HER2 K593) and is a glycine
with no side-chain room. The patent adds the structural reason: **"V_H99 was a highly conserved residue
in variable loop 3… position i+2 of a type II β-turn which is most commonly a Gly… likely to be a
structure-stabilizing residue."** Two independent lines of evidence, same conclusion.

*(Note one dissent: pmid 31891584 proposed `G99D/R` as a deliberate "size increase in the CDR-H3 loop"
to modulate electrostatics with HER2 D570. That is a real published idea, but it contradicts both the
β-turn argument and the 2.42 Å contact. Treat as high-risk/high-reward, not a default.)*

### Revised mutable set

| Tier | Positions (1N8Z numbering) | Rationale |
|---|---|---|
| **Preferred** | **D102** | Apex, 98% mutated in selection, D→W/F gives measured affinity gain |
| **Free** | G100, G101, A106, M107, D108, Y109 | No HER2 contact at any cutoff; absent from the Ala scan |
| **Usable with care** | F104 | Contact, but only 1.5–2× Ala effect |
| **Never** | W99, G103, Y105 | Two Ala-scan hotspots (50–5000×) + one structural Gly |

**7 preferred/free positions; 5 substitutions is comfortably achievable** (≥4 needed for <70%).

## 3. THE Tm BASELINE — ⚠ CORRECTED: it is literature hearsay, not a measurement

An earlier pass of this file stated "Herceptin Fab T_M = 83 °C **by DSC**". **That overstated it.**
A second independent mining pass found that in pmid 21526167 the 83 °C appears only as a comparator —
*"slightly lower than that of Hercepin Fab (83 °C)"* — **cited to an earlier study whose reference is
stripped from the retrieved text. It was not measured in that paper.**

What **was** measured there, by DSC (Microcal; 10 mM NaOAc pH 5 + 150 mM NaCl, 0.5 mg/mL, 1 °C/min):

| Fab | T_m | Relationship to Herceptin |
|---|---|---|
| Herceptin Fab | **83 °C** | **CITED, not measured — treat as hearsay** |
| bH1 | **77 °C** | measured |
| bH1-81 | **76 °C** | measured |
| bH1-44 | **74 °C** | measured |

**So we still do not have a trustworthy trastuzumab fragment Tm.** Other values found:
~73–76 °C for formulated Herceptin **IgG** by far-UV CD at 205 nm (pmid 31071101); WT trastuzumab IgG1
DSF curves exist in pmid 40653772 but **the numbers are figure-only**.

**No trastuzumab Fv or scFv T_m exists anywhere in the corpus.** pmid 32380351 (nanoDSF) and
pmid 37396255 (the `VL-R66G` +5 °C scFv result) are **both ABSENT** from the full-text corpus — a sweep
for "R66G" returned zero hits. The `VL-R66G` claim therefore rests on its **abstract alone** and is
unverified in detail.

### What the paratope rewrite costs in stability — now quantified

| System | Affinity change | Stability cost | pmid |
|---|---|---|---|
| Herceptin → bH1 → bH1-44 (Fab) | 20 nM → **0.2 nM (100×)** | T_m **77 → 76 → 74 °C** (~3 °C across the series; ~6–9 °C below the cited parent) | 21526167 |
| anti-HER2 **DARPin**, directed evolution | **>700×** | **−30 °C** apparent T_m | 30666176 |
| VH single domain, triple mutant | increased | **−18 °C** | 30666176 |

**The directly relevant one is the Herceptin lineage: ~100-fold affinity gain for ~3 °C.** The DARPin
and VH numbers show how bad it gets when affinity is pushed hard without stability counter-selection.

**Crucially, pmid 30666176 also reports the mitigation:** in a 12-mutation VH variant the dominant
affinity mutation **R100d, in CDR3, had *no* stability cost**, while the destabilising one (N72) was
outside CDR3. Its stated conclusion is that **CDR-H3-only libraries carry the lowest
affinity/stability trade-off risk.** That is direct support for confining our changes to CDRH3.

**And every non-trade-off success in the corpus worked *outside* the paratope:** +19.6 °C with unchanged
affinity and 27× yield from framework-only mutations (32326443); +4.6 °C and full survival of 2 h at
70 °C (32336757); DARPin `D17L` **+16.1 °C** with the paratope untouched (34793836).

## 4. Full thermodynamics of Herceptin–HER2 (pmid 21526167, ITC + SPR)

| Quantity | Value | Method |
|---|---|---|
| K_D | **0.5 nM** | SPR, BIAcore 3000, 30 °C, 1:1 Langmuir |
| k_on | **7 × 10⁵ M⁻¹s⁻¹** | SPR, 30 °C |
| k_off | **4 × 10⁻⁴ s⁻¹** | SPR, 30 °C |
| **ΔH** | **−13.5 kcal/mol** | ITC, VP-ITC, PBS pH 7.4, 30 °C |
| **−TΔS** | **+0.9 kcal/mol** | from ΔG(SPR) and ΔH(ITC) |
| ΔH corrected for protonation (~0.5 H⁺/event) | −9.2 kcal/mol | ITC buffer dependence |
| −TΔS corrected | −3.4 kcal/mol | ″ |
| ΔC_p | **−370 cal/mol·K** | ITC, ΔH vs T, 20–37 °C |
| ΔS_CONF | −70 cal/K·mol | derived |

**Binding is enthalpy-driven** (ΔH −13.5, −TΔS +0.9 — the entropy term is *unfavourable*). This
matters: bH1's rewritten paratope switched to **entropy-driven** binding (findings 01). A heavily
rewritten CDRH3 may bind by a different thermodynamic route — not necessarily worse, but different.

**⚠ K_D discrepancy to resolve:** **0.058 nM** (pmid 28750054, SPR) vs **0.5 nM** (pmid 21526167, SPR
30 °C) — about 10× apart, both by SPR. Likely IgG-vs-Fab and/or format/temperature differences. **Decide
which baseline we claim against, and state the method.** Related: pmid 28750054's best triple mutant
reached **K_D 1.9 pM, ~30× better than parent**, via ~6× k_on and ~5× k_off.

## 5. Interface size and composition

| Value | What | Source |
|---|---|---|
| **1350 Å²** | BSA of the HER2–trastuzumab complex (crystal) | 37569408 |
| ~1675 Å² | same, from MD-refined models | 21600050 |
| **910 Å²** | buried **apolar** area (bH1/HER2: 988 Å²) | 21526167 |
| **~60% apolar** | hydrophobic composition of the interface | 21526167 |
| ΔΔA_np = 30 Å² (3%) | agreement between ΔC_p-derived and observed apolar area → **binding is essentially rigid-body** | 21526167 |

**No per-CDR or per-chain BSA breakdown exists anywhere in the corpus** — that remains open, and is what
the structural agent is computing.

## 6. HER2 epitope — my measurement corroborated, with a functional test

pmid **34708303** states the crystal structure shows **19 amino-acid sidechains** of domain IV forming
the epitope. **My structural analysis independently measured exactly 19 contact residues.**

Three-loop definition (two sources, both citing Cho 2003): **loop 1 = 557–561, loop 2 = 570–573,
loop 3 = 593–603.** Loops 1 and 3 are primarily **electrostatic**; loop 2 is mostly **hydrophobic**.

**The functional proof:** the quadruple HER2 mutant **`P557S/P571S/P572S/F573S` completely abolishes
trastuzumab binding** (FACS on CHO transfectants). All four of those residues are in my measured contact
set, and `PPF 571-573` is the central patch. Human vs cynomolgus domain IV differs at only 2 of the 19
epitope positions, both peripheral, and **does not affect binding**.

## 7. Discrepancies to carry forward

- **Arg58 vs Arg59** (heavy chain) — pmid 21526167 writes `Arg50A+Arg58A`, pmid 23955267 writes
  `Arg50, Arg59`; US-8840896-B2 names both. Almost certainly the same Arg under Kabat vs sequential
  numbering, but **no source states the equivalence** — do not silently merge.
- **HER2 561 = Gln** (23955267's peptide sequence) **vs Asp** (25192037).
- **BSA 1350 Å² vs ~1675 Å²** — 25% apart, different methods.
- **D586 as an interface residue** (31891584) vs **"HER2 residues 581–590 are not visible in either
  3be1 or 1n8z"** (28750054) — so that contact is model-derived, not observed.
- **pmid 37569408's PRODIGY ΔG pair (−10.3 / −12.0 kcal/mol) is internally inconsistent with its own
  prose** and cannot be assigned as written. Do not use.
- **pmid 35808792 reports HER2 H-bond partners Asn89, Glu216, Tyr252, Cys289, Glu330** — domains I–III,
  irreconcilable with the crystallographic domain-IV epitope. **Docking artefact; discard.**

## 8. Still missing

- **pmid 12610629 (the 2.5 Å primary structure) is NOT in the corpus**, and **1N8Z's resolution is never
  stated** in any record we hold. The 2.5 Å figure comes from the structure's own PDB header.
- **No crystallographic H-bond/salt-bridge table or distance list** for the interface. The only explicit
  pairings (pmid 23955267) are **MD on a truncated 18-mer peptide**, not the crystal contact list.
- **No per-CDR BSA**, and **no per-residue ΔΔG for the paratope** — the Ala scan gives binned
  fold-changes, not ΔΔG.
- **No kink/bulge call or canonical-class assignment for trastuzumab's CDRH3.** Every corpus paper
  discussing H3 canonical clusters was checked individually and **none mentions trastuzumab, Herceptin,
  4D5 or 1N8Z**.
- The 19-residue epitope list in pmid 34708303 exists only in **Fig. 4A**, an image — our measured list
  is the usable version.
- **pmid 42968171 / PDB 7PKL is absent** from the corpus entirely.

## 9. ⚠⚠ THE SCORING RISK NOBODY FLAGGED: SPR may not be able to see our improvement

The competition measures binding **by SPR**. Two facts from the corpus make that a problem.

**1. Trastuzumab's off-rate is already near the instrument floor.**
Measured k_off for trastuzumab IgG1 is **4.4–5.8 × 10⁻⁵ s⁻¹** (pmid 31696314, six lots).
A Biacore T200's practical k_off limit is **~10⁻⁵ s⁻¹**. So the parent sits **within ~5× of the floor.**

**2. This has already defeated a published attempt on this exact molecule.**
pmid **31891584** engineered affinity-improved Herceptin variants and reported that they
*"were not distinguished via SPR from WT-like binders … these antibodies exhibit very slow off rates
that approach the koff limit (10⁻⁵ s⁻¹) of the SPR instrument used here (Biacore T200)."*
Their improved variants only showed up in **cell binding by flow cytometry**: apparent 0.02–0.04 nM
versus WT 0.3–0.9 nM on MCF7, i.e. ~10–20× tighter — invisible to SPR.

**Consequence for ranking:** a design that genuinely improves affinity by a k_off mechanism may score
as *indistinguishable from trastuzumab* on the SPR readout. **Prefer designs that improve k_on**, which
is measurable, over designs that improve k_off, which may not be. For reference, pmid 28750054's best
triple mutant got **~6× k_on and ~5× k_off**, and pmid 31891584 found **k_on varied only ~5× across all
20 variants** (6.7 × 10⁶ to 2.7 × 10⁵ M⁻¹s⁻¹) — so k_on is the harder thing to move but the thing that
will register.

**Also: measuring a Fab or scFv helps.** pmid 28750054 beat the ceiling by measuring **Fabs** (faster
k_off regime) with ≥900 s dissociation. Our designs are fragments, which works in our favour.

### Cross-assay spread on trastuzumab itself: ~30-fold

| K_D | Format / method | pmid |
|---|---|---|
| 58 pM | Fab, SPR 25 °C, HER2 immobilised | 28750054 |
| 63–91.6 pM | IgG1, SPR, anti-IgG1 capture | 31696314 |
| 0.283 nM | IgG, SPRi, antibody immobilised | 41538393 |
| 0.30 nM | IgG, BLI, Fc capture | 29593727 |
| 0.5 nM | Fab, SPR 30 °C | 21526167 |
| 1.8 nM / 0.16 nM | IgG, SPR vs SEC-SPR | 30328504 |

**A variant must beat the parent measured in the same run and same format.** Any "beats trastuzumab"
claim citing a literature K_D from a different assay is worthless — and this also resolves the
0.058 vs 0.5 nM discrepancy noted earlier: different formats, temperatures and orientations.

## 10. SPR protocol the corpus actually supports

- **Capture the antibody on an anti-Fc surface and flow HER2 ECD as analyte.** Both pmid 31891584 and
  31696314 do this. It gives genuinely monovalent affinity and avoids avidity inflation.
- **Do not covalently immobilise the IgG.** Random amine coupling can hit paratope lysines — pmid
  31891584 demonstrates exactly this failure chemically: variant `12-9-DM1` lost HER2 binding because
  amine/drug coupling hit the two engineered binding-site lysines. **Direct warning against introducing
  Lys into our CDRH3 if any amine conjugation is planned.**
- Buffer used identically by 28750054 and 31891584: **PBS + 3.4 mM EDTA + 0.05% Tween 20** (the EDTA is
  non-obvious; keep it for cross-comparability with the ADAPT numbers).
- **Single-cycle kinetics, 180 s association / 1800 s dissociation**, keep analyte R_max < 100 RU,
  double-reference, 1:1 fit. For sub-100 pM the long dissociation is mandatory.
- **Run SEC on every variant before SPR** to exclude aggregation artefacts (31891584 did).
- **Re-measure the parent in every batch and report fold-change against that same-batch parent.**
  pmid 28750054 shows this choice alone moves the headline numbers (104/46/32 vs 96/55/31).
- Avoid basing an affinity claim on a **two-state fit** — pmid 42909185 concedes its micro-rate
  constants are model-dependent.
- **Corroborate with an orthogonal readout**: cell binding by flow cytometry on matched high/low-HER2
  lines (SKOV-3 high, MCF-7 low), plus ITC or BLI.

## 11. The ADAPT numbering question — independently resolved, same answer

Two mining passes independently reconciled the two papers' numbering and agree:

- pmid 28750054 writes **`HD102F`** (sequential VH numbering).
- pmid 31891584, the same group's follow-up on the same campaign, writes **`H-D98F`** (Kabat).
- **These are the same CDR-H3 aspartate** — our **D102** in 1N8Z numbering.
- Decisive check: **there is no Asp at Kabat H102 in trastuzumab (it is Tyr)**, so `HD102F` cannot be
  Kabat. **Neither paper states its numbering scheme.**

**Canonical ADAPT affinity sites (Kabat):** VH **D31**, VH **D98**; VL **D28**, VL **S50**.
Full position list explored (pmid 31891584): heavy **D31, Y57, D98, G99**; light **D28, N30, F53, Y92**.
Note the substitution sets differ between papers (28750054: D31R, D28R/K; 31891584: D31K/M, D28K only).

**The best triple mutant's exact composition is in Fig 1C and is NOT in the retrieved text. Do not
guess it.** What is known: it is `HD102F` + a second single + (`LD28R` or `LD28K`), reached
**K_D 1.9 pM (~30× better than parent)**, and the best *double* (`HD102F` + `HD31R`) was
**"heterogeneous upon purification"** and was excluded — a developability warning about that pairing.

## 12. One contradiction between mining passes, resolved

The second pass reported that all 38 patent records have **empty `fulltext`** and concluded patents
contributed nothing. The first pass extracted extensive verbatim text from **US-8840896-B2**, including
the alanine-scan data. **Both are right:** patent text lives in the **`claims` and `description`**
fields, not `fulltext`. The second pass only checked `fulltext`.

**So the alanine-scan data in §2 stands** — it comes from `US-8840896-B2`'s description field, citing
Kelley & O'Connell, *Biochemistry* 32:6828-35 (1993). Anyone re-mining this corpus must read
`claims`/`description` for `AMPC_*` records.
