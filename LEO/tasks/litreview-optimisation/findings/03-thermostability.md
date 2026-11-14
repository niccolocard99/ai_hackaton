# Thermostability, stabilising mutations, affinity–stability trade-off

From the five `cache/biomedcore__*thermostability*/*stability*/*disulfi*` files. **Abstracts only.**

## Format arithmetic first *(inference, not abstract-stated)*

| Format | Approx. length | Verdict vs ≤250 aa cap |
|---|---|---|
| **VHH** | ~120–130 aa | Comfortable headroom |
| **scFv** | ~245–255 aa (VH 120 + linker 15–20 + VL 107) | **Right at the cap — tight** |
| **Fab** | ~435–445 aa | **Cannot fit** |

So the real choice is scFv (tight) or VHH (comfortable). Anything that *adds* residues (tags, longer
linkers) is budget-hostile for scFv but fine for VHH.

## The uncomfortable headline

**Trastuzumab's own Fv (clone 4D5) is used in this literature as the gold-standard *stable* scFv
framework.** Other groups graft their CDRs *onto* it to rescue fragments that failed from
"insufficient thermal stability" (pmid **10582696**), and it is a reference system for stability
engineering (pmid **10393549**, which calls such frameworks "very suitable recipients in CDR grafting").

**Consequence:** beating trastuzumab on Tm will not come from finding a better framework. It has to
come from (a) **format** — a VHH removes the VH/VL interface as a failure mode entirely, (b) an
**added engineered disulfide**, and (c) **not wrecking the loops**. *(Inference.)*

## Measured Tm numbers available (and the baseline gap)

| What | Number | pmid |
|---|---|---|
| Herceptin whole IgG, far-UV CD | **~75 °C** | 31071101 |
| Trastuzumab/pertuzumab bispecific MBS301 | Tm1 **72.0 °C**, Tm2 **84.8 °C**; trastuzumab's own values only "high and similar" | 30081724 |
| Trastuzumab + biosimilars, **nanoDSF, domain-resolved Tonset/Tm** | **no numbers in abstract** — most likely source of a real Fab baseline | **32380351** |
| HER2 nanobody NB46 (competes with trastuzumab) | retained HER2 affinity **up to 90 °C** | 42911746 |
| DSF vs DSC offset | DSF on average **3 °C lower** than DSC | 25385011 |

**The biggest hole in the whole review: no trastuzumab Tm for any *fragment* format.** The ~75 °C
IgG figure is a whole-molecule value and is **not** a valid scFv/Fab target. Pull pmid 32380351 full
text for the Fab transition.

## Stabilising toolkit, with numbers

| Strategy | Region | Effect | Affinity | pmid |
|---|---|---|---|---|
| Framework point mutations, thermal-challenge screen | VH/VL framework (positions not in abstract) | **+14 °C** single; **+12 °C** more in combination | "retain **full** antigen-binding affinity" | **20457695** |
| Framework replacement | **VH C–C′ loop** | **+19.6 °C**, +1.9 M chemical, 27× yield, stable at pH 2.0 | not reported | **32326443** |
| Rosetta-designed substitution *clusters* | two 2-residue clusters | higher Tm; full binding after **2 h at 70 °C**; 1.6 Å structure matched model | full binding retained | 32336757 |
| Consensus + H-bonds + surface hydrophilicity | 14 positions | **75% of mutations stabilising, none destabilising**; **+9.1 kcal/mol** combined | not reported | 16926023 |
| **Interdomain disulfide VH44–VL100** | VH44/VL100 | GdnHCl₅₀ **2.4 M → 4.2 M**; monomer-only by SEC, no aggregates | "**without affecting affinity**" | **21339972** |
| Intradomain Ig-fold disulfides (CH3 benchmark) | CH3 | **+10 °C**, **+5 °C**, **+15 °C** combined; **+7.5/+15.5/+19 °C** in a HER2-binding CH3 | FcRn binding unaffected | 22272277 |
| Mutual VH/VL stabilisation | VH **K66R**, **N52S** | "additive and transferable"; stabilising the partner domain stabilises the other | not reported | 9748318 |
| Framework mutations port across formats | VH/VL FW | same mutations improve scFv, Fab **and** IgG | not reported | 22763265 |
| Trastuzumab-specific FW design, **distal to CDRs** | trastuzumab FW (positions not in abstract) | stabilising mutants designed by MD | **preserved HER2 engagement and ADCC** — but one LC FW mutation **lost Fab- and Fc-mediated effector function** while retaining binding | **40653772** |

**Prefer clusters/combinations** — additivity is reported repeatedly (9748318, 16926023, 20457695, 22272277).

## The trade-off: our exact failure mode is documented

**pmid 36153698** CDR-grafted a VHH onto a different framework: "significant decrease in affinity
**and** thermal stability and a large conformational change in CDR3". Cause: loss of **bulky residues
near CDR3** and of **FR2↔CDR3 intramolecular contacts**. Framework **back-mutations restored both**.
**This is the mechanistic template for "my CDRH3 is rewritten — where do I buy stability back":
FR2 and the FW residues packing against CDR3, not generic consensus positions.**

Supporting evidence that loops, not framework, set fragment Tm:
- **NbThermo** (564 nanobodies with measured Tm, pmid **37042467**): "**no apparent differences in
  sequence patterns between the frameworks** of Nbs with lower and higher melting temperatures,
  indicating that the **highly variable loops** play a relevant role". *So for a VHH, a rewritten
  CDR3 is the dominant Tm variable and framework tinkering may not rescue it.*
- **AbMelt** (pmid **38851888**): CDR radius-of-gyration deviation at 400 K is the **highest-correlated
  descriptor with aggregation temperature** (rp = −0.68 ± 0.23); internal-contact deviation correlates
  with Tm (−0.69 ± 0.03). A floppier rewritten CDRH3 should be expected to drop both Tagg and Tm.
- **Aggregation moves independently of Tm** (pmid **21732420**): VH domains differing only in CDR
  sequence had "dramatically different aggregation propensities **in a manner not correlated with
  their conformational stabilities**", hotspots in/near CDR1. **A design can pass on Tm and still fail.**

**How big is the penalty? Unknown.** Honest answer: **no CDR-mutation paper in this cache reports a
ΔTm number.** Every hard ΔTm comes from the framework/disulfide side. *Inference:* the stabilising
toolkit's headroom (≈+5 to +19 °C) is comfortably larger than typical loop penalties, so the
trade-off is usually winnable — but **budget for an unquantified loss and over-engineer the
framework rather than aim for break-even.**

## Three results that break the trade-off framing

1. **CDR mutations can give affinity AND stability together** (pmid **23995618**): the authors
   "eschewed framework-based stabilizing mutations in favor of CDR mutagenesis" with a thermal-challenge
   selection, and got **>18-fold affinity and improved thermal stability from only 4 mutations in
   VL-CDR3**; final molecule formulated at 100 mg/ml, <2% HMW after 7 weeks. **CDR-L3 is a cheap place
   to recover both properties** — and it is not part of the epitope contact set we must preserve.
2. **A pre-stabilised CDR-H3 is what *enables* affinity optimisation** (pmid **41613601**): a **twin
   cysteine (paired)** motif in CDR-H3 "minimized structural perturbations caused by affinity-enhancing
   mutations", and "antibodies **lacking** H3 CDR stabilization showed **no affinity improvement**
   after the same optimization". For a team deliberately rewriting CDRH3, build rigidity into the new
   loop *before* optimising binding.
3. **4–6 CDR mutations is the sweet spot, with mixed outcomes per mutation** (pmid **28928732**):
   "some natural diversity mutations **simultaneously enhance all three** key properties (affinity,
   specificity, stability) while other mutations display trade-offs".

**Contested:** does maturation rigidify CDRs? 26132144 (VH rigidifies, VL/CDR-L2 become *more*
flexible) and 21287614 say yes; 29545810 (thousands of models + PDB + MD) finds only a "slight"
decrease, "not as drastic as previously reported". **Don't over-weight a rigidity objective.**

## Engineered disulfides: positions and the paired-Cys caveat

**Positions with support:** **VH44–VL100** is the industry standard, measured GdnHCl₅₀ 2.4→4.2 M with
unchanged affinity (21339972). Conserved VH/VL framework positions generally (8356052). Additional
**intradomain** disulfides in VH or VL framework, which "did not reduce antigen-binding affinity"
(25038323). **Paired Cys inside CDR-H3** (41613601).

**Three documented failure modes:**
1. **Never delete the canonical intradomain disulfides** (VH 22/96, VL 23/88). They were "**essential
   for the stability of all antigen binding fragments investigated**"; reduction gave irreversible
   denaturation (pmid **1736986**). Floor is 4 paired Cys in an scFv, 2 in a VHH; an engineered pair
   takes that to 6 (or 4), all paired.
2. **An odd cysteine does its damage in a CDR.** pmid **24598312**: a free Cys-VL85 mispaired with
   Cys-VL16, "linking light chain FR1 and CDR3", with "deleterious negative effects … likely due to
   **spatial hindrance in the light chain CDR3**"; Cys→Ser gave +30% protective efficacy in vivo.
   **This is the cleanest experimental justification for the no-unpaired-cysteine rule.**
3. **Even a correct pair can mispair *between* molecules** (pmid **41996728**): an "unexpected
   interchain disulfide bond linking two scdsFv segments". Prefer a **buried intradomain** pair over a
   solvent-exposed interface pair if aggregation is scored.

If a Cys must be removed, **Val and/or Ala pairs were preferred** for replacement (25038323).

## Predictor accuracy is coarser than the effects we design

TEMPRO (nanobody Tm): **MAE 4.03 °C, RMSE 5.66 °C** (pmid 39154093). AbMelt: **R² = 0.60 ± 0.06**
(38851888). pmid 30192891 warns of "the **limited applicability domain** of each model".
**A predicted +3 °C is inside the noise.** Stack several independently-supported stabilising changes;
do not rely on a predictor to certify that we beat trastuzumab.

## Read in full

| pmid | Why |
|---|---|
| **40653772** | Only trastuzumab-specific stabilising-FW paper — need the actual positions |
| **36153698** | Our exact failure mode in a VHH; FR2↔CDR3 levers that buy stability back |
| **20457695** | Largest clean scFv FW ΔTm (+14/+12 °C) with full affinity retained — need positions |
| **21339972** | VH44–VL100 reference implementation |
| **23995618** | 4 CDR-L3 mutations → affinity *and* stability; the thermal-challenge selection to copy |
| **10393549** | Uses the 4D5 (trastuzumab) scFv; weak-part triage (VL vs VH vs interface) |
| **32380351** | nanoDSF on trastuzumab — **likeliest source of a real Fab Tm baseline** |
| **32326443** | +19.6 °C framework swap localised to the VH C–C′ loop, away from the epitope |

## Gaps

- **No trastuzumab fragment Tm baseline anywhere** — the single biggest hole.
- **No quantified ΔTm cost of CDR/CDRH3 mutation.** Nothing answers "what does a <70% identity
  CDRH3 rewrite cost in °C?" We are blind on the size of the penalty to cover.
- **No VHH framework stabilising set with numbers**; 37042467 says FW sequence does not predict Nb Tm.
- **Linker length/composition unquantified**; no abstract compares VH-VL vs VL-VH orientation for stability.
- **Exact positions missing for the three biggest ΔTm results** (20457695, 32326443, 32336757).
- **Nothing connects engineered-disulfide placement to epitope preservation.** VH44/VL100 sits near
  VL P44 / VH L45, flagged as "of particular importance for Fv association" (28947772).

Two on-point records are conference abstracts with no pmid and unverified retraction status
(DOI 10.1158/1557-3265.d32026-a020; DOI 10.1158/1538-7445.am2026-2759).
