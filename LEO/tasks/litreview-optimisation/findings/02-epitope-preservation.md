# Epitope preservation — trastuzumab paratope / HER2 domain IV

Synthesised from abstracts only (`cache/biomedcore__trastuzumab-her2-extracellular-domain-iv-epitope-structure.json`,
`...-fab-her2-crystal-structure-contact-residues-1n8z.json`, `...-her2-epitope-mapping-competitive-binding-trastuzumab-pertuzu.json`).
**Not yet verified against full text.** Numbers marked "conference abstract" are not peer-reviewed.

## Headline for design decisions

1. **A VHH with no light chain and none of trastuzumab's CDRs binds trastuzumab's epitope at
   sub-nM affinity** (pmid 37632647). Epitope identity is a property of the docked surface, not of
   the CDRH3 sequence. The <70% CDRH3 identity rule is therefore **not inherently
   epitope-breaking** — this is the single most important finding for our track.
2. **"Binds HER2 domain IV" is not evidence of epitope preservation.** Domain IV contains at least
   two non-overlapping epitopes (HLX22 vs trastuzumab, pmid 38982548) plus a membrane-proximal one.
   A near-neighbour epitope can fail to block trastuzumab entirely (AACR abstract,
   doi 10.1158/1538-7445.am2026-4553). Only competition against trastuzumab counts.
3. **The epitope is discontinuous, spanning ~50 aa** of domain IV (SITC abstract,
   doi 10.1136/jitc-2024-sitc2024.0269). A ~50 aa discontinuous footprint cannot be covered by
   CDRH3 alone, so leaving the light chain and CDRH1/H2 untouched should anchor the docking
   register while CDRH3 is rewritten. *(Inference from geometry, not a measurement.)*
4. **Binding is electrostatics-dominant**, hydrophobics secondary (pmid 23955267). Preserve the
   net charge/polar character of the new CDRH3 face; hydrophobic-heavy rewrites are riskier.

## Interface map (what the abstracts actually establish)

| Claim | Evidence | pmid / doi | Type |
|---|---|---|---|
| Trastuzumab Fab binds the juxtamembrane region of HER2 domain IV | X-ray 2.5 Å, human HER2 ECD–Herceptin Fab | 12610629 | measured |
| Epitope discontinuous, ~50 aa | contrast with H2CasMab-2's short continuous epitope | doi 10.1136/jitc-2024-sitc2024.0269 | measured (conf. abstract) |
| Three HER2 contact loops: **991–995**, **1004–1007**, **1027–1037**; loops 1–2 dominate; Lys1027 least favourable; total ΔG −68.55 kcal/mol | AMBER 14 MD + MM-GBSA on 1N8Z | doi 10.1096/fasebj.31.1_supplement.775.23 | computational (conf. abstract) |
| Trastuzumab Ser50 (+0.116) and Asn30 (−0.999 kcal/mol) are weak/unfavourable contributors — **chain and CDR not stated** | same MM-GBSA decomposition | doi 10.1096/fasebj.31.1_supplement.775.23 | computational (conf. abstract) |
| Trastuzumab binding does **not** block the HER2–HER2 dimerisation interface | MD of HER2 dimers ± trastuzumab | 33557368 | computational |
| Trastuzumab IgG binds HER2 essentially **monovalently** (1:1, steric crowding); pertuzumab 2:1 | DLS + SEC-tetra detection | 37569408 | measured |
| No cooperativity between trastuzumab and pertuzumab sites | cryo-EM 4.36 Å ternary complex | 31042744 | measured |
| Better modern coordinates exist: HER2–m66–trastuzumab **3.2 Å**, HER2–r40–trastuzumab–pertuzumab **3.1 Å** | cryo-EM | 41538393 | measured |
| Interface is reproducible by solution covalent-labelling MS | FFAP footprinting + MS | 38698660 | measured |
| Trastuzumab's epitope lies in an evolutionarily conserved HER2 motif | PROSITE/PRINTS motif prediction | 25192037 | computational |

## Epitope-shift risk rules

**High risk**
- CDRH3 rewrite that changes **loop length or apex geometry** so it can no longer reach HER2 loops
  991–995 and 1004–1007 (the energetically dominant segments).
- **Flipping the electrostatic character** of the CDRH3 face — can still bind HER2 but re-dock onto
  an adjacent hydrophobic surface.
- Drifting toward the **HER2 dimerisation face** (the pertuzumab-like failure mode).
- Treating "good affinity for domain IV" as success.

**Lower risk**
- Large CDRH3 **sequence** change at constant length, backbone conformation class and net charge.
- Mutating the weak contributors (Ser50, Asn30 — chain unassigned, verify first).
- Leaving light chain + CDRH1/H2 intact as the docking anchor.
- Note: **avidity cannot rescue a weakened paratope** — trastuzumab engages HER2 monovalently.

**Cheap negative control:** trastuzumab's epitope requires properly folded, non-truncated domain IV.
A variant whose binding **survives DTT-induced misfolding or detects p95** has likely shifted epitope
(doi 10.1136/jitc-2024-sitc2024.0269).

**Competition assays are confounded.** Excess pertuzumab — a confirmed *non-overlapping* epitope —
*inhibited* radiolabelled trastuzumab binding in a cell assay (pmid 34757737). A "fails competition"
verdict should be confirmed structurally or by footprinting before discarding a design.

## Critical gaps (these block a per-CDR design budget)

- **No trastuzumab paratope residue list by CDR anywhere in the cache.** Not one record names a
  trastuzumab CDRH3 residue. No buried-surface-area or H-bond inventory per CDR.
- **No per-CDR or per-chain energy decomposition for trastuzumab.** No alanine scan, no CDR swap,
  no chain-shuffling. So "how much binding energy does CDRH3 carry?" is **unanswered**.
- **No experimental limit on how much CDRH3 can change** while retaining the epitope. Nothing
  supports or refutes the <70% threshold empirically.
- HER2 epitope residue numbers rest on **one un-peer-reviewed computational conference abstract**.
- **1N8Z is never confirmed** as the HER2–trastuzumab complex by a structural paper's own abstract —
  only by downstream docking papers citing it as input. Verify against pmid 12610629.

**How to close these:** read pmid 12610629 in full; recompute the contact set and a per-CDR energy
budget on the 3.1–3.2 Å cryo-EM coordinates (pmid 41538393) or on 1N8Z with FoldX/MM-GBSA/Rosetta.

## Read in full, in order

| pmid | Why |
|---|---|
| 12610629 | Primary 2.5 Å structure; the contact-residue list and per-CDR footprint are in the paper |
| 41538393 | 3.2/3.1 Å cryo-EM with trastuzumab — best coordinates for recomputing contacts |
| 37632647 | VHH at trastuzumab's epitope: existence proof + reusable epitope-mapping protocol |
| 38982548 | Two non-overlapping epitopes within subdomain IV — the exact failure mode to catch |
| doi 10.1096/fasebj.31.1_supplement.775.23 | Only per-residue energy decomposition; conference abstract, verify everything |
| 23955267 | MD/QM/MM-GBSA interface energies; establishes electrostatics dominance |
| 38698660 | FFAP footprinting — orthogonal epitope-preservation measurement |
| 34757737 | Competition-assay confounders before trusting any competition score |

No retracted records. Five load-bearing claims come from conference abstracts with no pmid.
