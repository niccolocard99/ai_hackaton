# Final ranked list, optimisation track

Written 2026-11-14 14:59 UTC. Every design passed the real abguard optimisation checker (assigned `novel`) and fails the de novo script, as expected: a design is never reclassified.

## What to expect, measured and not hoped

Held-out backtest (`backtest.py`, 20 repeats of 5-fold cluster-grouped CV on the 179 measured binders, test variants at least 3 mutations from anything trained on, which is the situation of a legal design). Inside the legal band (4 or more substitutions from trastuzumab) the ranker's top pick beat trastuzumab's measured 0.152 nM **1.0%** of the time, landed within 3x **14%** of the time, and the median measured KD of its top picks was about 9 nM. Binding is a long shot for every design. That is why the lead four are four bets on different axes and not the four best scores.

## The lead four

| # | design | add-ons | why |
|---|---|---|---|
| 1 | `WGVFGFYEMS` | D28K+fwTm | Literature D102F loop (3 other substitutions only) + VL D28K + two framework Tm nudges. Every lever with outside evidence in one molecule: the S-tier aspirant and the highest-variance design. |
| 2 | `WWVEGFYEMI` | D28K | The ranker's best bootstrap-robust CDRH3 with the light-chain affinity boost. Tests whether the literature boost stacks on the data-driven loop. |
| 3 | `WWSDGFYEMY` | fwTm | Epitope-conservative loop (D102 kept, contact risk 0.62) plus the framework nudge: the safest bet on the epitope axis, with a second axis from Tm. |
| 4 | `WKSQGFYEMI` | cdrh3only | A ranker-best loop with nothing added. If slots 1 and 2 beat it, the add-ons are doing the work. |

## Full list

| # | slot | CDRH3 | add-ons | substitutions | pred | p10 | epi risk | VH dTm C | H3 A, exact fold |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 1 all three levers | `WGVFGFYEMS` | D28K+fwTm | G101V/D102F/A106E/D108S | -0.02 | 0.92 | 1.00+L1 | +3.4 |  |
| 2 | 2 data loop + D28K | `WWVEGFYEMI` | D28K | G100W/G101V/D102E/A106E/D108I | -0.60 | 0.89 | 1.38+L1 | -0.3 |  |
| 3 | 3 epitope + Tm | `WWSDGFYEMY` | fwTm | G100W/G101S/A106E/D108Y | -0.02 | 0.81 | 0.62 | +2.9 |  |
| 4 | 4 control | `WKSQGFYEMI` | cdrh3only | G100K/G101S/D102Q/A106E/D108I | -0.43 | 0.90 | 1.38 | +0.1 |  |
| 5 | fill (d102c) | `WTVFGFYEME` | D28K | G100T/G101V/D102F/A106E/D108E | -0.02 | 0.85 | 1.38+L1 | +1.6 |  |
| 6 | fill (d102c) | `WGSFGFYEME` | fwTm | G101S/D102F/A106E/D108E | -0.01 | 0.96 | 1.00 | +4.2 |  |
| 7 | fill (data) | `WYSEGFYEMY` | cdrh3only | G100Y/G101S/D102E/A106E/D108Y | -0.34 | 0.86 | 1.38 | +0.5 |  |
| 8 | fill (epi) | `WKSDGFYHMI` | fwTm | G100K/G101S/A106H/D108I | +0.17 | 0.84 | 0.62 | +2.5 |  |
| 9 | fill (d102a) | `WWVFGFYEMI` | D28K+fwTm | G100W/G101V/D102F/A106E/D108I | -0.28 | 0.91 | 1.38+L1 | +2.1 |  |
| 10 | fill (data) | `WWVEGFYEMI` | cdrh3only | G100W/G101V/D102E/A106E/D108I | -0.60 | 0.89 | 1.38 | -0.3 |  |
| 11 | fill (d102c) | `WGVFGFYEMS` | cdrh3only | G101V/D102F/A106E/D108S | -0.02 | 0.92 | 1.00 | +0.8 |  |
| 12 | fill (epi) | `WWVDGFYNMI` | cdrh3only | G100W/G101V/A106N/D108I | +0.19 | 0.81 | 0.62 | -0.8 |  |
| 13 | fill (d102a) | `WGSFYFYEMI` | D28K | G101S/D102F/G103Y/A106E/D108I | -0.02 | 0.90 | 1.88+L1 | +1.3 |  |
| 14 | fill (data) | `WGSEAFYEMI` | fwTm | G101S/D102E/G103A/A106E/D108I | -0.34 | 0.89 | 1.88 | +4.2 |  |
| 15 | fill (d102a) | `WYSFGFYEMT` | cdrh3only | G100Y/G101S/D102F/A106E/D108T | +0.03 | 0.94 | 1.38 | +0.2 |  |
| 16 | fill (data) | `WKVYGFYEMY` | D28K | G100K/G101V/D102Y/A106E/D108Y | -0.26 | 0.81 | 1.38+L1 | -0.3 |  |
| 17 | fill (d102c) | `WTVFGFYEME` | cdrh3only | G100T/G101V/D102F/A106E/D108E | -0.02 | 0.85 | 1.38 | +1.6 |  |
| 18 | fill (epi) | `WWSDGFYEMY` | cdrh3only | G100W/G101S/A106E/D108Y | -0.02 | 0.81 | 0.62 | +0.3 |  |
| 19 | fill (data) | `WKSQGFYEMI` | D28K | G100K/G101S/D102Q/A106E/D108I | -0.43 | 0.90 | 1.38+L1 | +0.1 |  |
| 20 | fill (d102c) | `WGSFGFYEME` | D28K | G101S/D102F/A106E/D108E | -0.01 | 0.96 | 1.00+L1 | +1.8 |  |

## Read before trusting any number

- `pred` is the ranker bootstrap mean of log10 KD plus one standard deviation. It is a ranking, not a KD, and it was fitted on the trastuzumab light chain and framework, so it says nothing about any add-on.
- The bootstrap lower bound and an anchored estimate did **not** beat the plain ranker on held-out variants (Spearman +0.668 and +0.673 against +0.666). They are a conservative choice, not a validated gain.
- VL D28K and D102F/W come from the literature, not from our assay; the size of D28K alone is not in the retrieved text. D28 contacts HER2, so those designs change a contacting residue (`+L1`).
- VH dTm is NanoMelt, nanobody-trained and out of domain for a VH, error 4.1 C. A light-chain change is not visible to it. Framework gains of +1 to +3 C are inside that error.
- Boltz-2 does not rank binders (every metric at chance on measured data, EXP-004 and EXP-006). The base H3 column, when filled, is a pose-plausibility flag for the CDRH3 loop only, folded from the exact final sequence.
- Thermostability and epitope preservation are not predicted by the ranker. No developability scorer is installed on the VM.
