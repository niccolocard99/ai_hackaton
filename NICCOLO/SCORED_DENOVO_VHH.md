# Scoring the 12 de novo VHH designs

Date 2026-11-14 · Source `~/hackathon-team-3/LEO/tasks/denovo-vhh-epitope/` on the VM · Script
`eda/p21_score_leo_final.py` sibling analysis, output `eda/data/denovo_scored.csv`

**These are genuinely de novo.** All 12 have identity to trastuzumab 0.630–0.683, i.e. under the 0.70 gate;
`track_assigned = de_novo`; all PASS base validation. 127-residue VHHs, so a separate prize pot from the
optimisation list.

## Stage 1 does not apply here, and must not be used

The frozen pre-filter reads trastuzumab's 10-residue CDRH3 window on trastuzumab's paired framework. A VHH
is a different molecule — single domain, no light chain, different framework — and `check_parent` refuses it
by design. **There is no binding predictor for these designs.** That is the central fact about this list.

What *does* apply, and it is the mirror image of the optimisation track:

| Tool | Optimisation track | De novo VHH track |
|---|---|---|
| sequence pre-filter | the main signal (Spearman +0.69) | **inapplicable** |
| NanoMelt Tm | out of domain (VH, not a nanobody) | **in domain** (639 measured nanobody Tm) |
| fold confidence (ipTM) | measured at chance (AUROC 0.50) | the regime where it has signal (AF3 0.75) |
| epitope coverage | no handle at all | **measurable** against a defined patch |

## The scores

| rank | design | identity | epitope_cov | contacts_in_ep | NanoMelt Tm °C | Boltz-2 ipTM | min_ipae |
|---|---|---|---|---|---|---|---|
| 1 | fixedfw_89 | 0.630 | 0.68 | 0.87 | 69.7 | 0.575 | 7.06 |
| 2 | fixedfw_6 | 0.646 | 0.58 | 0.65 | **73.5** | **0.961** | **0.44** |
| 3 | fixedfw_13 | 0.638 | 0.47 | 0.60 | 69.0 | 0.664 | 8.71 |
| 4 | fixedfw_5 | 0.654 | 0.47 | 0.69 | 70.5 | 0.546 | 8.41 |
| 5 | fixedfw_48 | 0.648 | **0.79** | **0.94** | 65.6 | 0.437 | 15.03 |
| 6 | freefw_0 | 0.683 | 0.63 | 0.92 | 69.7 | **0.973** | **0.37** |
| 7 | fixedfw_42 | 0.651 | 0.63 | 0.92 | 70.5 | 0.673 | 8.55 |
| 8 | freefw_1 | 0.659 | 0.58 | 0.92 | 71.6 | 0.546 | 11.31 |
| 9–12 | fixedfw_0, _7, _3, _39 | 0.640–0.659 | 0.37–0.53 | 0.78–0.91 | 67.7–70.4 | not yet folded | — |

NanoMelt was run CPU-only (`CUDA_VISIBLE_DEVICES=""`) to avoid adding a third job to a GPU that already had
two Boltz-2 containers on it.

## Thermostability: real discrimination, no usable baseline

Tm spans **65.6 to 73.5 °C**, a 7.9 °C spread against NanoMelt's 4.1 °C MAE — so unlike the CDRH3 case,
where the whole spread sat inside the error bar, there is genuine signal here. `fixedfw_6` (73.5 °C) and
`freefw_1` (71.6 °C) are the stable end; `fixedfw_48` (65.6 °C) the weak end.

**The catch:** the tier score is improvement *over trastuzumab*, and trastuzumab is not a nanobody, so
NanoMelt cannot give a comparable baseline. We can rank these 12 against each other on Tm. We cannot say any
of them beats trastuzumab.

## Binding: Boltz-2 cannot endorse any of them

Only **2 of 8** folded designs reach ipTM ≥ 0.7. Put that against the published VHH–antigen benchmark
(`eda.py`, Smorodina et al.), which is exactly this regime:

- Boltz-2 gives ipTM **≥ 0.7 to 94% of shuffled, known non-binding** VHH–antigen pairs.
- Its real-vs-shuffled AUROC is **0.64** — the weakest of the three models (AF3 0.75, Chai-1 0.57).

So **6 of 8 designs score below the level Boltz-2 hands most known non-binders.** The two that score high
(`fixedfw_6` 0.961 / min_ipae 0.44, `freefw_0` 0.973 / 0.368) merely reach the level of a typical
non-binder; Boltz-2 cannot distinguish them from one. This is a weak-negative read on the batch, not an
endorsement of the top two.

**Recommendation: refold with Protenix v2 or OpenFold3, not Boltz-2.** Boltz-2 is the weakest model for this
task, and in my own calibration it had almost no dynamic range, while Protenix showed a real spread. Both
are on the VM, both take templates, and Protenix costs ~12.8 s/structure — 12 designs is about 3 minutes.

## Epitope coverage: the de novo track's one real advantage

4 of 12 clear both `epitope_cov ≥ 0.6` and `contacts_in_ep ≥ 0.7`: **fixedfw_89, fixedfw_48, fixedfw_42,
freefw_0**. `fixedfw_48` is the best on this axis (0.79 / 0.94) while being the worst on Tm (65.6 °C).

This matters strategically: the optimisation list has **no epitope handle at all**, so a de novo VHH is the
only design that can make a two-property claim backed by anything.

## Verdict

**12 valid de novo candidates.** On the three judged properties:

- **Binding: unpredicted.** No applicable model. Boltz-2 ipTM leans weakly negative on 6 of 8.
- **Thermostability: ranked, not baselined.** Real spread, in-domain predictor, but no trastuzumab comparison.
- **Epitope: 4 designs with a defensible claim** — the only epitope evidence anywhere in our submission.

Best single candidates on current evidence: **`freefw_0`** (epitope 0.63/0.92, the highest ipTM, median Tm)
and **`fixedfw_6`** (highest Tm 73.5 °C, second-highest ipTM, weaker epitope). **`fixedfw_48`** is the
epitope-purist pick at the cost of Tm.

The de novo track is worth entering precisely because it is the only place we have epitope evidence — and
because an empty track scores zero while the optimisation track is, per LEO's own backtest, a 1% shot on
binding.
