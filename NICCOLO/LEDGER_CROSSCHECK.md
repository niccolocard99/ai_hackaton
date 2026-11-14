# Cross-check of LEO's experiment ledger against the full SPR-768 table

**Date:** 2026-11-14 · **Checked:** `LEO/tasks/experiment-ledger` (EXP-001…008, TH-001…019, REVIEW.md)
**Scripts:** `eda/p15_ledger_check.py`, `eda/p16_position_test.py`

The ledger's EXP-008 / TH-019 asked for exactly this check as **IDEA-017**: its affinity conclusions came
from the 60-variant fold-calibration subset I built, which was *chosen* (30 binders spanning the KD range,
30 non-binders), not sampled. I have the full table locally, so here is the answer on all **131 legal
quantified variants** instead of 21.

**Verdict: the existence proof survives, the rate is 6x worse, and the design recipe is partly inverted.**
IDEA-008 was next in the ledger's queue to be handed to a GPU. It should not be run as written.

---

## 1. The hit rate is 1 in 131, not 1 in 21

| Baseline for "beats trastuzumab" | Legal variants beating it |
|---|---|
| 0.152 nM (the single row EXP-008 used) | **1 of 131 (0.8%)** |
| 0.135 nM (median of the 8 replicate rows) | 1 of 131 (0.8%) |
| 0.029 nM (best replicate) | **0 of 131** |

EXP-008 inferred 1 of 21 = 4.8%. The true rate is **0.8%**, six times rarer. The winner is the same single
variant, `WGGSGPYLFD` (4 substitutions, 0.095 nM).

Two things the ledger should add:

- **The baseline needs stating explicitly.** Trastuzumab's 8 replicate rows span 0.029–0.178 nM, a 6x
  range. Whether *anything* in the legal band beats trastuzumab depends on which replicate you pick.
- **None of the 131 is submittable.** They are all *in* the provided HER2 data, so all 131 fail the novelty
  rule, not just `WGGSGPYLFD`. They are evidence about the landscape, never candidates.

## 2. The A106/M107/D102 convergence is a base-rate artefact

EXP-008 found A106 mutated in 7 of 8 and M107 in 8 of 8 of its tightest legal variants, and read that as a
pattern to design toward. But those positions are mutated in almost *every* legal variant:

| Position | Mutated in all 131 legal | Mutated in the 8 tightest |
|---|---|---|
| M107 | **92%** | 75% |
| D102 | 86% | 100% |
| A106 | 66% | 75% |

M107 is mutated *less* often among the tightest than at baseline. Conditioning on winners without comparing
to the base rate produced a pattern that carries no information. This is the standard trap: a feature common
among winners is only interesting if it is rarer among losers.

## 3. Worse — the direction is backwards for two of the three

Mann-Whitney on log10 KD, mutated vs conserved, all 131 legal variants, Benjamini-Hochberg corrected across
all 10 positions:

| Position | WT | n mut / con | Median KD mutated | Median KD conserved | Effect | q (BH) |
|---|---|---|---|---|---|---|
| **Y105** | Y | 56 / 75 | 167.6 nM | 48.0 nM | mutating **hurts** | **<0.001** |
| **W99** | W | 71 / 60 | 157.1 nM | 47.3 nM | mutating **hurts** | **<0.001** |
| **F104** | F | 87 / 44 | 139.0 nM | 47.8 nM | mutating **hurts** | **0.028** |
| **M107** | M | 120 / 11 | 128.8 nM | 17.0 nM | mutating **hurts** | **0.040** |
| **G100** | G | 55 / 76 | 66.0 nM | 142.1 nM | mutating **helps** | 0.062 |
| **G101** | G | 68 / 63 | 90.0 nM | 147.0 nM | mutating **helps** | 0.077 |
| G103 | G | 42 / 89 | 156.5 nM | 93.0 nM | mutating hurts | 0.090 |
| A106 | A | 86 / 45 | 120.0 nM | 112.3 nM | none | 0.833 |
| D102 | D | 113 / 18 | 124.0 nM | 111.1 nM | none | 0.861 |
| D108 | D | 54 / 77 | 111.5 nM | 116.0 nM | none | 0.861 |

**IDEA-008's recipe was `D102 + A106 + M107 + one of F104/D108`.** Of those, M107 and F104 are measurably
*harmful*, D102/A106/D108 are inert, and the two positions that actually help — **G100 and G101** — are not
in it. Caveat: M107 conserved has n=11, so its effect size is uncertain even though the sign is clear.

## 4. What we were all missing: a better aggregate signal

Counting substitutions at the seven "mutating hurts" positions (W99, D102, G103, F104, Y105, A106, M107)
predicts affinity far better than mutation count alone:

| Predictor | Spearman vs log10 KD, 131 legal variants |
|---|---|
| substitutions at the "hurts" positions | **+0.505** (p < 1e-8) |
| total mutation count | +0.234 (p = 0.007) |

**Revised design guidance:** the identity rule forces ≥4 substitutions, so spend them where they cost least
— **G100 and G101 first** (both associated with *tighter* binding), then D102/A106/D108 (inert) — and
**conserve W99, Y105, F104 and M107**.

That this is partly discoverable from the ledger's own structural work makes it more credible, not less:
EXP-002 found A106/M107/D108 never contact HER2 in 8 crystal structures. Non-contact predicts *tolerance*
(A106 and D108 are indeed inert) but not *benefit* — and M107 shows non-contact positions can still be
costly to change. Contact analysis bounds the damage; it does not pick the substitutions.

---

## 5. What the ledger got right that this folder had missed

- **The ≥4 substitution floor.** EXP-007 derived it independently from 10-mer arithmetic (3/10 = exactly
  70%, and the rule is *strictly under* 70%). My pre-filter carried ≥3 for hours. Correctly caught.
- **EXP-002's structural contact analysis** across 8 crystal structures. Nothing in `NICCOLO/` does this,
  and it is the right way to bound which positions are safe to touch.
- **The framing that the constraint, not the model, costs the 4x.** EXP-007 states this more sharply than my
  own write-up did.
- **It flagged its own weakness** (IDEA-017) rather than leaving it for someone else to find. The error here
  is one of over-reading a small chosen sample, not of hiding it.

## 6. Suggested ledger actions

| Item | Action |
|---|---|
| **TH-019** | revise: hit rate 0.8% not 4.8%; 0 of 131 against the best replicate; state which baseline |
| **TH-009 / TH-015** | the "safe positions" claim needs splitting: tolerated (A106, D108, D102) vs beneficial (G100, G101) vs costly (M107, F104) |
| **IDEA-008** | **do not run as written.** Replace `M107 + F104` with `G100 + G101`; keep A106/D108 as filler |
| **EXP-008** | add a base-rate column to any "mutated in N of N winners" claim |
| **IDEA-017** | **answered** — see §1; no `scp` needed, the numbers are here |
| New standing rule | any "winners share X" claim must report how common X is among non-winners |

All figures reproducible with `python -I eda/p15_ledger_check.py` and `python -I eda/p16_position_test.py`.
