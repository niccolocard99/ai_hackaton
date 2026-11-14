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

## 4. The position table is a readout of the ranker, not new predictive signal

An earlier version of this file quoted "substitutions at the costly positions → Spearman **+0.505**, versus
+0.234 for mutation count" and called it a better aggregate signal. **That framing was wrong.** The +0.505
is an in-sample univariate correlation on the 131 legal variants, with no cross-validation, so it is not
comparable to the ranker's out-of-fold **+0.683**. It is lower, not higher.

Two tests (`eda/p17_is_it_new.py`):

| Test | Result |
|---|---|
| Ranker's implied per-position cost vs the measured effect | **Spearman +0.794**, p=0.006, agreeing on **8 of 10** positions |
| Ranker + costly-position counter, out-of-fold | +0.683 → +0.684, 95% CI [−0.005, +0.011] — **no difference** |

The ridge already encodes the per-position structure in its coefficients. The counter is redundant.

**What the position table is actually for**, then:

1. **Design guidance, which the ranker cannot give.** The ranker *scores* finished sequences; it does not
   tell you where to spend the ≥4 substitutions the identity rule forces when *generating* a design. The
   position table does, and that is what IDEA-008 needs.
2. **Independent corroboration that the ranker is not fitting noise.** Ridge coefficients from 179 variants
   and a rank test on 131 legal variants are different estimators, and they agree on 8 of 10 positions.
3. **The two disagreements are the interesting part — and they are exactly EXP-008's two positions.** At
   D102 and A106 the ranker says mutating helps while the measured test says inert. Both agree they are not
   costly, and neither supports building a recipe around them.

None of this weakens §1–§3, which are facts about the landscape and a methodological error, not Spearman
claims: the hit rate is 1 in 131, the "8 of 8" convergence is a base-rate artefact, and M107/F104 are
measurably costly to mutate.

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
| New standing rule | compare a new predictor against the existing model out-of-fold, never an in-sample univariate r against a CV score |

All figures reproducible with `python -I eda/p15_ledger_check.py` and `python -I eda/p16_position_test.py`.

---

# Appendix: JAHN's scoring pipeline — integration and effectiveness

Checked `JAHN/pipeline` (abpipe) on 2026-11-14. Scripts: `eda/p18_jahn_eval.py`, `eda/p19_agg_bootstrap.py`.

## Integration: careful, and I broke it

The wrapper (`abpipe/_niccolo.py`) is good work. It imports `prefilter.py` rather than reimplementing it,
handles its cwd dependency, fingerprints the source file into the cache key so a changed model cannot serve
stale scores, enforces the trastuzumab-window scope via `check_parent`, rejects unscorable rows instead of
passing them, and its tests pin my published numbers as a contract.

**I broke all of it.** Defaulting `PreFilter` to `use_facs=False` removed the `p_facs_high` column that
`AGENT_GUIDE.md` had documented as an output. `score_windows` read that key unconditionally, so the scorer,
the gate filter and `emit_facs` every one failed with `KeyError: 'p_facs_high'`. The 13 integration tests
did not catch it because they are gated on `abguard`, which is not installed locally, so they silently skip.

Fixed in this change:

- `score_windows` now follows the keys the fitted model actually emits, so both model variants work.
- `has_facs(fit)` added; the scorer raises a clear `PrefilterUnavailable` when `emit_facs=True` has no FACS
  model, instead of a `KeyError`.
- `find_splits(required=False)` and `get_prefilter(use_facs=...)`: the FACS clonotype splits are no longer a
  hard requirement, since the default model does not use them. The pipeline now runs off the SPR table alone
  — one 200 KB file instead of an 11 MB extract from a 2.5 GB zip.

Verified end to end: scorer, gate and the no-splits path all produce sensible numbers, and his 23 unit tests
still pass.

## Effectiveness, measured against real KD

179 variants with a measured KD. **The pre-filter columns are supplied out-of-fold** (clusters held out) so
the comparison is not rigged in its favour; his scorers never see KD at all.

| Signal | Spearman vs log10 KD | Top-20 median KD |
|---|---|---|
| `score__p_binder` (pre-filter, OOF) | **0.720** | 2.2 nM |
| `score__pred_log10_kd` (pre-filter, OOF) | 0.680 | 1.9 nM |
| `score__binder_distance` (his, KD-blind) | **0.595** | 8.6 nM |
| `score__parent_distance` (≈ mutation count) | 0.412 | 6.5 nM |
| `score__neighbor_binder_frac` (his) | **0.100 (p=0.18)** | 27.7 nM |

Set median is 86 nM. Two things stand out.

**`binder_distance` at 0.595 is a genuinely good KD-blind signal** — nearest-neighbour distance to a
FACS-measured high binder, using no affinity data at all. It does not add to the ranker, but it is the best
thing available for a design the ranker cannot speak to.

**`neighbor_binder_frac` does not work.** Spearman 0.100, p=0.18, and it actively dilutes: ranker +
that scores 0.508 against 0.680 for the ranker alone. It should be dropped or reworked, and it must not sit
in a default `borda` column set.

### Aggregation: combining the two pre-filter scores helps

| Ranking | Spearman, all 179 | Spearman, legal band (n=131) |
|---|---|---|
| `pred_log10_kd` alone (what `niccolo_loop.yaml` shipped) | +0.680 | +0.618 |
| **borda(`pred_log10_kd`, `p_binder`)** | **+0.729** | **+0.659** |
| his 3 scorers, no pre-filter | +0.431 | +0.235 |
| all 5 columns | +0.659 | +0.576 |

Cluster-bootstrapped difference for adding `p_binder`: **+0.049, 95% CI [+0.013, +0.094]** on all 179 —
real. Inside the legal band: +0.041, CI [-0.005, +0.092] — **not established**, though the sign is the same
and it never hurt in any slice. The top-20 median KD difference (0.58 nM vs 1.92 nM) is **not** significant,
CI [-0.91, +0.60] in log10 units; top-20 medians at this n are too noisy to rank on.

`niccolo_loop.yaml` now ranks on both columns, with that caveat recorded inline.

**Note for my own pipeline:** I used `p_binder` only as a gate and ranked on `pred_log10_kd`. Ranking on
both is better, and I had not tested it. That came out of his aggregator, not mine.

**Dilution is real:** all five columns (+0.659) is worse than the best two (+0.729). `borda` with no
`columns` set — the `baseline.yaml` default — ranks on every registered score equally, which drags the two
good signals down with the weak ones. Column sets should be chosen, not defaulted.

## Re-check after JAHN's selection fixes (commit `8cd70a4`)

He fixed two defects in `select_diverse` found in his first full 3-round run. Re-verified on my side:

**Nothing of mine broke.** His 25 tests pass (2 new regression tests), and the pre-filter scorer, the gate
filter and the no-splits path all still work. He did not touch my `_niccolo.py` fixes, and the
`niccolo_loop.yaml` column change survived.

**His factual claims about my data are exact.** 658 legal 10-mer variants, 131 with a KD, 9 of those below
1 nM, tightest 0.095 nM = log10 −1.02. His extrapolation warning in `LOOP.md` — that a predicted log10 KD
below −1.02 is tighter than anything measured in the legal band — is correctly derived and is a caveat I
had not written down myself.

**The fix is a measured improvement.** Ranking the 131 legal variants with measured KD by the pipeline's own
borda(`pred_log10_kd`, `p_binder`) using out-of-fold scores, then running both versions of `select_diverse`:

| n asked | old: median KD | new: median KD | plain top-N | old median pool rank | new |
|---|---|---|---|---|---|
| 10 | 10.5 nM | **3.9 nM** | 3.9 nM | 7 | 5 |
| 20 | 15.2 nM | **4.6 nM** | 6.5 nM | 13 | 10 |
| 30 | 16.5 nM | **12.0 nM** | 12.9 nM | 20 | 14 |

Cluster bootstrap on the difference at n=20: **−0.504 log10, 95% CI [−1.212, −0.041]** — entirely below
zero, so the new selector really does pick tighter binders, by a median factor of about **3.2x**. His
"starved walk" diagnosis is confirmed: the old version drew from deeper in the pool at every size.

At n=20 the new selector also beats a plain top-20 (4.6 vs 6.5 nM) while still enforcing the diversity cap,
so the constraint is not costing affinity here — it is dropping redundant near-duplicates.

**The trade-off to plan around:** the new version returns a shorter list rather than padding. Across
bootstrap resamples it returned a median of **12 of 20** requested. If we need 20 submissions, either the
candidate pool has to be larger or `--pool-frac` raised; asking for 20 will not reliably give 20.
