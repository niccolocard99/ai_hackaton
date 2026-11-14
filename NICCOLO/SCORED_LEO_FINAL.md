# Scoring LEO's 20 final designs with frozen stage 1

Date 2026-11-14 · Source `~/wt/leo-design-loop/LEO/tasks/design-loop/out/FINAL.csv` (md5 matches the copy
under `~/hackathon-team-3/`, so the two are the same file) · Script `eda/p21_score_leo_final.py` · Output
`eda/figs/leo_FINAL_scored.csv`

## First, they are not de novo designs

All 20 pass the abguard **optimisation** checks and fail the de novo script, and `FINAL.md` says so itself.
The `track_assigned = novel` column is abguard's *novelty_class* (set from `novelty_class` in
`NICCOLO/scripts/check_optimisation.py`), not a track assignment — so nothing is mis-tracked. They are
optimisation-track designs on trastuzumab's framework, which is exactly stage 1's scope.

## They do carry a different parent, and the pre-filter cannot see it

LEO's loop adds "addons" on top of the CDRH3 substitutions:

- **`D28K`** — a light-chain substitution. 8 of 20 designs have a VL that differs from trastuzumab's.
- **`fwTm`** — heavy-framework substitutions outside the CDRH3 window. 12 substitutions across 6 designs.

Stage 1's pre-filter reads **the 10-residue CDRH3 window only**. It was trained on a library where all 732
variants share exactly one VH framework and one VL, so it is blind to both addons by construction, and any
framework change puts the design outside its training distribution.

This is visible directly in the list: **6 of the 14 distinct windows appear twice** under different addons,
and the pre-filter gives each pair identical scores — correctly, since it sees nothing else.

## Scores, and where they agree with LEO's

| | |
|---|---|
| Spearman(LEO `kd_mean`, pipeline `pred_log10_kd`) | **+0.941**, p < 1e-9 |
| Spearman(LEO final rank, pipeline rank) | **+0.030**, p = 0.90 |

The two binding models agree almost perfectly. That is real mutual corroboration: two independently built
surrogates, one of them validated out-of-fold against 179 measured KDs, put these windows in the same order.

The *final ranks* do not agree at all, and they should not. LEO ranks on three properties; stage 1 sees one.
**Stage 1 must not be used to re-rank this list** — it has nothing to say about epitope risk or Tm, which is
where LEO's ranking differs from pure binding.

## Checks that stage 1 can contribute

**Interpolation: passes.** Every design is 3–4 substitutions from its nearest variant in the ranker's
179-variant training set; **none is at ≥5**, where the measured signal is Spearman −0.057. So these scores
fall inside the range where the ranker has demonstrated power — unusual, and worth knowing.

**Extrapolation: flagged on a third of the list.** 7 of 20 have `pred_log10_kd` below −1.02, i.e. predicted
tighter than *any* measured variant in the legal band (0.095 nM). LEO's own `kd_mean` flags 6. Both models
are extrapolating there; the ordering may still be useful, the absolute values are not.

## One thing worth deciding before submission

**`D28K` is scored as a no-op by all three of LEO's models.** For the three windows that appear with and
without it, `kd_mean`, `dtm_vs_wt` and `epi_risk` are *identical* — yet the pairs land at ranks 2/10, 5/17
and 19/4:

| window | addons | kd_mean | dtm_vs_wt | epi_risk | ranks |
|---|---|---|---|---|---|
| WKSQGFYEMI | `D28K` vs `cdrh3only` | −1.153 both | 0.05 both | 1.375 both | 19 and 4 |
| WTVFGFYEME | `D28K` vs `cdrh3only` | −0.493 both | 1.55 both | 1.375 both | 5 and 17 |
| WWVEGFYEMI | `D28K` vs `cdrh3only` | −1.091 both | −0.32 both | 1.375 both | 2 and 10 |

So three of twenty submission slots are spent on a light-chain change that nothing in the scoring can
distinguish, and the pairs are separated only by tie-breaking. Either the models should be made to see
`D28K` or it should not consume a slot.

`fwTm` is different and looks like a genuine play: those pairs differ by **+2.4 to +2.6 °C** in `dtm_vs_wt`,
which matches `NICCOLO/thermo/`'s finding that framework substitutions are the only thermostability lever
available — CDRH3 cannot move Tm at all.

## Net

Stage 1 corroborates LEO's binding ranking strongly (+0.941) and confirms the designs sit inside the
ranker's validated distance range. It cannot and should not re-rank the list. The two things to act on are
the `D28K` no-op slots and the seven extrapolated binding predictions.
