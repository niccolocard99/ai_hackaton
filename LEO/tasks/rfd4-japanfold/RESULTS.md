# RFD4-Proteina on JapanFold → base validation

**Date:** 2026-11-14 · **Agent:** Claude Opus 5 · **Compute:** JapanFold API (no local GPU) + VM for `abguard`
**Spend:** $0.047 of $1,098 JapanFold credit · 44 designs · ~50 s per design

Ran RFD4-Proteina through the JapanFold API because the team's single H100 was at
100% utilisation with a teammate's Protenix container. JapanFold serves `rfd4`
remotely, so this cost no GPU contention.

Then ran **base validation** (`abguard`) for the first time on any design in this
repo. It overturned both of this task's prior positions.

## Headline

| | |
|---|---|
| RFD4 designs generated | 44 (6 seeds) → **13 distinct** CDRH3s (30% diversity) |
| Framework held fixed | **44/44** — the contig was honoured, verified residue by residue |
| Eligible, optimisation track | **9/18** candidates |
| Eligible, de novo track | **0/18** — every design is 96–98% identical to trastuzumab |
| Geometry-driven designs eligible | **0/4** |

## 1. `abguard`'s CDRH3 is a 13-mer, and neither earlier guess was right

This task spent considerable effort on whether the CDRH3 window was a 10-mer or an
11-mer. It is **neither**:

```
SRWGGDGFYAMDY     13 residues (IMGT), reported in the report's `cdrh3` column
```

Consequences:

- **Identity denominator is 13.** 4 substitutions → 9/13 = 69.2% (passes `<70%`);
  3 → 10/13 = 76.9% (fails). So **≥4 substitutions** remains the floor, but for a
  different reason than either the 10-mer or 11-mer arithmetic gave.
- The pre-filter's window (`TZ = "WGGDGFYAMD"`, L=10) is **not** `abguard`'s window.
  Positions 97, 98 (`S`, `R`) and 109 (`Y`) are inside `abguard`'s CDRH3 and were
  never considered designable here. That is 3 unexploited positions.
- `EXP-007`'s correction (3 subs = exactly 70% on a 10-mer = fail) reaches the right
  answer — ≥4 — from the wrong denominator.

## 2. The binding constraint is novelty distance, not identity

Every geometry-driven design **passed** the identity rule at 69.2% and was still
rejected:

```
FAIL LEO-01  CDRH3 only 1 mutation(s) from tz_her2_facs_3classes_524346
             (SRWGGDGFYSLEY); needs >= 3
```

`cdrh3_min_distance` — Hamming distance to the nearest of the 525,063 known CDRH3s
— is the rule that actually binds. All four geometry designs sit at **distance 1**.

**This was predicted and then confirmed.** [`findings/07-rule-verification.md`](../litreview-optimisation/findings/07-rule-verification.md)
§4 argued that confining mutations to the few contact-free positions shrinks the
design space into exactly the region a 524,346-variant binder-selected library
covers densely. That is what happened. The geometric reasoning was sound and the
designs are still ineligible, because a smaller space is a more crowded one.

## 3. RFD4 chose the opposite positions, and it was right about eligibility

The structural analysis in the sibling task concluded `W99` and `Y105` were
untouchable (alanine-scan hotspots, 7/8 contact frequency) and that `A106`,
`M107`, `D108`, `Y109` were the only free positions. RFD4, given only the
structure and a contig, did close to the reverse:

| design | substitutions | min_dist | identity | p_binder |
|---|---|---|---|---|
| **RFD4-02** | W99Y G101L D102Y Y105W A106E | 3 | 61.5% | **0.843** |
| **RFD4-01** | W99Y G101P D102Y Y105W A106E | 3 | 61.5% | **0.805** |
| RFD4-11 | W99Y G101L D102H G103S F104A Y105W A106Y M107F | 3 | 38.5% | 0.679 |
| RFD4-05 | W99Y G101L D102G G103S Y105W A106E | 3 | 53.8% | 0.615 |
| RFD4-08 | G100D G101P D102F F104T Y105W | 3 | 61.5% | 0.546 |
| RFD4-03 | W99Y G100R G101S D102L G103Y F104D A106E | 3 | 46.2% | 0.321 |
| RFD4-10 | W99S G100T G101S D102S G103V F104S Y105W M107F | **4** | 38.5% | 0.295 |
| RFD4-13 ⚠ | 9 subs, `N-F-T` sequon | **4** | 30.8% | 0.154 |
| RFD4-12 ⚠ | 9 subs, `N-A-S` sequon | 3 | 38.5% | 0.082 |

**8 of 13 distinct designs put Tyr at position 99**, across independent seeds. RFD4
converges on mutating the hotspot that the contact analysis said to protect. Two
readings, not yet separable: either the geometric hotspot argument is wrong, or
RFD4 is trading measured affinity for the novelty distance the rules demand.
Nothing here decides it, and it is the most interesting open question in this task.

⚠ **RFD4-12 and RFD4-13 carry N-glycosylation sequons** (`N-A-S`, `N-F-T`). They
pass `abguard` and should still be dropped.

## 4. Rank on `p_binder`, never on predicted KD

The pre-filter's ranker **cannot resolve the regime we are scored in**. Leak-free
grouped CV ([`subnm_resolution.py`](subnm_resolution.py)):

| regime | n | Spearman | median err |
|---|---|---|---|
| full range | 179 | +0.687 | 0.46 log10 |
| ≥4 subs (legal band) | 131 | +0.628 | 0.39 |
| measured < 10 nM | 45 | +0.415 | 0.58 |
| measured < 1 nM | 17 | **+0.167** | 1.25 |
| < 1 nM **and** ≥4 subs | 9 | **−0.283** | 1.25 |

**The ranker's predicted floor is 1.135 nM. Trastuzumab measures 0.152 nM.**
Seventeen variants measure sub-nanomolar; it predicts none below 1.1 nM. Its
median error there is 17.7×, and its rank correlation inside the legal band is
*negative*.

Confirmed on controls it was trained on — a resubstitution check it fails:

| | measured | predicted |
|---|---|---|
| `WGGSGPYLFD` (best legal) | **0.095 nM** | 2.418 nM — ranked *worst* |
| trastuzumab | 0.152 nM | 1.800 nM |
| `WGAEGFYYWD` | 0.368 nM | 1.864 nM |

The **classifier** is a different story: AUROC 0.913 inside the ≥4 band, a real
validated signal for "binds measurably at all". So `p_binder` is usable and
`pred_kd_nM` is not. This is the second independent resolution wall in this
project, after trastuzumab's k_off sitting within ~5× of the Biacore T200 floor.
Both say the same thing: **do not stake the pitch on winning binding.**

## 5. A tension to resolve before ranking

Novelty margin and predicted bindability are **anti-correlated**. Seven of nine
eligible designs sit at `min_distance` exactly 3 — zero margin against a gate
whose official run happens at submission. The two with margin (distance 4) are
the two worst binders. There is no design here that is both safe and good.

## Reproduce

```bash
# 1. design on JapanFold (no local GPU; ~$0.008 per 8 designs)
CONTIG="A1-107,/0,B1-98,10-10,B109-120,/0,C489-580,/0,D591-607"
japanfold design --model rfd4 --protocol rfd4-scaffold \
  --structure apheris/rfd4/Data/12-her2/her2_fv_1n8z.cif \
  --contig "$CONTIG" --num-designs 8 --seed 1 \
  --param recipe=high/front_loaded --wait --out runs/

# 2. harvest, verify the framework, gate on novelty BEFORE scoring
python -I harvest_rfd4.py runs/

# 3. build the base-validation table and run both tracks
python -I build_designs_csv.py
python scripts/check_optimisation.py designs.csv -o reports/optimisation.csv \
    --provided data/trastuzumab_her2_spr_all_768.csv
python scripts/check_denovo.py designs.csv -o reports/denovo.csv
```

Needs `data/` symlinked to the SPR-768 table and the FACS splits, and `prefilter.py`
from `origin/nick:NICCOLO/eda/`. `abguard` is at `~/.venvs/abguard/bin/` on the VM.

## Files

| Path | What |
|---|---|
| `designs.csv` | 18 candidates in base-validation format |
| `reports/optimisation.csv` | **the authoritative verdict**, 9 PASS |
| `reports/denovo.csv` | 0 PASS, as expected — never reclassified |
| `rfd4_designs_gated.csv` | 13 distinct RFD4 CDRH3s with gates and scores |
| `leo_candidates_scored.csv` | the 4 geometry designs + measured controls |
| `japanfold/run-s1`, `run-more` | 44 design CIFs, job manifests, metrics |
| `harvest_rfd4.py` | CIF → CDRH3, framework check, novelty gate, scoring |
| `build_designs_csv.py` | parent VH/VL from the SPR table, CDRH3 swap |
| `subnm_resolution.py` | the resolution analysis in §4 |

## Open

1. **Positions 97, 98, 109 are inside `abguard`'s CDRH3 and unexplored.** Designing
   over 13 positions instead of 10 is the cheapest way to buy novelty distance
   without spending a contact.
2. **Resolve the W99Y disagreement.** RFD4 and the contact analysis disagree about
   the most important residue in the loop.
3. **Nothing is folded.** RFD4 orders designs structurally, not by confidence. But
   `EXP-003/004/006` show folding metrics barely rank HER2 binders (Boltz-2 AUROC
   0.645; decoys scored *higher*), so a refold is unlikely to settle anything.
4. **No thermostability prediction** for any of these, and the CDRH3 rewrite is
   where the bH1 series lost 6–9 °C.
