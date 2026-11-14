# Generate -> filter -> score loop with Niccolo's pre-filter

```bash
cd JAHN/pipeline
~/.venvs/abpipe/bin/abpipe iterate configs/niccolo_loop.yaml --rounds 3 --select 60 --out runs/loop --propose-next
```

`configs/niccolo_loop.yaml` is the whole loop; swapping a module is editing a `name:` line.

| Stage | Module | What it does |
|---|---|---|
| generate | `scored_evolution` | Proposes CDRH3 windows. Prior from the classifier and ranker weights plus FACS 'high' enrichment; from round 1 it also learns from the previous round's selected windows. |
| filter | `framework_integrity`, `window_novelty`, `prefilter_gate`, `dedupe` | Cheap first, and all of it O(1) per design. `window_novelty` is an exact vectorised pre-screen of the ">2 mutations from every provided window" rule; `prefilter_gate` is step 2 of Niccolo's cascade. |
| score | `niccolo_prefilter` | `score__p_binder`, `score__pred_log10_kd` from `NICCOLO/eda/prefilter.py`, unmodified. |
| rank | `borda` on `score__pred_log10_kd` | Gate on p_binder, rank on predicted KD (`NICCOLO/eda/AGENT_GUIDE.md`, section 2). |
| validate | `abguard_hard`, `new_liabilities`, `track_check` | BASE VALIDATION, **after ranking**, on the top `validate_top` (2000) only. abguard is authoritative but expensive — 37 designs/s with the novelty cross join, ~7.5 h per 1M — and only the top of the list can be submitted. Selection is restricted to its survivors, so this is a reordering, not a weakening (`tests/test_validate_stage.py`). See `FAST_NOVELTY.md`. |

## How a round feeds the next

1. Round r's new survivors are ranked **together with the previous round's selected windows** on their absolute
   scores (elitism). Only the single best candidate is guaranteed not to get worse: the selected set as a whole can
   shift, because a child within 2 substitutions of an elite displaces it when it scores better.
2. Selection walks the best `--pool-frac` (default 10%) of that ranking best-first and keeps a window only if it
   differs from every one already kept in >= 3 positions and no single substitution appears in more than 25% of the
   list that is **returned** (the cap is re-derived until the list and its cap agree). It returns **fewer** than
   `--select` rather than pad the list with mid-pool candidates or let one substitution dominate; the log says so.
3. Those windows become `parents` for the next generator: residue, mutated-position and substitution-count
   distributions are pulled towards them (`eta`), and children come from `sample` / `mutate` / `cross`.
4. Every window ever generated is `history` and is never proposed again.

Re-running with a larger `--rounds` resumes after the last finished round. `--propose-next` additionally writes
`round_NN/proposals.csv`: the next round's candidates after generation and all filters, not yet scored.

Outputs per round: `selected.csv` (scored, ranked, diversified), `ranked.csv`, `validated.csv` (the top of the
ranked list that cleared BASE VALIDATION), `02_rejections.parquet` and `05_validate_rejections.parquet` (why each
candidate was dropped), `manifest.json`. `loop_summary.json` has the per-round funnel and score statistics.

## Setup this needs

- The `abpipe` venv with `pandas` and `scikit-learn` (`uv pip install pandas scikit-learn`); `prefilter.py` imports them.
- FACS clonotype splits (`train.csv`, `val.csv`): `splits_dir` in the config, `PREFILTER_SPLITS_DIR`, or
  `~/nick-work/splits` on the team VM.
- `NICCOLO/eda/prefilter.py` on the branch you run from (`git pull`). Module versions include its hash, so a changed
  model never serves stale cached scores.

## What the loop does not know

Read `NICCOLO/eda/AGENT_GUIDE.md` sections 3-5 before trusting a ranking.

- It improves the **surrogate's** score. The classifier (AUROC 0.927, 0.913 in the >=4-substitution band) and ranker
  (Spearman +0.687, +0.628 in band) were validated on library variants with cluster-grouped CV. The loop deliberately
  pushes into the region the model has not seen, where that validation does not strictly apply.
- Both models are linear in a one-hot encoding. Their optimum under the constraints is nearly reachable in one shot,
  so extra rounds mostly concentrate and diversify around it. The loop earns its keep when you add scorers that are
  not additive (thermostability, expression, anything with feedback from the lab) to the config.
- `pred_log10_kd` is a ranking, not a KD; it compresses at the tight end. **Predictions beyond the measured range are
  extrapolation.** In the legal band (>= 4 substitutions) 658 variants were measured, 131 gave a KD, the tightest is
  0.095 nM (log10 -1.02) and only 9 are below 1 nM. A predicted log10 KD below -1.02 is tighter than anything measured
  there. Because the loop selects for the model's best scores, the top of its list is the part most likely to be
  model error rather than binding.
- Thermostability and epitope preservation, two of the three judged properties, are **not scored**. Nor is expression.
- A window that passes `window_novelty` and `track_check` here is a prediction of eligibility; abguard's own
  assignment at submission time is the real one.
