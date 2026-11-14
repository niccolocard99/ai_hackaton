# The final dataset: 32 best sequences under one scoring

**File:** `final_32.csv` · **Compiled:** 14 Nov 2026, after the hard stop · **Scoring:** identical for
every candidate — NICCOLO's frozen stage-1 pre-filter (`score__p_binder` from the classifier,
`score__pred_log10_kd` from the ranker, both fitted on the 768-variant SPR table) aggregated with the
frozen `borda(pred_log10_kd, p_binder)` — the combination measured at +0.729 Spearman out-of-fold
(PR #13, `NICCOLO/eda/p18-p19`).

## How the 32 were chosen

1. **Pool** every selection the campaign produced: the batch-1 run (2 rounds), the fast1 run
   (10 rounds), the four sat2 replicas (round 0), LEO's design-loop final 20, and LEO's
   JapanFold-RFD4 final 233 — 1,578 rows, **1,089 distinct** (dedup on heavy+light chain).
2. **Eligibility**: every candidate must be ≥3 substitutions in the CDRH3 window from every provided
   window (the >2-mutations rule; exact vectorised screen). 1,059 legal.
3. **One scoring pass** over all of them (same fitted models, same two-column borda).
4. **Diversity walk**: best-first, keeping a window only if it is ≥3 substitutions from every window
   already kept (the campaign's own rule) until 32.

## What is in it

| origin | count |
|---|---|
| fast1 (the self-optimising loop, 10 rounds, validate-after-select) | 25 |
| batch1 (the loop's first run, official checkers on every survivor) | 5 |
| sat2 (the four saturation replicas, Leo's JF-RFD4 elites leading) | 2 |

Top score line: predicted log10 KD −1.07, p_binder 0.97 (tightest predicted binder in the set;
−1.07 ≈ 0.085 nM is beyond the tightest measured variant in the legal band — extrapolation, see
`NICCOLO/eda/AGENT_GUIDE.md`).

## Honest caveats, carried from the measured evidence

- **All scores are the surrogate's**, trained on the parent framework's CDRH3 variants. The pre-filter
  was never validated on novel frameworks; LEO's RFD4 designs (novelty-rich, classifier-poor — the
  handoff's central tension) rank low under this scoring and took no slots, even though they carry
  the best novelty margins. They remain available, separately scored, in
  `leo_final20_scored.csv` and `bootstrap/leo_jf_final.csv` — the submission's final mix is a team
  decision, not a scoring one.
- **Only the sequence-level rules were applied here.** The 32 still need the official
  `check_denovo.py` / `check_optimisation.py` (abguard assigns the real track at submission time);
  the sat2 rows were validated by the loop's own `track_check` before selection.
- Everything in this dataset descends from bootstrapped team work (LEO's evaluated 686 windows as
  never-re-propose history; NICCOLO's 1M-screen batch and LEO's JF-RFD4 final 233 as elites) and is
  logged in `LEO/tasks/generation-ledger/` (batch1, fast1, fast1-cont, sat2 rows; the failed runs
  carry their reasons).

## Reproduce

```bash
cd JAHN/pipeline && ~/.venvs/abpipe/bin/python scripts/rank_table.py <any-designs.csv> -o scored.csv
# the compilation itself pooled runs/orchestrated-{batch1,fast1}/round_*/selected.csv,
# runs/sat2-s*/selected_all.csv, LEO's two lists, then the frozen scorer + diversity walk above.
```
