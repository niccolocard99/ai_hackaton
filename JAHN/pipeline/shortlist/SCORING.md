# How the 32 submitted designs were scored, and why that scoring is justified

**Submission:** `submission.csv` — one CSV, both tracks' format, every row carrying its
officially-assigned `track`. **Verdict:** 32/32 eligible on the optimisation track, 0/32 de novo
(`optim_report.csv`, `denovo_report.csv`, the stock team checkers wrapping `abguard`). Composition:
16 from the orchestrated loops, 16 from LEO's design-loop final 20, as instructed.

## The scoring, in one paragraph

Every candidate in the pool — 836 distinct loop designs plus LEO's 20 — was scored by the same
fitted models in one pass: **p_binder** (logistic classifier: will this CDRH3 variant give a
measurable KD at all?) and **pred_log10_kd** (ridge ranker on log10 KD) from `NICCOLO/eda/prefilter.py`,
fitted on the organiser-provided 768-variant SPR table, aggregated by **borda rank over both columns**.
Eligibility was applied before ranking, never inside it: every design keeps ≥3 substitutions in the
CDRH3 window from every provided window (the ">2 mutations" rule, exact vectorised screen). The
final 32 were then chosen half/half — 16 from each source — by walking each ranked list best-first
and keeping a window only if it differs from every already-kept window in ≥3 positions.

## Why each piece is justified

**The models, and their validation.** The classifier and ranker are regularised linear models over a
200-column one-hot encoding of the 10-residue CDRH3 window — the entire library shares one VH/VL
framework, so the encoding spans the design space exactly. Validated with clusters (single-linkage,
Hamming ≤2) held out, so near-duplicates never straddle folds: **AUROC 0.927** (classifier; 0.913
inside the ≥4-substitution legal band), **Spearman +0.687** (ranker; +0.628 in band). Linear beat
gradient boosting and forests on both tasks with 620 labels — the additive form is most of the
signal, measured, not assumed. The leakage audit is clean: no shared sequences across folds,
permuted-label control at 0.486, nested CV within 0.005 AUROC of the fixed choice.

**The aggregation.** Ranking on both columns rather than pred_log10_kd alone was measured on the
179 variants with a known KD, out-of-fold: **+0.729 Spearman vs +0.680**, bootstrap difference +0.049
with 95% CI [+0.013, +0.094] (`NICCOLO/eda/p18-p19`). Ranking on every available score dilutes
(+0.659), so exactly these two columns are used, equal weights, for every candidate.

**The eligibility gates.** The 70% CDRH3-identity rule binds before the mutation count: trastuzumab's
13-mer CDRH3 needs ≥4 substitutions to fall under 70% identity, so the legal band is 4-8 substitutions
— the generator's bounds, the policy's hard floor (`min_subs ≥ 4` is in the allowlist), and the
novelty screen (≥3 from every one of the 525k provided windows) all encode the same rules the
official checker enforces. The official checkers then ran on the final 32 unchanged: **32/32 PASS
optimisation**, and their `track_assigned` is what populates the CSV's `track` column — nobody
assigned a track by hand.

**Diversity, because only four designs get tested.** The selection keeps windows ≥3 substitutions
apart so the four that reach the lab are not four copies of one bet. Ranking after the diversity walk
is by the same borda — the walk only decides *membership*, never *order*.

**Provenance of the pool.** The loop did not start from nothing: NICCOLO's 1M-candidate screen
(batch 1: the top 2,000 ranked windows) and LEO's JapanFold-RFD4 final 233 were the warm-start
elites, and 686 windows the team had already evaluated were pinned as never-regenerate history —
the generation ledger (`LEO/tasks/generation-ledger/`) records every run, its seed, its settings and
its outcome, so no sequence was paid for twice. LEO's 16 came through his own pipeline first
(NanoMelt Tm, epitope risk, fold checks) and were then re-scored here under the identical pass, so
the two halves are directly comparable.

**Structure check, done once, honestly.** The top 5 were re-predicted by AF2-IG via JapanFold from
the 1N8Z initial guess: pLDDT 0.87-0.89 (they fold), ipTM 0.20-0.22, interface PAE 26.4-26.8 Å. As
`NICCOLO/eda/AGENT_GUIDE.md` §6 established on measured data (n=60, then n=180), fold metrics do
**not** separate binders from non-binders for CDRH3-window designs — Boltz-2 ipTM ≥ 0.7 passes 100%
of true non-binders — so the fold is recorded as provenance, never used to rank or gate.

## What this scoring cannot tell you (the honest limits)

- **These are surrogate scores.** The models were fitted on the parent framework's measured variants;
  predictions below log10 KD −1.02 (0.095 nM) are tighter than anything measured in the legal band —
  extrapolation, and the top of any model's list is where model error concentrates.
- **Two of the three judged properties are not scored.** Thermostability and epitope preservation
  have no validated predictor here (NanoMelt was ported but parked: its no-align call fails on a
  feature-count mismatch, and LEO's own Tm numbers travel with his 16 via their descriptions).
  Binding probability and predicted tightness are the only scored properties; the S97 hedge
  (novelty margin at an unmeasured position) and the framework Tm variants were chosen by the team's
  other pipelines, not by this scorer.
- **Only 4 of these 32 will be tested.** The list is ordered so the best-scoring, mutually-diverse
  designs surface first; the loop-side half and LEO's half were kept distinct because they represent
  different bets (surrogate-optimum windows vs. Tm/epitope-hedged portfolio picks).
- Local verdicts are predictions: the organisers' `abguard` at submission time assigns the real track.

## Reproduce the whole ranking

```bash
cd JAHN/pipeline
~/.venvs/abpipe/bin/python scripts/rank_table.py any-designs.csv -o scored.csv   # same fitted models, same borda
```

Every input, seed, setting and outcome is in `LEO/tasks/generation-ledger/`; the runs' full journals
(STEPS/REPORT per round, decisions with reasoning) are under `runs/orchestrated-*` on the VM and
described in `ORCHESTRATOR.md`.
