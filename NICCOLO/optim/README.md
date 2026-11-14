# Optimisation-track CDRH3 variants: generate -> score -> RF3 fold -> rank

Scope: the 10-residue CDRH3 window (`WGGDGFYAMD`, VH 99-108) on the trastuzumab VH/VL framework, i.e. the optimisation track.
The prefilter classifier (`NICCOLO/eda/prefilter.py`) does not apply to de novo designs (its own README says so; `score()` also raises on,
or silently mis-windows, anything that is not exactly 10 residues).

Working copy: `/home/hackathon/nick-work/optim` on the VM (scripts use absolute paths into it; it also holds the 173 MB RF3 archive and the
unzipped folds). This directory in the repo, `NICCOLO/optim/` on branch `niccolo/optim-variants`, is a snapshot of the code and the small outputs.

## Where it stands (2026-11-14)
| Step | State |
|---|---|
| 120 variants generated and scored; official `check_optimisation.py` in fast mode and `--provided` mode | done, 120/120 pass (the reordered file passes again, fast mode, exit 0) |
| RF3 on JapanFold: 127 folds = 120 designs + 7 SPR controls | done. Job `d9c4b0fcd6d2b102a7950e3a9219d51c`, 127/127 ok. 7 min 4 s from submit to finished (8 min 15 s with the 173 MB download), $0.347 |
| Folds scored against the 1N8Z crystal (epitope, pose, interface PAE) | done: `out/rf3_structure_scores.csv` |
| Re-scored with the current origin/main prefilter (`use_facs=False`) | done: p_binder identical (max difference 3e-8), pred_log10_kd Spearman 0.989 |
| Final order | `out/final_variants_ranked.csv`; `out/final_variants_ranked_abguard.csv` is the same order in submission format (`head -n N+1` = top N) |
| Boltz-2 (`jobs/boltz2_payload.json`) | **not run and not recommended**, see below |

## Result
- **9 of the 120 designs have RF3 ipTM >= 0.70** (0.720-0.795; wild type 0.594). All nine also have the crystal-like pose (Fv CA RMSD 1.6-4.8 A to 1N8Z after
  fitting HER2, 84-100% of the 19 crystal epitope residues contacted) and a tight interface PAE (7.0-10.2 A; the wild-type fold has 19.4 A).
  They are final ranks 1-9: `WWEEGFYAMM, WYNEGFYTMI, WTYDGFYEMG, WPYEGFYAME, WGERGFYEMM, WEFSGFYEMD, WIVTGFYHMD, WGEWWFYEMI, WAWEGFYHMD`.
- Why that tail matters (Jahn's 273 SPR-labelled RF3 folds, `JAHN/pipeline/data/calibration_rf3.csv`, same model and HER2 construct, WT ipTM 0.594 in both runs):
  RF3 ipTM has AUROC 0.82 for >= 4 substitutions, and all of it is in the tail. The highest non-binder of 136 is 0.622; all 7 variants above 0.70 are binders.
  Below 0.62 it says nothing: 103 of 235 variants with ipTM <= 0.6 bind (44%, against a 50% base rate).
- RF3 reproduces the crystal pose for wild type (Fv CA RMSD 2.2 A, 18 of 19 epitope residues contacted). 44 of 120 designs stay within 5 A of it; the rest dock
  elsewhere, but those are low-confidence folds, and pose recovery has AUROC 0.63 for binder vs non-binder in Jahn's set (a measured binder control, `WGGDGGNTYD`,
  docks 36 A away). Do not use pose as a veto.
- Classifier and RF3 are nearly independent (Spearman |rho| <= 0.22 between p_binder / pred_log10_kd / n_mut and ipTM, pair ipTM, interface PAE), so they are two separate lines of evidence.
  The nine RF3-confident designs sit at classifier ranks 7-106 of 120; the classifier is saturated inside this pool (p_binder 0.89-0.99), so it cannot order them well.

## Ranking rule (`rank_final.py`)
Tier 1 = RF3 ipTM >= 0.70. Tier 2 = the rest. Inside a tier, order by the team's frozen stage-1 criterion (`NICCOLO/PIPELINE.md`): borda over `pred_log10_kd` and
`p_binder` of the current prefilter. The classifier orders, RF3 promotes. The first selection (`gen_score.py`) gated at `p_binder >= 0.80` and ranked on `pred_log10_kd`
alone, which is why the borda order differs (13 of its first 20 are also in the borda top 20).

## Files
| Path | What |
|---|---|
| `gen_score.py` | sample 4M variants (k=4..7 substitutions) -> hard rules -> `PreFilter.score` (da498b8 version) -> gate `p_binder>=0.80` -> hole check -> diversify |
| `build_submission.py` | windows -> full paired VH/VL -> official `check_optimisation.py` (fast and `--provided`) -> JapanFold payloads |
| `jobs/submit_japanfold.py` | submit, wait, download; records the job id and timings. **Dry-run unless `--go`** |
| `score_structures.py` | folds vs the 1N8Z crystal: epitope recall/precision, CDRH3 contacts, Fv pose RMSD, interface PAE |
| `rescore_current_prefilter.py` | re-scores the 120 with origin/main's prefilter; borda rank; distance to the SPR training windows |
| `rank_final.py` | tiering and final order |
| `sync_repo.sh` | pull origin/main, copy the whitelisted files here, commit, push `niccolo/optim-variants` |
| `out/final_variants.csv` | the 120 designs in abguard input format, original order |
| `out/final_variants_ranked.csv`, `..._ranked_abguard.csv` | final order with all scores; same order in abguard format |
| `out/final_variants_scored.csv` | classifier scores at selection time (da498b8 prefilter) and the official check columns |
| `out/final_variants_rescored_current_prefilter.csv` | current-prefilter scores, borda rank, `min_hamming_to_spr_train` |
| `out/rf3_structure_scores.csv` | per fold: ipTM, pair ipTM (A-H, A-L), epitope metrics, pose RMSD, PAE; 127 rows, controls included |
| `out/pool_stats.json`, `out/pool_top5000_by_pred_kd_before_hole_check.csv` | counts at every generation stage; the ranked pool before the hole check |
| `jobs/rf3_payload.json`, `jobs/boltz2_payload.json`, `jobs/targets_manifest.csv` | `POST /v1/predictions` bodies (127 targets, chains `H`,`L`,`A`, HER2 domain IV 119 aa) and the target list |
| `jobs/results/rf3_payload/{job.json,results.json}` | the job record and per-fold scores. `archive.zip` (173 MB) and the folds stay on the VM |

## Rules the variants satisfy (checked by the official script, not by my code)
`NICCOLO/scripts/check_optimisation.py out/final_variants.csv --strict-exit`: 120/120 eligible, exit 0; `--provided` (whole VH+VL against every antibody in the SPR table) also passes
for all 120. IMGT CDRH3 identity to trastuzumab must be < 70% on the 13-mer, which needs >= 4 substitutions (3 gives 76.9%). Every design is >= 3 edits from all 525k provided
CDRH3s and from the abguard reference. No Cys, no new N-X-S/T sequon, no substitution in more than 25% of the list at any prefix. The checkers did not change between da498b8 and origin/main 5136e34.

## Read these numbers correctly
- `pred_log10_kd` is a ranking, not a KD. The ranker compresses at the tight end (trastuzumab at 1.8 nM vs 0.03-0.18 nM measured) and the team documents that it extrapolates
  below log10 KD -1.02 (2 of 120 designs: `WGSQGFYEMN`, `WSSEGFYEMF`) and is blind beyond ~4 substitutions from its training set (here 3-5: 56 designs at 3, 51 at 4, 13 at 5).
- `p_binder` = P(measurable KD) from a one-hot logistic model fitted on 620 SPR outcomes. The pool's top is the model's own favourites (an additive model rewards stacking A106E,
  D102Q/E, G101S), so expect optimism; it is a prefilter, not evidence of binding.
- RF3 ipTM is high-precision and low-recall. >= 0.70 is strong evidence; below ~0.62 is uninformative, so tier 2 means "not screened by RF3", not "predicted non-binder". The six SPR
  controls here (three binders, three non-binders) all sit at 0.581-0.592, which is exactly the uninformative zone. 7 of 7 is a small sample (95% lower bound about 0.6), the
  calibration set is balanced 136/136, and these designs are novel sequences, so do not read tier 1 as a probability.
- Boltz-2 has no screening power in Jahn's calibration (ipTM AUROC 0.53; non-binders reach 0.96 and 135 of 136 are above 0.62). Do not spend the run.
- JapanFold indexes `pair_chains_iptm` by SORTED chain id (A, H, L), not input order. PAE tokens follow the file order (H, L, A).
- Thermostability was deliberately ignored (CDRH3 does not move Tm per `NICCOLO/thermo/README.md`).

## Reproduce
```
cd /home/hackathon/nick-work/optim
/home/hackathon/.venvs/abmodels/bin/python gen_score.py                # ~40 s, seed 7
/home/hackathon/.venvs/abmodels/bin/python build_submission.py         # ~60 s, official checks sharded over 8 processes
cd jobs && python3 submit_japanfold.py rf3_payload.json --go           # spends credits (~$0.35, ~7 min); needs $JAPANFOLD_API_KEY; the idempotency key is the payload hash
mkdir results/rf3_payload/structures && unzip -q results/rf3_payload/archive.zip -x '.parts/*' -d results/rf3_payload/structures     # untrusted download: its own new directory
cd .. && /home/hackathon/.venvs/abmodels/bin/python -I score_structures.py jobs/results/rf3_payload/structures jobs/results/rf3_payload/results.json \
    jobs/targets_manifest.csv out/final_variants_scored.csv /home/hackathon/leo-work/jf/Data/her2_fv_1n8z.cif out/rf3_structure_scores.csv
mkdir eda_current && git show origin/main:NICCOLO/eda/prefilter.py > eda_current/prefilter.py && ln -s ../eda/data eda_current/data    # eda/ = prefilter.py at da498b8 + data symlinks
/home/hackathon/.venvs/abmodels/bin/python rescore_current_prefilter.py && /home/hackathon/.venvs/abmodels/bin/python rank_final.py
bash sync_repo.sh "what changed"                                       # pull + push
```
