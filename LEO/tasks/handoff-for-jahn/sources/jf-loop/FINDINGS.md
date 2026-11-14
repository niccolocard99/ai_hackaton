# Findings and evidence

**As of:** 2026-11-14, about 16:05 UTC · **Written for:** whoever must defend or extend these results (the pitch scores judgment at 40%: this is the reasoning, with its limits stated). Numbers are traceable to files in [snapshot/](snapshot/) or on the VM; each section says which. Overview: [HANDOFF.md](HANDOFF.md). Operations: [RUNBOOK.md](RUNBOOK.md). The first run's results are in [RESULTS.md](RESULTS.md) (its section 5 was rewritten after an error).

## 1. The measured facts the design rests on

SPR-768 table, restricted to the **legal band** (at least 4 substitutions from trastuzumab's 10-mer; n = 556, 24% have a measurable KD). "Measurable" means a quantified KD rather than no reportable KD. Source: `prefilter/` analyses and `final/` inputs.

| Position | Parent kept: n, measurable | Substituted: n, measurable | Lift |
|---|---|---|---|
| **Y105** | 135, 56% | 421, 13% | **4.2x** |
| W99 | 112, 54% | 444, 16% | 3.4x |
| F104 | 74, 59% | 482, 18% | 3.3x |
| G103 | 242, 37% | 314, 13% | 2.7x |
| G100 | 216, 35% | 340, 16% | 2.2x |
| M107 | 26, 42% | 530, 23% | 1.9x |
| A106 | 145, 31% | 411, 21% | 1.5x |
| D108 | 302, 25% | 254, 21% | 1.2x |
| G101 | 290, 22% | 266, 26% | 0.8x |
| D102 | 90, 20% | 466, 24% | 0.8x |

Keeping Y105 holds for 8 of the 9 legal-band variants measured below 1 nM; keeping **both** W99 and Y105: n = 40, 88% measurable, median 17 nM, 7 of 9 sub-nM. **This is association, not causation**: the library was built by someone else's design process, so position and context are confounded. It is the best evidence available, and it is why Y105 is pinned.

## 2. The rules as measured

- abguard's CDRH3 is the 13-mer `SRWGGDGFYAMDY`: positions 97-109 in the sequential numbering used throughout this task (the IMGT CDR3). Under 70% identity needs **at least 4 substitutions** (9/13 = 69.2%; 3 give 76.9%).
- Novelty: at least 3 from **every** reference loop, by Levenshtein (indels count). The reference is the HER2 tables' 525,063 CDRH3s plus prior art (610,109 loops, per the `abguard_fast.py` docstring). It is the binding constraint, not identity.
- My first four geometry-driven designs passed identity at 69.2% and **failed novelty at distance 1** (`RESULTS.md`).
- Positions 97, 98 and 109 are inside abguard's CDRH3 and **constant across the library**.

## 3. Arms and yields (JapanFold pool at about 16:05 UTC; 1,432 designs, 393 distinct, 262 eligible)

Computed from `snapshot/pools/jf_pool.csv` (refresh with the RUNBOOK commands; the loop adds about 300 designs per cycle).

| Arm | Designs | Distinct | Distinct % | abguard pass | Eligible | Margin (>=4) | Best p_binder |
|---|---|---|---|---|---|---|---|
| `jy` (Y105) | 458 | 112 | 24% | 101 | 97 | 5 | 0.942 |
| `jy97` | 224 | 74 | 33% | 57 | 53 | 3 | 0.906 |
| `jw97` (W99+Y105, S97 designed) | 227 | 49 | 22% | 30 | 29 | 0 | 0.990 |
| `jfree` (control) | 120 | 54 | 45% | 46 | 38 | 24 | 0.972 |
| `jfree97` (control, S97 designed) | 32 | 28 | 88% | 25 | 19 | 10 | 0.980 |
| `jygf` | 80 | 15 | 19% | 9 | 8 | 0 | 0.964 |
| `jyf97` | 45 | 14 | 31% | 8 | 7 | 1 | 0.982 |
| `jym` | 22 | 9 | 41% | 9 | 8 | 3 | 0.671 |
| `jyf` | 18 | 12 | 67% | 4 | 3 | 1 | 0.974 |
| `jwy` (W99+Y105) | 37 | 17 | 46% | **0** | 0 | 0 | 0.983 |
| `jan4` (W99+G103+F104+Y105) | 164 | **5** | **3%** | **0** | 0 | 0 | 0.984 |
| `jyg`, `jygm` | 2, 3 | 1, 3 | | 0, 1 | 0 | 0 | |

What this says:
- **Generic RFD4 with W99 and Y105 both pinned and S97 left alone gives zero legal designs in 201** (`jwy`, `jan4`), and collapses (`jan4`: 164 designs, 5 distinct). The sequences it returns sit inside the dense part of the 525k library.
- **Designing S97 changes that** (`jw97`: 29 legal designs that all keep W99 and Y105). S97 is +1 novelty distance from every known variant for free.
- The Y105-only arms are productive but mode-collapse (`jy`: 24% distinct).
- The `jfree` controls have the best margin rates (24 of 38, 10 of 19) and mostly lose Y105 (section 4). That the two go together is shown; that one causes the other is not.

## 4. Keeping Y105 and having novelty margin pull against each other

296 eligible designs (JapanFold pool plus the frozen H100 pool), `min_dist` = abguard's novelty distance:

| | margin (>= 4) | edge (= 3) |
|---|---|---|
| Y105 kept (240) | 19 (8%) | 221 |
| Y105 lost (56) | 34 (61%) | 22 |

`min_dist` over all 296: **3 for 243, 4 for 51, 5 for 2.** The edge designs pass the local abguard exactly at the threshold. An earlier version of this claim ("margin and bindability are anti-correlated") was wrong: it rested on a misread column (section 10, item 3). The relationship that holds is **Y105 retention vs margin**.

## 5. The S97 lever

- In the 1N8Z parent complex S97 makes **0 heavy-atom contacts to HER2** within 4.5 A (R98: 0, Y105: 32, CDRH3 total 53). R98's salt bridge to HER2 D585 is not visible in 1N8Z because that structure does not model HER2 581-590 (it shows in 8PWH and 8Q6J); R98 is held fixed anyway.
- **Hedge test** (`snapshot/final/s97_hedge.csv`, from `compare_lists.py`): 95 variants of 24 top parents (S97A, G, T, V). **All 95 stay legal; the novelty distance rises by +1 in 92 of 95** (the 3 exceptions already differed at 97). For **23 of 24 parents some S97 substitution reaches margin >= 4.**
- The model itself chose Ala at 97 when allowed (the `jw97` design `WGSTSHYAFD` was `S97A`). Ser to Ala at that position (Kabat 93) is a common germline-like change in general antibody knowledge; **that is not measured here**.
- **Unknown, and the reason this is a decision for a human:** every variant in the SPR and FACS data keeps S97, so there is **no measurement** of whether S97X changes binding or stability.

## 6. Three lists on one yardstick (at 15:30 UTC; **predates the S97 arms**)

`compare_lists.py` recomputes everything from the sequences. "mine" = the then top 20 of this task; the other two are the sibling `design-loop` session's lists (`out/portfolio.csv`, `out/FINAL.csv`; the latter carries VL D28K and framework add-ons). Source: `snapshot/final/comparison.csv`.

| List | n | legal | at distance 3 | margin | W99 | G103 | F104 | Y105 | M107 | anchors kept of (W99, G103, F104, Y105) | median p_binder | near binders | CDRH3 liabilities vs parent |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mine | 20 | 20 | 19 | 1 | 5 | 1 | 1 | 20 | 9 | **1.35** | 0.80 | 15 | -1.95 |
| sibling portfolio | 20 | 20 | 20 | 0 | 18 | 17 | 20 | 20 | 20 | **3.75** | 0.94 | 20 | -0.65 |
| sibling final | 20 | 20 | 20 | 0 | 20 | 18 | 20 | 20 | 20 | **3.90** | 0.94 | 20 | -0.50 |

At that time the sibling's designs kept far more of the measured anchors and scored higher on the classifier; both sets sat at the novelty edge. The refreshed top 20 (HANDOFF) has since moved toward anchor-rich designs, so **re-run `compare_lists.py` before comparing**. The contact and epitope columns in `comparison.csv` are model-dependent (mine from RFD4's own complex, the sibling's from its Boltz-2 on the H100) and must not be compared across lists.

## 7. Where the current best designs come from

See the HANDOFF table. The top 3 (`WGLYGFYEFD`, `WGPAGFYEFD`, `WGPHGFYEFD`) keep W99, G103, F104 and Y105, have 4 substitutions (the legal minimum), score p_binder 0.97-0.98, dock at 1.15-1.28x the parent's CDRH3-HER2 contacts, and sit at novelty distance 3. Only 3 of 240 Y105-keeping eligible designs keep all four anchors.

## 8. Tooling checks, throughput, cost

- **Fast checker** (`fast_check.py`, rapidfuzz Levenshtein): on 60 designs, 32 of them rejected, matched the standard checker's stored verdicts on passes, novelty distance and identity, **60/60**, and rejected the same 32; 2.65 s against about 90 s. (Whole-sequence distance to the provided variants was the dominant cost.)
- **Patched abguard** (`abguard_fast.py`, Niccolo): not measured here. Its docstring claims identical verdicts and a 6,300x faster novelty stage, with distance capped at 4.
- **Throughput and cost, measured:** 38 jobs (304 designs) per cycle, about 6.5 min, balance $1,095.28 to $1,094.70 over the first full cycle = **$0.087/min**, about 47 designs/min. The shared H100 gave about 25 s per fresh design and 9 s per refined one.
- **VM load:** about 70 on 16 cores before the throttling (unthrottled 86-thread checkers), about 2.4-3 after.
- **Ledger:** 234 runs and 1,262 variant rows on the VM copy (230 runs mine). All runs launched before the rule was seen are logged retroactively and flagged.

## 9. Caveats on every metric (read before quoting anything)

- **`p_binder`** is a logistic classifier on one-hot CDRH3 (AUROC 0.913 inside the legal band, grouped CV, NICCOLO). It is a ranking, **not a calibrated probability for novel designs**. My calibration table (in-sample) was flattering: bins at 0.8-0.9 had 100% measurable binders, median 19 nM, and none beat trastuzumab.
- **The affinity ranker is blind where we are scored.** Its predicted floor is 1.135 nM against trastuzumab's measured 0.029-0.178 nM (8 replicates, 6.2x spread). Spearman by regime (grouped CV): full range +0.687, legal band +0.628, below 10 nM +0.415, below 1 nM **+0.167**, below 1 nM and legal band **-0.283** (n = 9). Median error 17.7x below 1 nM (`prefilter/subnm_resolution.py`). It interpolates only within about 4 substitutions of a measured variant (NICCOLO's AGENT_GUIDE section 3). **No `pred_kd` here is a KD.**
- **`n_bind_le4` and `best_kd_le4`** count measured binders within 4 substitutions. A design with the legal-minimum 4 substitutions has many near-parent measured variants within reach, so large counts mostly reflect *closeness to trastuzumab*, and `best_kd_le4` is often the same 0.023 nM variant (3 substitutions from trastuzumab). **It is locality, not quality.**
- **`contact_ratio`, `epi_overlap`, `y105_patch`** describe the design model's own complex. They are not validated against binding, and fold models rank binders at chance on this task (NICCOLO: Boltz-2 and Protenix, n = 180). `epi_overlap` is an unvalidated proxy for epitope preservation, which has no labelled data.
- **Thermostability is not predicted.** NICCOLO's NanoMelt analysis found CDRH3 variants cannot move predicted Tm beyond its 4.1 C error. The only developability signal here is abguard's liability scan.
- **Base rate for the whole approach** (sibling's held-out backtest): the ranker's top legal-band pick beats trastuzumab's KD about 1% of the time and lands within 3x about 14% of the time. **Nothing here is predicted to beat trastuzumab.**
- **The novelty margin** is zero for almost every top design (section 4), and the organisers' reference may differ from the local one.

## 10. Mistakes made in this session (corrected; do not repeat)

1. **Launched generation runs before reading the updated AGENTS.md** (two H100 runs and a JapanFold loop). Fixed: ledger protocol in both loops, 49 runs backfilled and flagged RETROACTIVE.
2. **Misread "use api for japanfold" as "in addition to the H100".** The user corrected it ("i told you use api not h100"); the H100 loop was stopped at 14:49 UTC.
3. **Read `nn_dist` (distance to the nearest *measured* variant) as `min_dist` (novelty margin)** and wrote that margin and bindability were anti-correlated. Rewritten in `RESULTS.md` section 5.
4. **Said "W99 + Y105 gives zero legal designs".** True only for JapanFold's generic RFD4 without S97 designed; the H100 `her2-novel` weights and the S97 arms produce legal designs.
5. **An `IncompleteRead` escaped the download handler** and left 5 ledger runs `started`. Fixed (retry, zip test, `reconcile()`).
6. **Ran unthrottled validation** (86 threads) on the shared VM, adding to a load of about 70. Fixed (fast checker, thread limits, `nice`).
7. **Self-matching process patterns** (`pgrep -f`, `pkill -f`) twice: a false "still running", and a status script that killed itself.
8. **Read "the Y" as Y105 without confirming.**
9. Earlier-phase corrections (literature review: a free-set claim that did not survive 8 structures, an 83 C Tm that was a citation) are in `LEO/tasks/litreview-optimisation/AUTORESEARCH_BRIEF.txt`.

## 11. Which file backs which claim

| Claim | File |
|---|---|
| Position lifts, W99+Y105 | `prefilter/subnm_resolution.py` family and `final/` inputs; recompute from `~/leo-work/prefilter/data/trastuzumab_her2_spr_all_768.csv` |
| Arm yields, Y105 vs margin | `snapshot/pools/jf_pool.csv`, `h100_pool.csv` |
| Shortlist, tiers, reasoning per design | `snapshot/final/candidates.csv` |
| BASE VALIDATION of the top 20 | `snapshot/final/reports/optimisation.csv`, `denovo.csv`, `liabilities.json` |
| Three-list comparison, S97 hedge | `snapshot/final/comparison.csv`, `s97_hedge.csv` |
| Every run, seed, contig, recipe, ledger id | `snapshot/jobs/`, `snapshot/ledger-vm/` |
| Per-ingest abguard verdicts | `snapshot/abguard/` |
| The structures of the top 20 | `snapshot/final/cif/` |
