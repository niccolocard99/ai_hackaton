# Making the novelty gate fast enough for a 1M-candidate screen

The optimisation-track screen is cheap except for one stage. This is the measurement,
the fix, and the limits.

## 1. Where the time actually goes (measured, not assumed)

Profiled on the team VM, `~/.venvs/abguard/bin/python`, `POLARS_MAX_THREADS=4`, `nice -n 19`
(so a busy machine — the real numbers on an idle VM are better). Candidates are the
`scored_evolution` output of a 50k pilot run; `bench/bench_abguard.py` reproduces it.

| `abguard.annotate` call | Rate | Extrapolated to 1M |
|---|---|---|
| numbering only (`reference=None`, the `abguard_hard` path) | **1,435–1,564 designs/s** | ~11 min |
| numbering **+ novelty gate** (the `track_check` path) | **37.2 designs/s** | **~7.5 h** |

Both are linear in batch size (checked at n=1k/10k/50k and n=100/500/2000). So novelty
is ~97% of the cost: **26 ms/design**, against 0.67 ms for numbering.

The cause is one line, `abguard/novelty.py:375` in `_nearest`:

```python
queries.join(known, how="cross")
```

Every distinct design CDRH3 is compared against **all 610,109** reference loops. The
per-comparison cost is already good (Rust, parallel, ~23M Levenshtein/s) — the
*algorithm* is the problem. Numbering is not: `immunum.polars.number_domains` is a
single vectorised expression over the whole column, which is why it is 40x faster.

## 2. The fix

The gate never needs the true nearest neighbour. It asks one question —
`cdrh3_min_distance >= required_cdrh3_mutations`, for a `required` of 2 or 3
(`MUTATIONS_REQUIRED = {novel: 3, de_novo: 2}`). A bounded question is answerable in
time independent of the reference size.

`abpipe/novelty_index.py` (`NoveltyIndex`):

1. **Deletion-neighbourhood index (SymSpell).** If `levenshtein(q, r) <= d`, the sets of
   subsequences obtained by deleting up to `d` characters from `q` and from `r` must
   intersect. A few hash lookups replace 610,109 comparisons.
2. **Exact verification.** The converse does *not* hold — an intersection only proves
   `distance <= 2d` — so every candidate is verified with a real Levenshtein
   (`rapidfuzz.process.cpdist`, bit-parallel, multi-threaded). This is what makes the
   result exact, and it is also why hashing keys is safe: a collision adds a candidate,
   and verification drops it.
3. **Lazily-built length buckets.** `levenshtein(q, r) >= |len(q) - len(r)|`, so a 13-mer
   query can only match references of length `13 ± d`. Buckets are indexed on first use,
   which keeps the reference's long tail (loops up to 60 residues, where the number of
   deletion combinations explodes) out of the index entirely.

`abpipe/abguard_fast.py` swaps **only** `_nearest`. Numbering, hard checks, the liability
scan and all grading logic are abguard's own, untouched.

```python
from abpipe import abguard_fast
abguard_fast.install()      # idempotent; uninstall() restores the original
```

## 3. Result

Reference: the bundled 610,109 loops. Queries: 1,000,000 trastuzumab CDRH3 variants with
4–8 substitutions in the 10-residue window (998,027 distinct). `bench/bench_novelty_index.py`.

| | Rate | 1M candidates |
|---|---|---|
| abguard cross join | 37.2 designs/s | 7.47 h |
| `NoveltyIndex`, cap 1 | **233,578 designs/s** | **4.3 s** + 3.4 s build |

**~6,300x on the novelty stage.** It is also doing real work: of the 1M variants, 317 are
at distance 0 from a known loop (literal prior art) and 11,771 at distance 1 — 12,088 that
the official rule rejects.

## 4. The other half: abguard runs after ranking

Speeding novelty up left numbering as the top cost (~0.7 ms/design, ~60% of what remained).
The fix there is not a faster algorithm but a better place in the funnel.

abguard is the *authoritative* gate, but only the top of the ranked list can ever be
submitted — 20 designs at most, of which 4 are guaranteed to be tested. Numbering 360,000
survivors to submit 20 is work thrown away. So the abguard modules moved out of `filters`
and into a new `validators` stage that runs **after `rank`**, capped by `validate_top`:

```yaml
filters:                 # O(1) per design, runs on everything
  - framework_integrity
  - window_novelty
  - prefilter_gate
  - dedupe
scorers:
  - niccolo_prefilter
validators:              # authoritative and expensive, after ranking
  - abguard_hard
  - new_liabilities
  - track_check
validate_top: 2000
```

**This is a reordering, not a weakening.** The one real hazard is that `abpipe.iterate`
picks its elites from `03_scored.parquet`, which is *upstream* of the validate stage —
so without care, designs abguard rejects would seed the next round and reach
`selected.csv` with no symptom until submission. `iterate` now restricts selection to
validate survivors, and `tests/test_validate_stage.py::test_the_loop_never_selects_a_design_that_failed_validation`
fails if that is ever undone (it was written first, and did fail).

Effect on a 1M screen, using the measured per-stage rates:

| | abguard in `filters` | abguard in `validators` |
|---|---|---|
| designs numbered | ~359,000 (the gate survivors) | **2,000** |
| numbering time | ~4 min | **~1.4 s** |
| novelty time (patched) | ~2 s | ~0.01 s |

Combined with §2–3, the 1M screen goes from ~7.6 h to **~2.5 min**, and what dominates is
now generation — which is the thing actually producing candidates.

Raise `validate_top` if validation rejects so much that fewer than `--select` designs
survive; the round log and `manifest.json` report `rows_in` and `rows_out` for the stage.

## 5. What is verified, and what is not

**Verified** — `tests/test_novelty_index.py`, 14 tests, all passing (5 min):
- Exact agreement with brute-force rapidfuzz Levenshtein over the **full real reference**
  at `max_distance` 1, 2 and 3 — both the distances and the per-threshold verdicts.
- The official rule (`min_distance >= 2`) on 200 boundary-adjacent queries, asserting the
  test set contains both verdicts so the comparison cannot pass trivially.
- A reported nearest reference really is at the reported distance.
- Queries are deliberately boundary-adjacent: real reference loops mutated by 0–5
  substitutions/insertions/deletions, plus random noise.

**Not yet verified** — `tests/test_abguard_fast.py` is written but has **not been run**. It
needs abguard (a Linux-only compiled extension; `immunum` ships a 17 MB `.so`), so it
cannot run on a laptop. Until it runs on the VM, *"patched abguard gives identical
verdicts"* is a well-founded expectation, not a measured fact. Run it with:

```bash
cd ~/nick-1m/pipeline     # isolated rsync copy; does not touch the shared repo
ABPIPE_DATA_DIR=~/hackathon-team-3/apheris/hackathon-data POLARS_MAX_THREADS=4 \
  nice -n 19 ~/.venvs/abpipe/bin/python -m pytest tests/test_abguard_fast.py -q
```

## 6. Limits (do not read past these)

- **Distances saturate at the cap.** Built at `max_distance = max(MUTATIONS_REQUIRED) = 3`,
  every *gate* comparison is exact and so is every *rejection message* (a failing design is
  within the cap by construction). What is lost is the true distance of designs that
  comfortably pass — reported as 4 rather than, say, 7 — and their nearest-neighbour
  fields. The gate never uses that. `ABGUARD_FAST_MAX_DISTANCE` raises the cap; cost grows
  with it.
- **Cost grows with the cap, not the reference.** Cap 1 indexes 14 keys per loop; cap 3
  indexes 378. Cap 3 over 610k loops is ~230M keys — fine on the VM, heavy on a laptop.
- **This is a speed fix, not a science fix.** It changes nothing about what the gate
  *means*. Passing novelty still says nothing about binding, epitope preservation or
  thermostability (`NICCOLO/eda/AGENT_GUIDE.md` §7).
- **abguard's own build is the authority.** This patch is for *screening at volume*. The
  shortlist should still be gated by unpatched `abguard` / `NICCOLO/scripts/check_*.py`
  before submission — that run is small, so its cost does not matter.
- `window_novelty` in the loop config is a *substitution-only* pre-screen and is already
  exact for its rule; `scored_evolution` enforces `min_dist_provided: 3` at generation,
  which is why both dropped 0 rows in the 50k pilot. This index is the Levenshtein-correct
  gate that runs *after* them, at the same speed.
