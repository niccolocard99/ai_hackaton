# De novo gate: what the official checker actually measures

Everything below marked **measured** was read out of the installed code or produced by
running it on 2026-11-14. Nothing here is inferred from the hackathon docs.

## How to run it

```bash
# official verdict only
/home/hackathon/.venvs/abguard/bin/python \
  /home/hackathon/hackathon-team-3/NICCOLO/scripts/check_denovo.py \
  designs.csv -o report.csv --strict-exit          # exit 0 only if every row PASSes

# official verdict + our margin columns
/home/hackathon/.venvs/abmodels/bin/python \
  /home/hackathon/hackathon-team-3/NICCOLO/denovo/gate.py designs.csv -o out.csv
```

`gate.py` needs pandas + numpy; the official checker needs abguard. No venv on this box
has both, so `gate.py` resolves an abguard interpreter itself and reaches it only through
`subprocess`. Run as a CLI from an interpreter without pandas it re-execs into one that
has it (`PANDAS_PY` / `ABGUARD_PY` override the search).

## The reference CDRH3 set: 610,109 entries

**Measured.** File:
`/home/hackathon/.venvs/abguard/lib/python3.12/site-packages/abguard/reference_data/abguard_cdrh3.parquet`
(2,028,925 bytes; a 1-table Parquet bundle with `reference_data/manifest.json` beside it,
located at runtime by `abguard.reference.bundled_reference()`).

- 610,109 rows, columns `cdrh3, source, id`.
- All 610,109 `cdrh3` values are **distinct** — so 610,109 is both the row count and the
  number of unique loops the gate compares against. The manifest's `n_rows` agrees.
- CDRH3 length: min 1, median 13, mean 13.2, max 60.
- Sources: `tz_her2_facs_3classes_524346` 524,308 · `uspto` 38,308 · `plabdab` 38,197 ·
  `pdb` 5,943 · `plabdab_nano` 1,123 · `trastuzumab_her2_spr_all_768` 751 ·
  `proteinbase_public` 571 · `target_tfr1` 391 · `target_tnfa` 373 · `swissprot` 108 ·
  `therasabdab` 34 · `thpdb` 2.

**The brief's two candidate counts are both wrong as a description of the reference set.**
610,109 is the real count. 525,063 is close to the two trastuzumab assay sources added
together (524,308 + 751 = 525,059) and is presumably a pre-dedup count of those tables
alone, not of the reference.

Every reference loop is **re-derived** from the depositor's full sequence by
`abguard.reference.cdrh3_from_sequences`, so both sides of every distance are immunum
IMGT CDR3 regions, never a depositor's own CDRH3 definition.

## How identity is computed, and against what

**Measured**, `abguard/novelty.py::_identity_to_trastuzumab`.

It is **only ever compared to trastuzumab** — confirmed, not assumed. The only references
in the expression are the module constants `TRASTUZUMAB_VH` (120 aa) and `TRASTUZUMAB_VL`
(107 aa, DrugBank DB00072). No other antibody takes part in the identity step. The
reference database is used *exclusively* for the CDRH3 distance, never for identity.

The number is **not** an alignment identity. It is one pooled edit ratio:

```
identity = 1 - (d_H + d_L) / (max(|VH|, |ref VH|) + max(|VL|, |ref VL|))
```

- `d_H`, `d_L` are **Levenshtein** distances (`polars-distance`).
- Each chain is first reduced to its **immunum IMGT Fv**: `fr1 cdr1 fr2 cdr2 fr3 cdr3 fr4`
  concatenated, with `prefix`/`postfix` dropped — which is what makes identity blind to an
  scFv linker or a purification tag. The trastuzumab side is segmented by the same code
  (`_reference_fv`), not written down, so the two sides always use the same boundaries.
- A `vhh` has no light-chain domain, so the light terms drop out of **both** numerator and
  denominator and its identity is decided purely by its VH. Brief correct.
- For `paired` and `scfv` both chains count, length-weighted. For an `scfv` the VL is found
  as the second numbered domain inside `heavy_sequence`, so a whole-scFv-in-heavy row is
  still pooled over two chains. Brief correct.
- Chains are picked as the *first* domain of chain class `H` and the first of `K`/`L`, from
  the pooled domain list across whichever sequence columns exist.
- Null identity (no numbered domain at all) is excluded from the class split; such a row has
  already failed hard validation.

Class split (`abguard.novelty.NoveltyClass`): `identity >= 0.70` → `novel`,
`identity < 0.70` → `de_novo`. **Measured boundary behaviour:** a design at exactly
0.700000 was classed `novel` and FAILED, so 0.70 is a hard exclusive ceiling.

`check_denovo.py` rejects anything not classed `de_novo` outright — it does not reclassify
it onto the optimisation track. Note the failure text says "identity to a known antibody",
which is misleading: it is identity to trastuzumab and nothing else.

## CDRH3 novelty: the exact metric

**Measured.** **Levenshtein** (not Hamming), via `polars-distance`, between the design's
immunum IMGT CDRH3 and every reference CDRH3, reported as the single minimum
`cdrh3_min_distance`.

- Trastuzumab's own CDRH3 is simply **prepended as one more reference row**
  (`trastuzumab_reference()`), so "at least N from trastuzumab AND from every reference
  CDRH3" collapses into one minimum. The brief's two-part phrasing is behaviourally right
  but there is only one comparison in the code.
- **Measured redundancy:** `SRWGGDGFYAMDY` is *already in the bundle* as `pdb:1fvc_B`, so
  the trastuzumab row adds nothing when the bundle is used.
- Required minimum, `abguard.novelty.MUTATIONS_REQUIRED`: `de_novo` → **2**,
  `novel` → 3. `check_denovo.py` hard-codes its own `MIN_CDRH3_MUTATIONS = 2` and only
  reaches the CDRH3 test for rows already classed `de_novo`.
- Trastuzumab's IMGT CDRH3 is `SRWGGDGFYAMDY`, 13 aa — **measured** via
  `trastuzumab_cdrh3()`. Note the HER2 assay tables name the 10-mer core `WGGDGFYAMD`
  instead; abguard deliberately does not use that.

## Shared hard checks

**Measured**, `abguard.config.ErrorKind` / `abguard.annotate.with_errors`. A row with any of
these fails before the novelty rules are even considered, reported as `hard check: <kind>`:

`invalid_amino_acids`, `sequence_too_long` (`MAX_SEQUENCE_LENGTH = 250`),
`unpaired_cysteine`, `numbering_failed` (`MIN_CONFIDENCE = 0.1`, scheme `imgt`),
`wrong_domain_count`, `missing_sequence`, `unexpected_sequence`, `duplicate_sequence_id`,
`unknown_kind`, `insufficient_cdrh3_novelty`.

## Report field names

**Measured** from `_track_common.finish`. `check_denovo.py -o report.csv` writes exactly:

```
sequence_id, kind, track_assigned, identity_to_trastuzumab, cdrh3,
cdrh3_min_distance, cdrh3_identity_to_trastuzumab, passes, failures
```

- **`passes`** is the boolean verdict — in CSV it is the lowercase string `true` / `false`.
- **`failures`** is a list column; in CSV the reasons are joined with `"; "`, empty on pass.
- `track_assigned` is `novelty_class`, i.e. `de_novo` or `novel`.
- `cdrh3_identity_to_trastuzumab` is a *different* quantity from `identity_to_trastuzumab`:
  it is `1 - lev(cdrh3, SRWGGDGFYAMDY) / max(len)`, CDRH3 only, added by the track script.
- stdout prints `PASS`/`FAIL` per row plus `N/M design(s) eligible`; `--strict-exit` returns
  1 unless N == M, and the report is written **before** that exit, so a failing run still
  produces a readable report.

## Fixture results, as measured

`fixtures/fixture_designs.csv` → `fixtures/fixture_official_report.csv` (official, exit **1**)
and `fixtures/fixture_gate_report.csv` (enriched).

| sequence_id | kind | official | track | identity | cdrh3 | len | d(tz) | d(min ref) | internal |
|---|---|---|---|---|---|---|---|---|---|
| `fx_tz_wt_paired` | paired | **FAIL** | novel | 1.000000 | SRWGGDGFYAMDY | 13 | 0 | 0 | FAIL |
| `fx_tz_cdrh3_2sub` | paired | **FAIL** | novel | 0.991189 | SRWGGDGFYKLDY | 13 | 2 | 1 | FAIL |
| `fx_vhh_denovo` | vhh | **PASS** | de_novo | 0.611940 | AAQLPHNTRIEKSGWVDFYTGNQLPSHE | 28 | 24 | 17 | PASS |

Both failures report the identity rule, not the CDRH3 rule:
`not de novo: 100% / 99% identity to a known antibody (needs < 70%); it would be an
optimisation design and is rejected here, not reclassified`.

Note `fx_tz_cdrh3_2sub`: its CDRH3 is 2 edits from trastuzumab's but only **1** edit from
the nearest *reference* loop — the 524k trastuzumab HER2 assay variants in the bundle cover
single-mutant neighbourhoods of `SRWGGDGFYAMDY` densely. A CDRH3 graft on trastuzumab's
framework is therefore doubly dead: identity and prior art.

### Where the fixture sequences came from

- Trastuzumab VH (120 aa) and VL (107 aa) were extracted from
  `apheris/rfd4/Data/12-her2/her2_fv_1n8z.cif` chains **B** and **A**. **Measured:** they
  are byte-for-byte identical to abguard's hardcoded `TRASTUZUMAB_VH` (120 aa) and
  `TRASTUZUMAB_VL` (107 aa), which is why the WT row measures identity exactly 1.000000.
- `fx_tz_cdrh3_2sub` is that VH with `SRWGGDGFYAMDY` → `SRWGGDGFYKLDY` (A10K, M11L), WT VL.
- `fx_vhh_denovo` is the caplacizumab VHH **framework** from
  `apheris/rfd4/Data/03-vhh/pmvk_caplacizumab.cif.zst` with synthetic non-caplacizumab
  CDR1/CDR2/CDR3. FR1 `EVQLVESGGGLVQPGGSLRLSC`, FR2 `WFRQAPGKGRELVA`,
  FR3 `YPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYC`, FR4 `WGQGTQVTVS`;
  CDR1 `AASGRTFSTYTMG`, CDR2 `AISWSGGHTA`, CDR3 `AAQLPHNTRIEKSGWVDFYTGNQLPSHE` (134 aa).
  CDR3 length was chosen by measurement, not taste: at CDR3 lengths 15/22/24/26/**28**/30
  identity measures 0.653 / 0.633 / 0.631 / 0.621 / **0.612** / 0.603 — a short CDR3 clears
  the official 0.70 bar easily but sits above our internal 0.62 ceiling.

## `gate.py`

`run_gate(designs, *, output=None, reference_db=None, abguard_py=None, checker=None,
identity_max=0.62, cdrh3_distance_min=4) -> pandas.DataFrame`. `designs` is a submission
path or a DataFrame in that shape. CLI: `python gate.py designs.csv -o out.csv`
(`--strict-exit` gates on our stricter bar, not the official one).

Columns: `sequence_id, kind, official_pass, official_track_assigned, official_reason,
identity_to_trastuzumab, margin_vs_identity_cap, cdrh3, cdrh3_len,
cdrh3_dist_trastuzumab, cdrh3_dist_min_reference, official_cdrh3_min_distance,
internal_pass, internal_reason, design_description`. `df.attrs` additionally carries
`official_exit_code`, `official_stdout`, `reference_path`, `reference_entries`,
`abguard_python`.

`internal_pass` = `official_pass` AND `identity <= 0.62` AND `cdrh3_dist_trastuzumab >= 4`
AND `cdrh3_dist_min_reference >= 4`.

`cdrh3_dist_min_reference` is recomputed independently of the official run (numpy only, no
abguard) so it can be used as a design-time filter and as a cross-check. A composition
lower bound — an edit changes the amino-acid multiset by at most one, so
`lev >= max(excess_ref, excess_query) >= |len difference|` — screens all 610,109 entries in
one vectorised pass; buckets of equal bound are then aligned in ascending order and the
scan stops as soon as the bound exceeds the best exact distance found. The result is the
exact minimum, not an approximation.

**Verified:** on 200 independently generated VHH candidates, `cdrh3_dist_min_reference`
matched the official `cdrh3_min_distance` (polars-distance over the full cross join) on
**200/200 rows**, distances spanning 4–18.

**Measured speed:** 200 designs in 52 s end to end (~0.25 s/design, one subprocess call for
the whole batch); 3 designs in 3.6 s warm, 5.9 s on the first-ever run including cache
build. CPU only.

`cache/reference_cdrh3_v1.npz` (7.5 MB) is derived data — a numpy-readable copy of the
reference's CDRH3 column, built on first use because the only interpreter that can read
abguard's Parquet is the one without pandas. Delete it and it rebuilds; it is gitignored.
`DENOVO_GATE_CACHE` relocates it, `--refresh-reference` rebuilds it.

## Premises in the brief that turned out wrong

1. **Reference count.** Neither 525,063 nor the doc's alternative describes the set. It is
   **610,109** unique CDRH3s in `abguard/reference_data/abguard_cdrh3.parquet`. (The brief
   asked which of the two was right; the answer is the larger one, 610,109.)
2. **`pmvk_caplacizumab.cif.zst` does not contain caplacizumab's sequence.** It is an
   RFD4-Proteina *design input* (`data_pmvk_3CH4__VHH_caplacizumab_cdr123.cif.zst`) in which
   **CDR1, CDR2 and CDR3 are masked as `X`**: chain H is 83 framework residues only, and
   chain A is PMVK (phosphomevalonate kinase, 128 aa), not an antibody at all. The full
   caplacizumab VHH is **not** recoverable from it. Only the framework was used, which is
   what option (iii) allowed.
3. **Identity is not an alignment identity** and is not computed "against the nearest known
   antibody" despite the failure message saying so — it is a single pooled Levenshtein ratio
   against trastuzumab's two variable domains and nothing else.
4. **The two CDRH3 rules are one rule.** Trastuzumab's CDRH3 is just another reference row,
   and it is already present in the bundle (`pdb:1fvc_B`), so "from trastuzumab AND from
   every reference CDRH3" is a single minimum with a redundant first term.
5. Minor: the brief implies trastuzumab's real VH/VL had to be recovered from the CIF.
   They are also hardcoded in `abguard/novelty.py` as `TRASTUZUMAB_VH`/`TRASTUZUMAB_VL`,
   and the CIF chains B/A are **byte-for-byte identical** to them (measured) — so the CIF
   extraction was a confirmation, not a necessity.
6. Everything else in the brief checked out: the CLI invocation, the input schema,
   `SRWGGDGFYAMDY` (13 aa), VH-only identity for a VHH, pooled identity for paired/scFv,
   the `< 0.70` cap, the `>= 2` de novo CDRH3 floor, length `<= 250`, the canonical-20 and
   even-cysteine checks, and `min_confidence = 0.1`.
