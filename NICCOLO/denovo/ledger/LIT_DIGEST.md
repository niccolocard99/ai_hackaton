# Literature digest — de novo VHH binders to HER2 domain IV

One section per design decision we actually have to make. Each gives the
evidence, then the number or rule we are adopting.

**Every number is marked `CITED` (with a PMID or a file path) or
`OUR JUDGEMENT`. These are never blurred.** Where a threshold is ours, the
reasoning is given so it can be argued with, and the queued AMASS query that
would replace it with a cited one is named in
[`research/pass-001.md`](research/pass-001.md).

Sourced **cache-only**: `AMASS_API_KEY` was unset and credits unclaimed, so
zero AMASS network calls were made. Records come from the 32 cached responses
in [`cache/`](cache/) (585 unique records). Platform numbers come from the
Apheris docs under `/home/hackathon/hackathon-team-3/apheris/`.

---

## 1. CDR3 length range for a de novo VHH binder

### Evidence

| Source | What it says | Status |
|---|---|---|
| `apheris/rfd4/Data/03-vhh/pmvk_caplacizumab.cif.zst` | The supplied VHH example's own ranges, read from `annotation_expsegmin/max_1_residue`: **CDR1 13-13, CDR2 9-10, CDR3 12-21** | CITED (file) |
| PMID:19062174 (Proteins 2009, 188 cites) | Median CDR-H3 backbone RMSD of prediction by loop length: **1.6 Å (4-6 res), 1.9 Å (7-9), 2.4 Å (10-11), 3.1 Å (12-14), 6.0 Å (17-22)** | CITED |
| PMID:23255066 (Proteins 2013) | Ab initio H3 prediction benchmarked over 53 loops, length **4-22** residues | CITED |
| PMID:8918938 (JMB 1996) | The 1230-fold-matured anti-c-erbB-2 scFv had a **20 aa** VH CDR3 ("Due to its length of 20 amino acids, four V(H) CDR3 libraries were constructed") | CITED |
| PMID:25484051 (mAbs 2015) | Crystal structure with a **24 aa** protruding VH CDR3, "critical to the predicted binding mode" | CITED |
| PMID:18687338 (JMB 2008, 75 cites) | A soluble, stable human VH domain library spanned CDR3 **7 to 24 residues**, all purifiable from E. coli periplasm | CITED |
| PMID:41193805 (Nature 2025) | De novo VHHs bound 4 epitopes with cryo-EM-confirmed CDR accuracy — but **no length numbers in the abstract**, and the cache has no fulltext | CITED (existence only) |
| `apheris/rfd4/` (whole doc set) | **No allowed min/max CDR length is stated anywhere.** Only per-example values. | CITED (negative) |

### Decision

**Design CDR3 at 12-16 residues; keep CDR1 at 13 and CDR2 at 9-10 as supplied.**
Carry **one control arm at 17-21** so the ceiling is tested, not assumed.

- Floor of **12**: CITED, the platform example's own lower bound.
- Ceiling of **16**: **OUR JUDGEMENT.** Reasoning: our entire funnel ranks on
  predicted-structure quantities (min_ipae, iptm, binder_plddt). PMID:19062174
  shows H3 prediction error roughly doubles from 3.1 Å at 12-14 residues to
  6.0 Å at 17-22. A 6 Å backbone error makes an interface-PAE ranking
  meaningless, so designing past ~16 means generating molecules we cannot
  triage. Note honestly that this is a 2009 homology-modelling benchmark
  extrapolated to a 2026 co-folding model — hence the control arm.
- The counter-evidence is real and worth stating in the pitch: 20 aa and 24 aa
  CDR3s demonstrably bind, and 7-24 aa VH domains express fine. **Our bound is
  about our ability to rank, not about what can bind.**
- Queued Q2 would give a natural camelid CDR3 length prior, which we lack.

---

## 2. How many CDRs to design — CDR3 only, or all three?

### Evidence

| Source | What it says | Status |
|---|---|---|
| `apheris/rfd4/02-miniprotein.md:30-31` | "Antibody examples keep framework sequence and internal geometry, generating **three CDRs for VHH** or six for Fv." | CITED (file) |
| `apheris/rfd4/15-antibody-design.md:62` | "A VHH spec has one antibody chain and **three CDRs**." | CITED (file) |
| `apheris/rfd4/12-her2.md:127-128` | "**Only CDRH3 changes.** To change other CDRs as well, redesign them with the antibody model." | CITED (file) |
| PMID:41193805 (Nature 2025) | For an scFv, cryo-EM verified "the atomically accurate design of the conformations of **all six CDR loops**" | CITED |
| PMID:8918938 (JMB 1996) | Sequential **VL-CDR3 then VH-CDR3** diversification gave the 1230-fold gain to KD 1.3e-11 M — multi-CDR, not CDR3 alone | CITED |

### Decision

**Design all three CDRs.** Not a judgement call so much as a structural fact of
the route: the HER2 fine-tuned weights (`her2-2026-11-07`,
`her2-novel-2026-11-09`) only redesign trastuzumab's CDRH3 **on the Fv**, and
there is **no HER2 VHH weight and no HER2 VHH example** on this platform. A de
novo VHH therefore has to run on `antibody-2026-10-30`, whose VHH mode
generates three CDRs by construction (CITED). CDR3-only is not an option we
have for this arm.

Supporting rather than deciding: PMID:41193805 designed all six CDRs of an scFv
to atomic accuracy, and PMID:8918938's largest affinity gain came from
diversifying two CDR3s in sequence, so multi-CDR design is well precedented.

---

## 3. Hotspot count, and whether multi-patch beats single-patch

### Evidence

| Source | Hotspot atoms used | Status |
|---|---|---|
| `apheris/rfd4/Data/12-her2/her2_cdrh3.yaml` | **4 atoms across 2 chains**: `sel('D/*/593/NZ')`, `sel('C/*/573/CD1')`, `sel('D/*/591/CG2')`, `sel('C/*/572/CD')` — "C_HOT marks the HER2 atoms nearest CDRH3 in **PDB 1N8Z** as contact hotspots" | CITED (file) |
| `apheris/rfd4/02-miniprotein.md:35-39` | Miniprotein example: **3 atoms** | CITED (file) |
| `apheris/rfd4/Data/03-fv/design.json` | Fv example: **2 atoms** | CITED (file) |
| `apheris/rfd4/Data/03-vhh/*.cif` | VHH example: **1 atom** (`condition_hotspot_1_atom`, chain A, ARG 132, CB) | CITED (file) |
| `apheris/rfd4/15-antibody-design.md:30` | "**The epitope**: a few antigen atoms the antibody should contact (hotspots)." | CITED (file) |
| `apheris/rfd4/09-parameters.md` | `resampleHotspots` default false; true "replaces annotated hotspots with **2-3 sampled nearby atoms**". "CDR lengths and hotspot selections are **input conditions, not workflow flags**." Do not combine partial diffusion with hotspot resampling. | CITED (file) |
| PMID:12610629, PMID:31042744 | The 1N8Z and cryo-EM structures give the contact set — **contacts, not per-residue energetics** | CITED |
| Whole Apheris doc set | **The term "multi-patch" appears in no file.** No hotspot count limit documented. | CITED (negative) |

### Decision

**Use the 4 `C_HOT` atoms from the supplied HER2 spec. Single patch.**

- **4 atoms, those 4**: CITED. They are the only epitope definition in the repo
  anchored to a real co-crystal (1N8Z), and they are the ones the HER2 example
  itself uses. Using them keeps our epitope claim traceable to a structure
  rather than to a guess.
- **Single-patch, not multi-patch**: **OUR JUDGEMENT**, and the honest framing
  is that **there is no evidence either way.** "Multi-patch" is not a documented
  concept on this platform; multiple atoms across multiple chains are
  demonstrably supported (the HER2 spec spans chains C and D), but a second
  *distinct* epitope patch per run is neither documented nor forbidden. We go
  single-patch because our scored property is **competition with trastuzumab**,
  which is a claim about one epitope — splitting hotspots across two patches
  would dilute the only thing we are being measured on. **No one should cite
  "multi-patch beats single-patch" or its converse: we have no source for it.**
- `resampleHotspots: false` (CITED default). Leave it off: with only 4
  cited atoms, replacing them with 2-3 sampled neighbours would discard the
  structural anchoring that is the whole reason we trust them.

---

## 4. Sequence liability bans — thresholds, each labelled

**Framing.** RFD4-Proteina optimises structure; it has **no developability
objective**. So these are post-generation filters, run before refold GPU is
spent, since all are computable from sequence alone (Biopython is installed in
`/home/hackathon/nick-work/nmenv/`).

**A cache-wide scan found zero records stating a numeric developability cutoff
for any of these.** The *existence* of each liability class is well cited; the
*numbers* mostly are not. That split is marked per row.

### HARD KILL — discard the design, do not re-rank

| Liability | Threshold | Status | Source / reasoning |
|---|---|---|---|
| Free / unpaired cysteine | **0 Cys in CDRs**; only the 2 canonical framework Cys allowed | **CITED** | `abguard` rejects unpaired cysteines outright — `apheris/challenge/ai-antibodies-hackathon-apheris.md:236-246`. A hard submission gate, not our preference. |
| N-glycosylation sequon | **No `N-X-S/T` (X != P) anywhere in the three CDRs** | **CITED for the ban, OUR JUDGEMENT for the scope** | PMID:24492299 (mAbs 2014, 29 cites): affinity-improving mutations at heavy-chain position 53 "introduced glycosylation sites in heavy chain CDR 2 (CDR H2) that **impaired binding** of antibodies expressed in mammalian cells", and the one KD<10 pM variant kept was the combination that **avoided** CDR H2 glycosylation. That is a measured case of a CDR sequon destroying a good binder. Extending it from CDR H2 to all three CDRs is ours. The `N-X-S/T, X != P` motif definition itself is textbook, not from this cache. |
| Total length | **<= 250 aa** | **CITED** | `abguard`, same file. A VHH is ~125 aa, so this is slack for us. |
| Novelty | **< 70% sequence identity to any known antibody** | **CITED** | `abguard`, de novo track rule, same file. Note this is the **de novo** rule; the 70%-CDRH3-identity-to-trastuzumab rule applies to the *optimisation* track and does not bind us. |
| `immunum` numbering | **confidence >= 0.1** | **CITED** | `abguard`, same file. |

### SOFT RE-RANK — record, penalise, keep

| Liability | Threshold | Status | Source / reasoning |
|---|---|---|---|
| Deamidation / isomerisation motifs (`NG`, `NS`, `DG`) in CDRs | **<= 1 across all three CDRs** | **OUR JUDGEMENT — no citation** | The chemistry is textbook; **no cached record gives a per-design count**, and nothing in the cache bans these in CDRs. Soft, not hard, precisely because we cannot defend the number. Queued Q3 would fix this. |
| Arginine fraction | **<= 10% of CDR residues, <= 2 Arg in CDR3** | **OUR JUDGEMENT for the numbers; CITED for the direction** | PMID:28778924: scFv selected under weak specificity pressure were Arg-enriched and poorly specific, and alanine scanning showed their affinity depended heavily on those Arg, whereas genuinely specific variants did not. Mechanism and direction are cited; **the source publishes no cap**. Keep it soft and record Arg fraction for every design, so we can test whether Arg fraction correlates with min_ipae rank — if it does, our *metric* is rewarding non-specific binding (TH-003, IDEA-002). |
| Hydrophobic patch | **<= 3 consecutive residues from {F,W,Y,L,I,V,M} in a CDR** | **OUR JUDGEMENT — no citation, and our crudest number** | Deliberately a sequence-run proxy, because no structure-based patch calculator is installed: no PyRosetta, no FoldX, no CamSol, no propka. A real SAP/patch score would need tooling we do not have. Queued Q4. Do not present this as a developability assessment. |

### Not adopted

Rigidity, B-factor and loop-entropy filters: **deliberately not built** (TH-004,
CITED PMID:29545810 vs PMID:26132144 — the literature is split, so any such
ranking is indefensible). This saves effort rather than costing it.

---

## 5. Thermostability proxies for a VHH — and is TEMPRO actually obtainable?

### Evidence

| Source | What it says | Status |
|---|---|---|
| PMID:39154093 (Sci Rep 2024) | **TEMPRO**: nanobody Tm from protein embeddings, **MAE 4.03 C, RMSE 5.66 C** — an error bar on exactly the quantity Adaptyv measures | CITED |
| PMID:41564239 (ACS AMI 2026) | **NBsTem**: `NBsTem_Tm` trained on **514 experimental Tm** values, **Pearson 0.83** on external test; `NBsTem_Q` theoretical model **accuracy 0.84** from MD on **704** nanobody structures; proposes screening criterion **"Tm > 65 C and Qclass IV"**; applied to INDI, found **~12%** thermostable nanobodies | CITED |
| PMID:36683173 | PLM + Rosetta-feature CNN for scFv — trained on **TS50, a proxy, not Tm**; recovered **18 of 20** known thermostable positions and 5 identical mutations on a held-out mAb | CITED (and **not applicable** — we are a VHH) |

### Is TEMPRO obtainable? **Yes — verified by installing and running it, not by reading about it.**

- Repo: **`github.com/Jerome-Alvarez/TEMPRO`** (Jerome Alvarez is the first
  author of PMID:39154093). Not the `BiomedicalMachineLearning` URL the brief's
  premise implied — **that was wrong**.
- **43 MB**, 5 archives at repo root; it is a Jupyter-notebook distribution, not
  a package (no `setup.py`, no module).
- **Pretrained weights ship in the repo** (`user.zip` → `saved_ANNmodels_1500epoch/`:
  `ESM_650M.keras`, `ESM_3B.keras`, `ESM_15B.keras`). No separate download. The
  Tm head is tiny: 1280 → 1 MLP, **328,833 params**.
- Deps: Python **3.9** + TensorFlow **2.10.1** + Keras 2.10.0 + fair-esm 2.0.0.
  The Python 3.9 pin is a non-issue (`uv` fetches CPython 3.9.25; TF 2.10.1 has
  a cp39 wheel).
- **Ran end to end**: a real VHH sequence → ESM-2 650M embedding (1280-d, on the
  H100) → **predicted Tm 75.55 C**.
- **Two-environment split is required**: embeddings in `nmenv` (torch 2.15 /
  CUDA 13.2), Tm head in a TF 2.10 env. TF 2.10 pins numpy 1.23 and fights torch
  2.15 in one env.
- **Caveats to carry.** (1) The repo has **no LICENSE file** — all rights
  reserved by default. Fine for internal ranking; flag it before redistributing
  anything derived from it. (2) Only the `ESM_650M` head is usable; `ESM_15B`
  needs ~58 GB of embeddings the author never uploaded. (3) We can reproduce
  TEMPRO's **Tm** but **not** the `Qclass IV` half of PMID:41564239's criterion,
  which needs their NBsTem_Q model or web server.

### Decision

**Rank on TEMPRO (ESM-2 650M head). Gate at Tm > 65 C, treated as necessary but
not sufficient.**

- TEMPRO's **MAE 4.03 C** is CITED and is what makes a *ranking* defensible —
  a predictor without an error bar cannot rank.
- The **Tm > 65 C** gate is CITED (PMID:41564239), but we can only satisfy half
  of the published criterion, so it must be described as a partial gate. Saying
  "we applied the published criterion" would be an overclaim.
- Because MAE is 4.03 C, **do not treat a 2 C gap between two designs as a
  ranking**. Gap must exceed ~4 C to mean anything (OUR JUDGEMENT, read directly
  off the cited MAE).
- Free cross-checks available with what is installed (`nmenv`): ESM-2
  pseudo-likelihood, AbLang / AntiBERTy likelihoods as humanness proxies, and
  ANARCI germline identity + VHH hallmark residues at IMGT 42/49/50/52 (ANARCI
  functionally verified, HMMER bundled in `nmenv/bin/`).
- **The scFv predictor question does not arise.** Leo's TH-007/IDEA-003 was
  blocked because the answer was format-dependent; we submit a VHH, so we get
  the Tm-trained side. **IDEA-003 is closed.**

---

## 6. How many designs per arm

### Evidence

| Source | What it says | Status |
|---|---|---|
| `apheris/rfd4/16-refolding.md:40` | "Upstream reports binder pass rates of **7.7% to 10.6%**, so plan **100 to 200 designs per target**. On an A100, refolding takes about **25 s per design per model**, plus 20 s for the apo run." | CITED (file) |
| `apheris/rfd4/15-antibody-design.md:78-79` | "Plan 100 to 200 designs per target: upstream reports binder pass rates around **8 to 11%**. On an A100 a VHH takes about **24 s per design** and a VH/VL about 36 s." | CITED (file) |
| `apheris/rfd4/16-refolding.md:37-38` | "The cut-offs were calibrated on **AlphaFold3** and can shift by model, so compare pass rates **within one model**." | CITED (file) |
| `apheris/rfd4/sources.md` | "upstream" = pinned `RosettaCommons/RFD4-Proteina` | CITED (file) |
| `apheris/challenge/...apheris.md:190` and `apheris/rfd4/12-her2.md:84` | `her2-novel-2026-11-09`: 100 designs, 90% distinct, 4% known binders, **86% pass the challenge filter**, interface error 6.8 Å. **"HER2 novel was tested on one seed."** | CITED (file) |
| AMBC_6poBHWLhkEhP77QCbUh5yjX60jv (ECS Meeting Abstracts 2026) | "experimental success rates of only about **10-20%** for RFdiffusion and other de novo design methods" — **a meeting abstract, not peer-reviewed** | CITED (weak) |
| PMID:41193805 (Nature 2025) | Initial designs at **tens to hundreds of nanomolar KD**; OrthoRep maturation → **single-digit nM** | CITED |
| `apheris/rfd4/09-parameters.md`, `18-vm.md` | `numReplicates` 1-**200**; `refoldTopN` default = top **20%** rounded up | CITED (file) |
| Whole Apheris doc set | **No numeric min_ipae / iptm / binder_plddt cutoff documented anywhere.** For the HER2 example `pass` is permanently **empty** (`18-vm.md:169-172`). | CITED (negative) |

### Decision

**200 designs on the primary arm; 100 on each control arm.**

- **100-200 per arm**: CITED, the docs' own planning number, and 200 is
  `numReplicates`' ceiling.
- **200 rather than 100 on the primary arm**: **OUR JUDGEMENT.** Arithmetic:
  200 designs → refold top 20% = ~40 refolds → at 7.7-10.6% pass, **~15-21
  designs clear**, which fills a 4-design submission with margin. At 24 s per
  VHH design that is ~80 min of A100 for generation, plus ~25 s per refold per
  model. The asymmetry decides it: under-sampling risks an empty funnel with no
  second chance before the deadline; over-sampling costs 80 min of shared GPU.
- Measure our own rate on the first 100 and correct (IDEA-010). The cited rate
  is upstream's, on other targets, with AF3-calibrated cut-offs.

### Two premises in the brief that need correcting

1. **"The RFD4 docs state upstream binder pass rates of 8-11%."** Nearly right,
   and worth stating precisely: the literal string `8-11%` appears in **no
   file**. The docs say "around **8 to 11%**" (`15-antibody-design.md`) and,
   more precisely, "**7.7% to 10.6%**" (`16-refolding.md`). Both are attributed
   to **upstream**, neither is HER2-specific, neither is per-weight, and no
   Apheris-measured binder pass rate exists in local docs. Quote **7.7-10.6%
   (upstream)**.
2. **"The Apheris HER2 page reports 86% novelty pass for one weight."**
   Confirmed, with a caveat that changes how much weight it carries: 86% is
   `her2-novel-2026-11-09`'s share of **100 designs on a single seed** passing
   the **optimisation track's novelty filter**. It is **not** a binder pass
   rate and **not** a refold pass rate. It also **does not transfer to us**:
   that weight only redesigns trastuzumab's CDRH3 on the Fv, we are on the de
   novo VHH route with `antibody-2026-10-30`, and the de novo novelty rule is a
   different rule (<70% identity to any known antibody).

### Our pass rule

Since **no numeric min_ipae / iptm / binder_plddt threshold is documented
anywhere**, and the HER2-style spec leaves the automated `pass` column empty
(relative rather than absolute coordinates, so the scorer cannot measure
displacement), **any threshold we apply on those three is OUR JUDGEMENT and
must be labelled as such in the pitch.** Rank within one folding model only
(CITED warning). The only numeric binder gates in the docs are the Adaptyv
miniprotein campaign's custom ones (strict: pose RMSD <= 2.5 Å, interface PAE
<= 1.5 Å, binder pLDDT >= 80), explicitly disclaimed as "custom strict and
relaxed filters, not the module's automated pass rate" — borrow their *shape*,
not their authority.

---

## 7. Bugs in the upstream `research_pass.py`

Both claims in the brief were checked against
`LEO/tasks/experiment-ledger/research_pass.py`. **Both are real.** Not fixed in
Leo's folder, which is read-only for us.

### Bug 1 — trialcore records are scored on fields that do not exist. **REAL.**

`score_record` builds its match blob from `title` + `abstract`:

```python
blob = f"{rec.get('title','')} {rec.get('abstract','')}".lower()
...
quotes = [s for s in sentences(rec.get("abstract", "")) if NUMBER_RE.search(s)]
```

A trialcore record has neither key. Verified on
`cache/trialcore__trastuzumab-her2-breast-cancer-resistance.json` (20 records):

```
{'has_title': 0, 'has_abstract': 0, 'has_briefTitle': 20, 'has_briefSummary': 17}
```

The fields present are `briefTitle`, `briefSummary`, `officialTitle`,
`primaryOutcomeMeasures`, `nctId`, `amassId`. So `blob` is `" "` for every
trialcore record, `hits` is empty, `len(hits) < MIN_KEYWORD_HITS` (3), and the
record scores **0** and can never surface. It fails **silently** — no crash,
because `.get('title','')` returns `''` for an absent key rather than `None` —
which is exactly why it presented as the "20 cached records returned empty
titles" blocker in Leo's `STATE.md`. **Fix shape**: resolve title and abstract
per core, e.g. `rec.get('title') or rec.get('briefTitle') or rec.get('officialTitle')`
and `rec.get('abstract') or rec.get('briefSummary')`.

### Bug 2 — the cache key collides. **REAL, and live already.**

```python
def slug(text): ...; return s[:60]
def cache_path(core, query): return os.path.join(CACHE, f"{core}__{slug(query)}.json")
```

Two independent problems:

- **Truncation at 60 chars.** **12 of the 32** existing cache filenames are
  already at exactly 60 characters, so there is zero headroom. Demonstrated:
  `trastuzumab alanine scanning mutagenesis her2 binding hot spot residues` and
  `trastuzumab alanine scanning mutagenesis her2 binding hot spot energetics 2026`
  both slug to `trastuzumab-alanine-scanning-mutagenesis-her2-binding-hot-sp`
  and resolve to **one file**. The second query silently returns the first
  query's cached results.
- **The key omits request-determining parameters.** `cache_path` takes only
  `core` and `query`. `fetch` also sends **`limit`**, which is not in the key —
  so `limit=20` and `limit=200` already share one cache file. The brief framed
  this as "collides once filters are added"; it is **worse than that**, because
  it does not wait for filters.

**Fix shape**: key on a hash of the full request — `core`, `query`, `limit` and
any future filters — e.g. `f"{core}__{slug(query)[:40]}__{sha1(canonical_params)[:12]}.json"`,
keeping a readable prefix for humans and a collision-free suffix for
correctness.

### A third thing, not a bug

`research_pass.py` has a working cache-only mode: `--max-new 0` makes zero
network calls and `read_key` is only called `if allowed_new`, so it does not
even look for a key. That is how a cache-only pass is reproducible here with
`AMASS_API_KEY` unset. Worth knowing before anyone "fixes" that branch.
