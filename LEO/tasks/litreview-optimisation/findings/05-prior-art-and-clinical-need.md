# Prior art + clinical need (pitch: "user fit" 20%, "path to value" 20%)

From `cache/patentcore__*.json` (38 distinct publications) and `cache/trialcore__*.json` (20 trials).
**Abstract-level only.** PatentCore search returns `hasClaims` as a *boolean* — no claim text, no
sequence listings. `fetch_fulltext.py` requests `include=claims&include=description` per patent, so
claim text should be available once that fetch completes; **re-verify every collision call against it.**

> Not legal advice and not a freedom-to-operate opinion. An FTO determination needs claim-set review
> and a sequence search by a qualified patent professional.

## Strategic headline

**The HER2 target is extremely crowded, but our specific route is comparatively open — and the
obvious alternative route is not.**

- **Trastuzumab's composition-of-matter IP has lapsed.** Biosimilars are in trials
  (`NCT03755141`, Herzuma®). The scaffold is free to build on. *(Inference from biosimilar
  existence, not a stated legal fact.)*
- **The "find a different HER2 epitope" design-around is already claimed by at least four separate
  groups** — AbClon `US-2016053011-A1`, GC Cell `US-12116417-B2` (granted 2024), MAb Discovery
  `US-11091561-B2` and `US-2022204644-A1`. Each defines its invention as binding an epitope *other
  than* trastuzumab's.
- So **rewriting CDRH3 while preserving trastuzumab's epitope is the less crowded strategy**, which
  happens to be exactly what our track's rules force. This is a pitch asset, not just a constraint.
- **Caution:** `US-12116417-B2` (GC Cell, granted) expressly claims superior killing of
  trastuzumab-resistant / low-sensitivity cells. **If our story is "a better antibody for
  trastuzumab-resistant disease", that narrative is already occupied.** Position on the molecule
  and the interface, not on the resistant-patient claim alone.

## Patents to check before asserting novelty

| Publication | Assignee (as stated) | Why it matters |
|---|---|---|
| **US-2021388063-A1** | Fusion Antibodies plc | **Most on-point.** Method of building variant libraries from somatic-hypermutation motifs, demonstrated *on trastuzumab*, and "the variants produced are also claimed" — product *and* method. Still an application; scope unknown |
| **US-8840896-B2** | Genentech (Lowman/Gerstner/Carter) | "anti-HER2 antibody variants having substitutions at positions within the variable domains" — genus-level, reads onto our concept. Priority 2002-04-10, so likely lapsed ~2022 *(arithmetic, no legal-status field)* |
| **US-10167342-B2** | Shokri/Jeddi-Tehrani/Amiri ("hersintuzumab") | Granted, anti-HER2 **mutated on CDRs and framework**. Collision depends on unseen SEQ ID NO: 1/5 |
| **US-12227591-B2** | Zymeworks BC | Granted 2025. Claim shape to fear: "modifications in framework and/or CDRs vs **wild-type parent**, increasing affinity **≥2-fold**" — parent-referenced and *functional*, not sequence-locked. Parent is 2C4/pertuzumab (ECD2), so likely outside, but assume a trastuzumab-parented analogue exists outside this 38-record slice |
| **US-10584178-B2** | Hoffmann-La Roche | "novel HER2 antibody variants" — broad phrase, one-line abstract, genuinely unassessable |
| **US-2025277053-A1** | Genentech / Roche | Only record naming the **4D5 epitope** (trastuzumab's). Claims a *regimen*, not a composition — published 2025-09-04 off a 2008 priority. The regimen space around trastuzumab is still actively prosecuted |

**Crowding evidence:** 36/40 records carry CPC `C07K16/32` (anti-ErbB2); `C07K2317/565` (CDR subject
matter) on 16; `C07K2317/92` (affinity characterisation) on 22. Genentech/Roche hold 26% of
publications and the four most-cited records. Two queries overlapped on only 2 publications, so the
true family count is far larger than 38 — **we are seeing a thin slice.**

**Claim forms used in granted precedent** (what our CDR set would have to clear): "same CDRs as
[clone X]"; "CDR H1, H2 and H3" triplets; "parent + CDR modifications + ≥2-fold affinity gain".

## Clinical need (all record-stated)

- **Resistance to the ADC generation is now the live problem**, not resistance to naked trastuzumab.
  Dedicated mechanism studies: `NCT03829306` (T-DM1, KATIA/GEICAM) and `NCT06330012` (T-DXd, spatial
  transcriptomics, recruiting to 2027).
- **Patients past two HER2-directed regimens have no defined option** — `NCT03755141` states
  efficacy of re-use beyond progression "has not been tested yet".
- **"Treatment beyond progression" is a twenty-year unresolved argument** — `NCT00444587`,
  `EUCTR2006-002860-25-BG` (PANDORA), `EUCTR2007-005021-30-IT` (TOP) all randomise continuation vs
  discontinuation. Clinicians are *recycling* agents because nothing better exists
  (`JPRN-UMIN000014189`, lapatinib reinduction then trastuzumab retreatment).
- **High-risk residual disease after neoadjuvant therapy**: `EUCTR2020-003681-40-DK`, 1700 patients.
- **Immunogenicity is tracked as a liability of every current HER2 biologic** — ADA monitoring in
  `NCT04784715`, `EUCTR2020-003681-40-DK`, `NCT03755141`. *A de novo or heavily rewritten variable
  domain carries real regulatory risk here and should be designed against from the start.*
- **Resource asymmetry (inference, but quantitative):** resistance-biology trials enrolled
  **10, 11, 32, 48, 50** patients (`NCT01265927`, `NCT03135171`, `NCT03829306`, `NCT03947242`,
  `NCT06330012`); ADC races enrolled **1157 and 1700** (`NCT04784715`, `EUCTR2020-003681-40-DK`).
  The field answers resistance by changing the **payload, not the binder**.

**Path to value — the defensible positioning:** not replacing naked trastuzumab in first line (ADCs
are winning that), but supplying a better **targeting module**. Every ADC, CAR and bispecific in this
landscape (`US-10160812-B2`, `US-2021087292-A1`, `US-12116417-B2`, `US-2023212312-A1`) needs a
binder, and a binder engineered for the resistant setting is the input they lack.

## Gaps that matter

1. **No claim text or sequence listings in the search cache** — so sequence identity against claimed
   CDRs **cannot be computed**. `include=claims` fetch pending; this is the top gap to close.
2. **No legal status, term, expiry or abandonment data.** Every expiry statement above is arithmetic
   on `priorityDate`.
3. **Continuation chains invisible** — `parentPublicationNumbers`, `childPublicationNumbers`,
   `citedPatents`, `priorityClaimNumbers` are **empty in all 40 records**. Pending continuations are
   the usual source of surprise claims, and records off 2000/2008 priorities prove they exist here.
4. **Jurisdictional blind spot:** 37/40 US, **zero EP, zero JP**.
5. **No results data.** 5 trials have `hasResults: true` but the schema carries only outcome
   *measure names* — so **the pitch cannot yet quantify how badly trastuzumab fails** (no "X%
   progress within Y months").
6. **No molecular resistance mechanism is actually reported** in these records. They state mechanisms
   are being *investigated*. **Do not cite p95HER2, PIK3CA/PTEN, HER2 heterogeneity or antigen loss
   to these records** — nothing in this cache names them.
7. **No query run on**: trastuzumab biosimilars, ML/computational antibody design as claimed art,
   CDRH3 libraries as such, or HER2-low. We have no read on whether "ML-designed CDRH3" is itself
   claimed by a platform company — **worth one query before the pitch.**
8. Query noise: 7/20 records in the affinity-maturation file are not HER2 at all (CD47, c-Met,
   Ang-2). Effective HER2 prior-art set is ~31 of 38.
9. All 20 trials are **breast cancer only**, though patents repeatedly claim HER2+ gastric and colon.
