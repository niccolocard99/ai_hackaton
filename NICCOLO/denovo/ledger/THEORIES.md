# Theories — de novo VHH binders to HER2 domain IV (trastuzumab epitope)

Campaign: de novo VHH binders to HER2 domain IV at trastuzumab's epitope,
generated with RFD4-Proteina. This is the **de novo track**, not the
optimisation track, so the novelty rule that binds us is "under 70% sequence
identity to any known antibody", not the CDRH3 rule.

Same schema as [`../../../LEO/tasks/experiment-ledger/THEORIES.md`](../../../LEO/tasks/experiment-ledger/THEORIES.md),
validated by the adapted [`validate.py`](validate.py) in this folder.
Append-only ids. `confidence`: `high` | `medium` | `low` | `retired`.

**Id policy.** Theories imported from Leo's ledger **keep their original ids**
(TH-003 … TH-008) so a claim is traceable to one address across both ledgers.
Ids new to this campaign start at **TH-009**, which is the first id Leo's
ledger has not used. Every imported block carries `imported_from` and
`imported_date`.

**Every theory below states what it implies for either (a) a generation arm
parameter or (b) a funnel filter threshold.** That is the `implies` field. A
theory with no `implies` is not actionable and does not belong here.

**CITED vs OUR JUDGEMENT** is marked inline on every number. A number with no
PMID or file path behind it is labelled `OUR JUDGEMENT` and must stay labelled.

Seeded 2026-11-14 by [`research/pass-001.md`](research/pass-001.md), a
**cache-only** pass: `AMASS_API_KEY` was unset and credits unclaimed, so zero
network calls were made and zero searches were spent. Queries we wanted but
could not run are queued in that pass file.

---

### TH-003 — Arginine-rich designed CDRs score well on predicted affinity but fail epitope specificity, because Arg-driven affinity is non-specific

- **confidence:** medium
- **status:** active
- **property:** binding; epitope
- **imported_from:** LEO/tasks/experiment-ledger TH-003
- **imported_date:** 2026-11-14
- **falsified_by:** High-Arg designs compete with trastuzumab for HER2 as well as low-Arg designs of matched min_ipae.
- **evidence_for:** PMID:28778924 (CITED) — scFv selected under weak specificity selection were enriched in Arg CDR mutations and were poorly specific; alanine scanning showed their affinity depended heavily on those Arg residues, whereas genuinely specific variants did not. Title: "Arginine mutations in antibody complementarity-determining regions display context-dependent affinity/specificity trade-offs"
- **evidence_against:**
- **implies:** **(b) funnel filter.** Cap Arg fraction over the three designed CDRs at **<= 10% of CDR residues, and no more than 2 Arg in CDR3** (OUR JUDGEMENT — PMID:28778924 establishes the *direction* and the mechanism but publishes no numeric cap; no cached record gives one). Applied as a *soft* re-rank, not a hard kill, with the count recorded per design so the correlation between Arg fraction and min_ipae can be measured. If high-Arg designs dominate the top of the min_ipae ranking, the ranking is rewarding non-specific binding and the metric, not the designs, is the problem.
- **query_terms:** arginine CDR mutation polyspecificity nonspecific binding antibody; charged residue CDR affinity specificity trade-off; antibody net charge nonspecificity developability
- **history:**
  - 2026-11-14 imported at Leo's confidence `medium`, id preserved (pass-001). Mechanism was shown on an anti-amyloid-beta scFv, not HER2 and not a VHH, so transfer is still an assumption. Numeric cap added here as OUR JUDGEMENT, since the source has none.

### TH-004 — Selecting designs on predicted loop rigidity or B-factor is not a useful filter, so we should not build one

- **confidence:** medium
- **status:** active
- **property:** binding
- **imported_from:** LEO/tasks/experiment-ledger TH-004
- **imported_date:** 2026-11-14
- **falsified_by:** Rigidity or B-factor ranking separates our strong from weak binders in the returned assay data.
- **evidence_for:** PMID:29545810 (CITED) — "Repertoire Analysis of Antibody CDR-H3 Loops Suggests Affinity Maturation Does Not Typically Result in Rigidification": rigidity theory over thousands of repertoire models found no clear delineation between naive and antigen-experienced antibodies, and only a slight CDR-H3 flexibility decrease across PDB structures
- **evidence_against:** PMID:26132144 (CITED) — QSFR analysis found rigidity emerges during affinity maturation in three antibody systems; PMID:21287614 argues the same. Smaller sample than PMID:29545810.
- **implies:** **(b) funnel filter — a negative.** Build **no** rigidity, B-factor, RMSF or loop-entropy filter. Zero engineering effort, zero GPU. This theory's value is the effort it saves: the literature is genuinely split, so any rigidity ranking we built would be indefensible in the pitch. **Note the scope limit we are importing:** this says nothing about pLDDT or PAE, which are *confidence* metrics, not rigidity metrics. We still rank on min_ipae / iptm / binder_plddt.
- **query_terms:** CDRH3 loop flexibility rigidity affinity maturation repertoire; antibody paratope conformational entropy binding; B-factor rigidity antibody affinity prediction
- **history:**
  - 2026-11-14 imported at Leo's confidence `medium`, id preserved (pass-001). Kept because a settled negative is cheap to carry and expensive to rediscover.

### TH-005 — Binding and stability can be decoupled by localisation: CDRs carry the binding risk, framework changes buy back stability

- **confidence:** medium
- **status:** active
- **property:** thermostability; epitope
- **imported_from:** LEO/tasks/experiment-ledger TH-005
- **imported_date:** 2026-11-14
- **falsified_by:** Framework-stabilised variants of our designs lose HER2 competition, or gain no predicted Tm.
- **evidence_for:** PMID:40653772 (CITED) — on trastuzumab itself, computationally designed framework mutations **distal to the CDRs** were stabilising and preserved HER2 engagement and ADCC; PMID:32336757 (CITED) — computational engineering of thermostabilised antibody fragments; PMID:22272277 (CITED) — engineered intradomain disulfides stabilised IgG1 Fc
- **evidence_against:** PMID:40653772 (CITED) — the same work found a light-chain framework mutation that kept antigen binding but lost effector function, so framework changes are not free of non-local consequence (harmless for our three scored properties, which exclude effector function)
- **implies:** **(a) generation arm — but CONSTRAINED, and this is where the import is weakest.** For this campaign the framework is **fixed, not designable**: every documented RFD4-Proteina antibody route keeps framework sequence and internal geometry and designs only CDR loops (`/home/hackathon/hackathon-team-3/apheris/rfd4/02-miniprotein.md:30-31`, `15-antibody-design.md:22-24`; CITED to those files). Framework redesign is **undocumented** on this platform. So the actionable form here is **framework *choice*, not framework design**: pick a VHH framework already known to be stable and keep it. The supplied example framework is caplacizumab (a marketed, humanised VHH), extracted from `/home/hackathon/hackathon-team-3/apheris/rfd4/Data/03-vhh/pmvk_caplacizumab.cif.zst`: `EVQLVESGGGLVQPGGSLRLSCXWFRQAPGKGRELVAXYPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCXWGQGTQVTVS` (86 positions, the three `X` at label_seq_id 23/38/76 are the CDR sentinels). **Use it as-is** (OUR JUDGEMENT: a clinically validated framework is the cheapest stability we can buy, and it costs nothing to adopt). Any post-hoc stabilising mutation must not leave an unpaired cysteine, which `abguard` rejects outright.
- **query_terms:** framework mutation antibody thermostability preserve antigen binding trastuzumab; VHH framework stabilizing mutation melting temperature; consensus framework engineering single domain stability
- **history:**
  - 2026-11-14 imported at Leo's confidence `medium`, id preserved (pass-001). **Re-scoped on import.** Leo's version assumed we could design framework residues; the RFD4-Proteina docs show we cannot, so the theory is retained in its "choose a stable framework and keep it fixed" form. Not lowered in confidence, because the underlying claim (binding risk and stability risk are spatially separable) is unaffected — only our lever on it is.

### TH-006 — A VHH can hold trastuzumab's own epitope at sub-nanomolar affinity, so the de novo VHH route is viable without giving up epitope preservation

- **confidence:** medium
- **status:** active
- **property:** binding; epitope
- **imported_from:** LEO/tasks/experiment-ledger TH-006
- **imported_date:** 2026-11-14
- **falsified_by:** Designed VHHs that rank well on min_ipae fail competitive binding against trastuzumab, indicating the epitope needs a paired VH/VL interface.
- **evidence_for:** PMID:37632647 (CITED, 2023, Immunologic Research, 9 cites) — "A high affinity and specificity anti-HER2 single-domain antibody (VHH) that targets trastuzumab's epitope". A llama-derived, phage-displayed anti-HER2 VHH with **sub-nanomolar** affinity measured by **biolayer interferometry**, which **"binds to the same epitope as trastuzumab"**, confirmed by epitope mapping plus in-silico modelling and docking.
- **evidence_against:**
- **implies:** **This is the campaign's feasibility premise — it is why the arm exists at all, and it sets the target, not a filter.** (a) generation arm: the whole VHH arm is justified by this one existence proof. Target = **sub-nanomolar KD**, i.e. the bar a *discovered* VHH already cleared at this exact epitope (CITED). (b) funnel: **we must run the competition check, not assume it.** The paper's VHH was *selected* by phage display against a real antigen; ours are *generated*, so epitope identity is a design intent, not a measured fact, until the competition assay returns. Honest caveats to carry into the pitch: n=1 existence proof, affinity by BLI not SPR, "same epitope" by mapping plus docking rather than a co-crystal, and the VHH sequence itself is in the paper, so a design too close to it is a novelty risk (see IDEA-007).
- **query_terms:** VHH nanobody HER2 trastuzumab epitope competitive binding; single domain antibody domain IV HER2 affinity; nanobody epitope mapping overlap monoclonal
- **history:**
  - 2026-11-14 imported at Leo's confidence `medium`, id preserved (pass-001). Kept at `medium`, not raised: a single existence proof justifies the route but predicts nothing about our hit rate. Not lowered either — it is primary, measured, and at exactly our epitope, which is the strongest kind of prior we have.

### TH-008 — The 1N8Z contact set tells us where to put hotspots, but a contact is not a hot spot, so it cannot tell us which designs preserved the epitope

- **confidence:** low
- **status:** active
- **property:** epitope
- **imported_from:** LEO/tasks/experiment-ledger TH-008
- **imported_date:** 2026-11-14
- **falsified_by:** Designs that conserve the 1N8Z contact set still lose competitive binding, or designs that break it still win it.
- **evidence_for:** PMID:12610629 (CITED) — structure of HER2 with the Herceptin Fab, the source of the domain IV contact set; PMID:31042744 (CITED) — cryo-EM of the HER2–trastuzumab–pertuzumab complex
- **evidence_against:**
- **implies:** **(a) generation arm, and an explicit refusal of (b).** The contact set is **good enough to place hotspots** and that is already how the supplied HER2 spec does it: `/home/hackathon/hackathon-team-3/apheris/rfd4/Data/12-her2/her2_cdrh3.yaml` marks **4 HER2 atoms across 2 chains** as `C_HOT` — `sel('D/*/593/NZ')`, `sel('C/*/573/CD1')`, `sel('D/*/591/CG2')`, `sel('C/*/572/CD')` — with the comment "C_HOT marks the HER2 atoms nearest CDRH3 in PDB 1N8Z as contact hotspots" (CITED to that file). Reuse those 4 atoms for the VHH arm: they are the only epitope definition in the repo that is anchored to a real co-crystal. **Build no "contacts preserved" scoring filter** — that would convert an unweighted contact list into a pretended energy, which is exactly the error this theory names. Epitope preservation is scored by competition, by the assay, not by us.
- **query_terms:** trastuzumab HER2 contact residues 1N8Z interface paratope; alanine scanning trastuzumab HER2 hot spot energetic contribution; buried surface area shape complementarity HER2 antibody interface
- **history:**
  - 2026-11-14 imported at Leo's confidence `low`, id preserved (pass-001). Stays `low`: no cached record gives a per-residue energetic contribution at this interface, so "load-bearing" is still unmeasured. The import is still useful because it fixes our hotspot atoms to a cited source instead of a guess.

### TH-009 — A de novo VHH CDR3 must be long enough to reach the epitope but short enough that our own ranking oracle is still accurate, and that upper bound is tighter than the platform default

- **confidence:** medium
- **status:** active
- **property:** binding
- **new_in:** this campaign
- **falsified_by:** Designs with CDR3 of 17-21 residues rank and validate as well as 12-16 residue designs, i.e. the predicted-structure error at long H3 does not propagate into min_ipae reliability. Also falsified if nothing in the 12-16 window reaches the epitope at all.
- **evidence_for:** PMID:19062174 (CITED, Proteins 2009, 188 cites) — "Toward high-resolution homology modeling of antibody Fv regions": median CDR-H3 backbone RMSD of prediction degrades monotonically and then collapses with loop length — **1.6 Å (4-6 residues), 1.9 Å (7-9), 2.4 Å (10-11), 3.1 Å (12-14), 6.0 Å (17-22)**. PMID:23255066 (CITED) — ab initio H3 prediction benchmarked over 53 loops of length 4-22. PMID:41193805 (CITED, Nature 2025, 122 cites) — de novo RFdiffusion VHHs bound 4 disease-relevant epitopes with designed CDRs confirmed atomically accurate by cryo-EM, i.e. de novo CDR design does reach real epitopes. Platform default range for a VHH, read out of the supplied annotated CIF `/home/hackathon/hackathon-team-3/apheris/rfd4/Data/03-vhh/pmvk_caplacizumab.cif.zst`: **CDR1 13-13, CDR2 9-10, CDR3 12-21** (CITED to that file).
- **evidence_against:** PMID:8918938 (CITED) — the 1230-fold-improved anti-c-erbB-2 scFv had a **20 aa** VH CDR3, and PMID:25484051 (CITED) — a crystal structure showed a **24 aa** protruding VH CDR3 "critical to the predicted binding mode". Long CDR3s demonstrably work in nature and in libraries; the constraint we are imposing is about *our ability to rank them*, not about whether they can bind. PMID:18687338 (CITED) — a soluble, stable human VH domain library spanned CDR3 lengths **7 to 24 residues**, all purifiable from the E. coli periplasm, so length alone does not predict expressibility.
- **implies:** **(a) generation arm parameter, directly.** Set the designed CDR3 length range to **12-16** instead of the platform's 12-21, keeping CDR1 at 13 and CDR2 at 9-10 as supplied. The 12 floor is the platform's own (CITED); the 16 ceiling is **OUR JUDGEMENT**, reasoned as: our entire funnel ranks on predicted-structure quantities (min_ipae, iptm, binder_plddt), PMID:19062174 shows H3 prediction error roughly doubles from 3.1 Å at 12-14 to 6.0 Å at 17-22, and a 6 Å backbone error makes an interface-PAE ranking meaningless. We would rather design in the window where our ranking is trustworthy than design long loops we cannot triage. Carry **one deliberate long-CDR3 arm (17-21)** as a control so this bound is tested rather than assumed.
- **query_terms:** CDRH3 loop length conformation canonical structure antibody; nanobody CDR3 length distribution camelid repertoire; H3 loop modelling accuracy RMSD length dependence
- **history:**
  - 2026-11-14 created, confidence medium (pass-001). Medium not high because the key inference — that model-era H3 prediction error at long loops propagates into *our* min_ipae ranking — is an extrapolation from a 2009 homology-modelling benchmark to a 2026 co-folding model, and is untested. The long-CDR3 control arm exists to test exactly that.

### TH-010 — Sequence liabilities must be filtered post-generation, because RFD4-Proteina optimises structure and has no developability objective at all

- **confidence:** medium
- **status:** active
- **property:** thermostability; binding
- **new_in:** this campaign
- **falsified_by:** Our generated CDR sets contain no liability motifs at any meaningful rate, making the filter a no-op; or liability-carrying designs score no worse on the returned assays.
- **evidence_for:** PMID:24492299 (CITED, mAbs 2014, 29 cites) — "In vitro affinity maturation of a natural human antibody overcomes a barrier to in vivo affinity maturation": affinity-improving mutations at heavy-chain position 53 **introduced N-glycosylation sites in CDR H2 that impaired binding of antibodies expressed in mammalian cells**, and the one high-affinity variant kept (KD < 10 pM) was specifically the combination that **avoided** CDR H2 glycosylation. This is a measured case of a CDR N-glycosylation motif destroying an otherwise-good binder. `abguard` hard rules, CITED to `/home/hackathon/hackathon-team-3/apheris/challenge/ai-antibodies-hackathon-apheris.md:236-246`: no unpaired cysteines; max total length 250 aa; `immunum` numbering confidence >= 0.1; under 70% sequence identity to any known antibody for the de novo track.
- **evidence_against:** No cached record gives a numeric threshold for Arg fraction, hydrophobic-patch size, or a CDR deamidation-motif ban. A cache-wide scan found 0 records stating any such cutoff. The thresholds below are therefore **unsupported by our literature** and labelled accordingly — the *existence* of the liability classes is cited, their *cutoffs* are not.
- **implies:** **(b) funnel filter thresholds, split by whether they are hard or soft.** HARD KILL (a design carrying one is discarded, not re-ranked): unpaired/free cysteine — count Cys in the full VHH, require the 2 framework Cys only, 0 in CDRs (CITED, `abguard` rejects outright); N-glycosylation sequon **N-X-S/T where X is not P** anywhere in the three CDRs (motif definition is textbook; the *ban* is justified by PMID:24492299, CITED; extending it from CDR H2 to all CDRs is OUR JUDGEMENT); total length <= 250 aa (CITED, `abguard`); <= 70% identity to any known antibody (CITED, `abguard`, de novo track). SOFT RE-RANK (recorded, penalised, not killed): `NG`, `NS`, `DG` motifs in CDRs — allow at most **1 across all three CDRs** (OUR JUDGEMENT, no citation; the chemistry is textbook but no cached source gives a per-design count); Arg fraction per TH-003 (OUR JUDGEMENT); hydrophobic patch — no more than **3 consecutive residues from {F,W,Y,L,I,V,M}** in a CDR (OUR JUDGEMENT, no citation, deliberately crude because we have no structure-based patch calculator installed — no PyRosetta, no FoldX, no CamSol). Every one of these is computable from sequence alone with Biopython, which is installed.
- **query_terms:** antibody developability sequence liability deamidation isomerization glycosylation CDR; nanobody aggregation hydrophobic patch prediction; N-glycosylation sequon CDR affinity loss expression
- **history:**
  - 2026-11-14 created, confidence medium (pass-001). The *classes* are well cited; the *numbers* mostly are not, and the `implies` field says which is which line by line. A pass that finds a real developability-threshold source should raise this and replace the judgement calls.

### TH-011 — Per-arm design counts must be set from the ~8-11% upstream binder pass rate, not from GPU convenience, or the funnel starves before the assay

- **confidence:** medium
- **status:** active
- **property:** binding
- **new_in:** this campaign
- **falsified_by:** Our own refold pass rate on a first 100-design batch is far from 8-11% — materially higher (the budget is wasteful) or near zero (the arm is broken and more designs will not fix it).
- **evidence_for:** CITED to Apheris docs, with the exact wording, because the rounded form is what our brief carried and the precise form is what the docs say: `/home/hackathon/hackathon-team-3/apheris/rfd4/16-refolding.md:40` — "Upstream reports binder pass rates of **7.7% to 10.6%**, so plan 100 to 200 designs per target."; `/home/hackathon/hackathon-team-3/apheris/rfd4/15-antibody-design.md:78-79` — "Plan 100 to 200 designs per target: upstream reports binder pass rates around **8 to 11%**. On an A100 a VHH takes about **24 s per design**." Independent and weaker: AMBC_6poBHWLhkEhP77QCbUh5yjX60jv (ECS Meeting Abstracts 2026, **a meeting abstract, not peer-reviewed**) — "experimental success rates of only about **10-20%** for RFdiffusion and other de novo design methods".
- **evidence_against:** Both Apheris figures are attributed to **upstream** (pinned RosettaCommons/RFD4-Proteina), not measured by Apheris, and neither is HER2-specific or per-weight. `16-refolding.md:37-38` warns "The cut-offs were calibrated on AlphaFold3 and can shift by model, so compare pass rates within one model." For the HER2 example the automated `pass` column is documented as permanently **empty** (`18-vm.md:169-172`), so no pass rate is even computable on that route.
- **implies:** **(a) generation arm parameter, directly.** **100-200 designs per arm** (CITED — the docs' own planning number, and `numReplicates` maxes at 200 per input). At 24 s per VHH design on an A100 (CITED), 200 designs is about **80 min of GPU per arm**, which is what makes 3-4 arms affordable in one afternoon. Refold the top 20% (`refoldTopN`, the documented default) — so 200 designs yields ~40 refolds, and at a 7.7-10.6% pass rate ~15-21 designs clear, which is enough to fill a 4-design submission with margin. **OUR JUDGEMENT:** run **200, not 100,** on the primary arm and 100 on each control arm, because the only quantity we actually need is 4 submittable designs and the cost of under-sampling (an empty funnel, no second chance before the deadline) is far worse than the cost of over-sampling (80 min of shared GPU).
- **query_terms:** de novo binder design experimental success rate hit rate RFdiffusion; nanobody de novo design yeast display screening hit rate; computational binder design pass rate refolding filter
- **history:**
  - 2026-11-14 created, confidence medium (pass-001). Medium because the number is real but borrowed: it is upstream's, on other targets, and our brief's "8-11%" turned out to be the docs' rounding of 7.7-10.6%. The first 100-design batch measures our own rate and should move this.
