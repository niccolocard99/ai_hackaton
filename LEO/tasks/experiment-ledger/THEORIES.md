# Theories — what we currently believe

Each theory is a falsifiable claim about how to beat trastuzumab on the three
measured properties. Confidence moves only when an experiment verdict or a
literature pass brings something new, and every move appends a dated line to
that theory's `history` citing the `exp_id` or pass id that caused it. Claims
are never silently edited — supersede them instead.

`confidence`: `high` | `medium` | `low` | `retired`.
`query_terms` drives [`research_pass.py`](research_pass.py); prefix a term with
`patent:` to send it to patentcore instead of biomedcore.

Seeded 2026-11-14 from the 30 cached AMASS queries in
[`../litreview-optimisation/`](../litreview-optimisation/) (pass-001). That
review was still running when these were written, so every confidence here is
provisional and several theories are deliberately stated more strongly than the
evidence yet supports, so that a pass can knock them down.

---

### TH-001 — Diversifying CDRH3 wholesale raises anti-HER2 affinity rather than destroying it, so our <70% CDRH3 identity constraint is an opportunity, not a tax

- **confidence:** high
- **status:** active
- **property:** binding
- **falsified_by:** A design series whose CDRH3 is >30% substituted shows systematically worse SPR KD than trastuzumab, with no high-affinity tail.
- **evidence_for:** PMID:8918938 — sequential VL-CDR3 then VH-CDR3 library mutation of an **anti-c-erbB-2 (HER2)** scFv gave 1230-fold KD improvement to 1.3e-11 M, with a 20 aa VH CDR3, i.e. long enough to tolerate library-scale change; PMID:29702641 — RosettaAntibodyDesign replaced individual CDRs with new **lengths and clusters** and improved affinity **10-50 fold, verified experimentally** on both a lambda and a kappa antibody-antigen complex; PMID:26225765 — targeted CDR3 randomisation of an anti-HER2 antibody optimised therapeutic targeting; PMID:35545451
- **evidence_against:** EXP-002 — wholesale CDRH3 change is no longer epitope-free: 4 substitutions are needed for the identity gate and only 3 CDRH3 positions never contact HER2, so an affinity-motivated redesign necessarily touches a contact (TH-009). The affinity case is strong; the *freedom* is not
- **query_terms:** CDRH3 library randomization affinity gain anti-HER2 scFv; CDR3 diversification picomolar affinity maturation; long CDRH3 20 residue loop tolerance substitution
- **history:**
  - 2026-11-14 created, confidence medium — one strong primary result (8918938) but on C6.5, not trastuzumab, and KD measured by a 1996 method (pass-001)
  - 2026-11-14 medium → **high** (pass-003). PMID:29702641 is a second, independent, experimentally verified result: 10-50 fold affinity gains from replacing whole CDRs including length changes. Combined with PMID:8918938, which is on an anti-HER2 scFv, two primary results now support wholesale CDR change as an affinity route. Note what did *not* move: EXP-002 added counter-evidence that the route is not epitope-free, so this theory is now "CDRH3 redesign can win affinity", not "CDRH3 redesign is safe".

### TH-002 — Affinity gains bought in CDRH3 cost thermostability, so the two prize properties are in direct tension and must be co-optimised, not optimised in sequence

- **confidence:** low
- **status:** active
- **property:** binding; thermostability
- **falsified_by:** Our top-affinity designs show no Tm penalty relative to trastuzumab, or the penalty is uncorrelated with the number of CDRH3 substitutions.
- **evidence_for:** PMID:30666176 — review of affinity/specificity/stability/solubility trade-offs, states improvements in one property commonly cause deficits in another; PMID:21540796
- **evidence_against:** PMID:40653772 — stabilising framework mutations distal to the CDRs preserved HER2 engagement, which suggests the tension can be *sidestepped* rather than traded (see TH-005); **EXP-005 — CDRH3 is thermostability-inert** (41 variants of the identical VH differing only in CDRH3 span 4.63 °C with 0 of 40 beyond the 4.1 °C MAE; Spearman(predicted Tm, measured KD) = +0.001), so an affinity gain bought in CDRH3 cannot cost Tm in the first place; AMBC_JwGSFsFYuumYi1zgXk0571g2F5A — an AI platform reports **Tm up to ~10 °C higher** with aggregation under accelerated thermal stress falling from 100% to <5% *while* optimising affinity; PMID:23995618 — CDR-restricted engineering of native human scFvs produced highly stable, soluble bifunctional antibodies formulated at 100 mg/ml with <2% high-molecular-weight species after 7 weeks
- **query_terms:** antibody affinity stability trade-off CDR mutation destabilizing; CDR mutation melting temperature penalty scFv; co-optimization affinity stability antibody engineering
- **history:**
  - 2026-11-14 created, confidence medium — the supporting source is a review, not primary measurement; the counter-evidence is primary and trastuzumab-specific (pass-001)
  - 2026-11-14 medium → **low** (EXP-005, pass-006). The mechanism this theory asserts does not exist in our design space: EXP-005 shows CDRH3 changes do not move Tm at all, with a control proving the predictor is not blind. Two further cases show affinity and stability being improved together rather than traded (AMBC_JwGSFsFYuumYi1zgXk0571g2F5A, PMID:23995618). **Downgraded, not retired**, and the distinction matters: the Tm figures in EXP-005 are *predicted* by a model that is out-of-domain for a VH (TH-016), and we have no *measured* Tm for any CDRH3 variant. The field-level trade-off (PMID:30666176) may still hold for framework and whole-molecule changes — it just does not apply to the CDRH3-only redesign our tools actually perform.

### TH-003 — Arginine-rich designed CDRH3 loops will score well on predicted affinity but fail epitope preservation, because Arg-driven affinity is non-specific

- **confidence:** medium
- **status:** active
- **property:** binding; epitope
- **falsified_by:** High-Arg CDRH3 designs compete with trastuzumab for HER2 as well as low-Arg designs of matched predicted affinity.
- **evidence_for:** PMID:28778924 — scFv selected under weak specificity selection were enriched in Arg CDR mutations and were poorly specific; affinity of the *specific* variants depended much less on their Arg residues, i.e. over-reliance on Arg for affinity reduces specificity
- **evidence_against:**
- **query_terms:** arginine CDR mutation polyspecificity nonspecific binding antibody; charged residue CDR affinity specificity trade-off; antibody net charge nonspecificity developability
- **history:**
  - 2026-11-14 created, confidence medium — mechanism shown on an anti-amyloid-beta scFv, not HER2; transfer is an assumption (pass-001)

### TH-004 — Selecting designs on predicted CDRH3 loop rigidity is not a useful filter, so we should not spend effort on it

- **confidence:** medium
- **status:** active
- **property:** binding
- **falsified_by:** Rigidity or B-factor based ranking separates our strong from weak binders in the returned assay data.
- **evidence_for:** PMID:29545810 — rigidity theory over thousands of repertoire models found no clear delineation between naive and antigen-experienced antibodies, and only a slight CDR-H3 flexibility decrease in PDB structures; PMID:26132144 and PMID:21287614 argue the opposite and are the minority position by sample size
- **evidence_against:** PMID:26132144 — QSFR analysis found rigidity emerges during affinity maturation in three antibody systems
- **query_terms:** CDRH3 loop flexibility rigidity affinity maturation repertoire; antibody paratope conformational entropy binding; B-factor rigidity antibody affinity prediction
- **history:**
  - 2026-11-14 created, confidence medium — this is a useful *negative* that saves effort if it holds; the literature is genuinely split (pass-001)

### TH-005 — Stability and binding can be decoupled: put all binding risk in CDRH3 and buy back thermostability with framework mutations that do not touch the epitope

- **confidence:** high
- **status:** active
- **property:** thermostability; epitope
- **falsified_by:** Framework-stabilised variants of our designs lose HER2 competition against trastuzumab, or gain no Tm.
- **evidence_for:** PMID:40653772 — on trastuzumab specifically, computationally designed framework mutations distal to the CDRs were stabilising and preserved HER2 engagement and ADCC; PMID:32336757 — computer-based engineering of thermostabilised antibody fragments; PMID:22272277 — engineered intradomain disulfides stabilised IgG1 Fc
- **evidence_against:** PMID:40653772 — the same work found one light-chain framework mutation that kept antigen binding but lost effector function, so framework changes are not free of functional consequence (harmless for our three properties, which exclude effector function)
- **query_terms:** framework mutation antibody thermostability preserve antigen binding trastuzumab; scFv framework stabilizing mutation melting temperature; consensus framework engineering Fv stability
- **history:**
  - 2026-11-14 created, confidence medium — strongest theory in this set because the evidence is primary, computational-plus-experimental, and on trastuzumab itself (pass-001)
  - 2026-11-14 medium → **high** (EXP-005). The decoupling now rests on a control, not one paper. 41 variants of the identical trastuzumab VH differing *only* in CDRH3 span 4.63 °C of predicted Tm with **0 of 40 beyond the predictor’s 4.1 °C MAE**, and Spearman(predicted Tm, measured KD) = **+0.001** — CDRH3 is inert. Framework point mutations in the same VH do move Tm (mean −0.76 °C, range −6.96 to +2.05, 3 beyond MAE) against CDRH3’s −0.09 °C and none beyond; |ΔTm| median 0.63 vs 0.19 °C, Mann–Whitney **p = 0.002**. So the model is not blind and the two regions genuinely differ. Design consequence: spend CDRH3 entirely on affinity and epitope. Tempered by TH-017 — the framework gains actually available are small and contaminated.

### TH-006 — A VHH can hold trastuzumab's own epitope at sub-nanomolar affinity, so the de novo track is viable without giving up epitope preservation

- **confidence:** medium
- **status:** active
- **property:** binding; epitope
- **falsified_by:** Designed VHHs that score well on min_ipae fail competitive binding against trastuzumab, indicating the epitope needs a paired VH/VL interface.
- **evidence_for:** PMID:37632647 — an anti-HER2 VHH with sub-nanomolar affinity binds the same epitope as trastuzumab, confirmed by epitope mapping and docking
- **evidence_against:**
- **query_terms:** VHH nanobody HER2 trastuzumab epitope competitive binding; single domain antibody domain IV HER2 affinity; nanobody epitope mapping overlap monoclonal
- **history:**
  - 2026-11-14 created, confidence medium — a single existence proof, which is enough to justify the route but not to predict our hit rate (pass-001)

### TH-007 — scFv thermostability is predictable from sequence well enough to rank our designs, which closes the open "no thermostability predictor" decision

- **confidence:** retired
- **status:** retired
- **property:** thermostability
- **falsified_by:** The chosen predictor's ranking shows no correlation with the measured Tm values Adaptyv returns.
- **evidence_for:** PMID:39154093 — TEMPRO predicts single-domain antibody Tm from protein embeddings with **MAE 4.03 °C, RMSE 5.66 °C**, i.e. a stated error bar on the quantity Adaptyv actually measures; PMID:41564239 (2026) — synergistic deep-learning nanobody thermostability prediction proposing a concrete screening criterion, **Tm > 65 °C and Qclass IV**; PMID:36683173 — PTLM and Rosetta-feature CNN classify thermostable scFv variants, and on an independent mAb with 20 characterised thermostable mutations the TS50-trained models recovered **18 of the residue positions and 5 identical amino-acid mutations**
- **evidence_against:** PMID:36683173 — trained on TS50, a proxy, not Tm, and out-of-distribution generalisation is explicitly qualified; the two strongest predictors (PMID:39154093, PMID:41564239) are **nanobody/sdAb-specific**, so they only apply to a VHH submission, not to an scFv
- **query_terms:** scFv thermostability prediction machine learning validation Tm; antibody melting temperature predictor benchmark accuracy; TEMPRO nanobody Tm protein embedding prediction
- **history:**
  - 2026-11-14 created, confidence low — this is the team's named open decision and the candidates are all trained on proxies; needs a pass of its own (pass-001)
  - 2026-11-14 low → **medium** (pass-002). Two things changed. TEMPRO now gives us a numeric error bar on predicted Tm (MAE 4.03 °C), which is what a ranking needs; and PMID:36683173's 18-of-20 position recovery is stronger transfer evidence than the TS50 caveat alone suggested. Not `high`, because the predictor choice turns out to be **format-dependent**: TEMPRO and PMID:41564239 are nanobody-only, so an scFv submission still rests on TS50-trained models. This splits IDEA-003 into a VHH answer and an scFv answer.
  - 2026-11-14 **retired** (EXP-005). The claim was that *scFv* thermostability is predictable well enough to rank, closing the open decision. The decision is now closed — with the opposite answer for scFv. Neither platform ships a Tm predictor (Nebius `stab-ddg` is StaB-ddG, *binding* ddG; Foundry’s AQ workflows are small-molecule tools), cheap features reach only Pearson +0.40 / MAE 7.6 °C against 8.3 °C for guessing the mean, and the one usable tool — NanoMelt, MAE 4.1 °C — is trained on 639 **nanobody** measurements and shows domain artefacts on a VH (11 of its top 12 single mutations are →Glu; `V37F` is a VHH hallmark position). Retired rather than downgraded, because the subject of the claim is exactly the case it cannot serve. Superseded by TH-016.

### TH-008 — The 1N8Z contact residues define which CDRH3 positions are load-bearing for the epitope, so a design that preserves those contacts preserves the epitope

- **confidence:** retired
- **status:** retired
- **property:** epitope
- **falsified_by:** Designs that conserve the 1N8Z contact set still lose competitive binding, or designs that break it still win it. — **Falsified before any design, by structure: see EXP-002.**
- **evidence_for:** PMID:12610629 — structure of HER2 with the Herceptin Fab, the source of the contact set; PMID:31042744 — cryo-EM of the HER2-trastuzumab-pertuzumab complex
- **evidence_against:** EXP-002 — 1N8Z alone labels 7 CDRH3 positions contact-free, but across 8 trastuzumab-HER2 structures only 3 never contact HER2. 1N8Z has unmodelled HER2 regions (102-110, 303-305, 361-364, 581-590) and several CDRH3 positions contact residues in 579-593, so a single-structure contact set under-counts contacts in the direction that flatters the design
- **query_terms:** trastuzumab HER2 contact residues 1N8Z interface paratope; alanine scanning trastuzumab HER2 hot spot energetic contribution; buried surface area shape complementarity HER2 antibody interface
- **history:**
  - 2026-11-14 created, confidence low — we have the structure but have not yet extracted a per-residue energetic contribution, and a contact is not the same as a hot spot (pass-001)
  - 2026-11-14 **retired** (EXP-002). Not downgraded but retired: the claim names 1N8Z as the authority, and 1N8Z is demonstrably insufficient. Superseded by TH-009 (which positions are safe) and TH-010 (how much of the interface CDRH3 even covers). Kept rather than deleted because the failure mode — concluding from one crystal structure with gaps — is the lesson.

### TH-009 — The <70% CDRH3 identity gate cannot be met without mutating at least one HER2-contacting position, and the cheapest one to sacrifice is Y109

- **confidence:** high
- **status:** active
- **property:** epitope; binding
- **falsified_by:** A design substituting only A106/M107/D108 passes abguard's CDRH3 identity check (it should not: 3/11 = 72.7%), or a design that also substitutes Y109 loses competitive binding against trastuzumab while one sparing Y109 retains it.
- **evidence_for:** EXP-002 — CDRH3 is WGGDGFYAMDY (11 aa, 99-109), so <70% identity needs 4 substitutions (4/11 = 63.6%); across 8 trastuzumab-HER2 structures at 4.5 A only A106, M107 and D108 never contact HER2 (0/8 each), while Y109 contacts in 1/8, G101 2/8, G100 3/8, F104 4/8, W99 5/8, D102 6/8, G103 7/8, Y105 7/8
- **evidence_against:** Contact frequency across crystal structures is a geometric proxy, not an energetic one. A 1-of-8 contact may still be energetically important, and no alanine scan has been run to rule that out (IDEA-008)
- **query_terms:** antibody paratope contact conservation across crystal structures; CDR position energetic hot spot alanine scan contribution; crystal packing artefact antibody antigen contact variability
- **history:**
  - 2026-11-14 created at confidence **high** (EXP-002). Unusually high for a new theory because the claim is arithmetic plus measured geometry across 8 independent structures, not an inference: 4 substitutions are required and only 3 safe positions exist. The design recipe A106/M107/D108 + Y109 follows directly. What is *not* high-confidence is the assumption that sacrificing a 1-of-8 contact is harmless, which is why the theory names Y109 explicitly so it can be falsified.
  - 2026-11-14 **no change**, and that is the finding (pass-003). Three targeted searches plus the full cache turned up nothing on the energetic contribution of trastuzumab CDRH3 Tyr109, or of any individual trastuzumab paratope residue. The literature cannot settle whether sacrificing Y109 is safe, so IDEA-010 is escalated to high priority: if we want that answer we have to compute it. The theory stands on geometry alone, which is exactly as strong as it was before.
  - 2026-11-14 **independently corroborated** (EXP-007). Our 4-substitution minimum came from the 11-residue CDRH3 `WGGDGFYAMDY` (4/11 = 63.6%). A teammate reached the same minimum from the 10-mer identity arithmetic (3/10 = exactly 70%, which fails “strictly under”) and corrected the pre-filter's design band from >=3 to >=4. Two independent routes, one answer. Both of us had carried >=3 for hours: root `AGENTS.md` warns about an off-by-one, and neither of us applied that warning to the *identity* rule — only to the mutation-count rule it was written about.
  - 2026-11-14 **strengthened and softened at once** (EXP-008). Strengthened: across the 8 tightest legal variants in the measured set, **M107 is mutated in 8 of 8 and A106 in 7 of 8** — exactly the two positions this theory identified as never contacting HER2. Measured affinity and structural contact analysis converge on the same positions from independent directions. Softened: the single legal variant that beats trastuzumab (`WGGSGPYLFD`, 0.095 nM) also changes **D102 and F104**, which contact HER2 in 6 of 8 and 4 of 8 structures. So touching a contact is a risk to price, not a veto — and D102 is mutated in 8 of 8 tight legal variants, making it a better fourth substitution than Y109, which the SPR table cannot even see.

### TH-010 — CDRH3-only redesign can address at most a third of the trastuzumab-HER2 interface, so CDRH3 is the wrong place to look for large binding gains

- **confidence:** high
- **status:** active
- **property:** binding
- **falsified_by:** CDRH3-only designs show SPR improvements over trastuzumab comparable to what CDR-wide or framework-inclusive redesign achieves.
- **evidence_for:** EXP-002 — per-CDR share of interface contacts across 8 structures: H3 12-34% (26% in 1N8Z), H2 14-47%, L3 11-29%, L1 5-17%, L2 0-30%, framework up to 23% (6BGT); heavy-chain CDRs carry 55% and light-chain CDRs 45% in 1N8Z
- **evidence_against:** Interface share counts contacts, not binding energy; a minority of contacts can carry a majority of the affinity. PMID:8918938 got a 1230-fold KD gain from VL-CDR3 and VH-CDR3 libraries alone, which shows CDR3s can dominate affinity even when they do not dominate contact count
- **query_terms:** antibody interface contact share CDR contribution binding energy; CDR H2 L3 contribution affinity antigen; framework contact antigen antibody interface contribution
- **history:**
  - 2026-11-14 created at confidence **high** for the measurement, which is unambiguous across 8 structures (EXP-002). The *strategic* reading is weaker and the counter-evidence says so: PMID:8918938 is a direct case of CDR3 libraries delivering a 1230-fold gain. Practical consequence: all four HER2 RFD4 weights redesign CDRH3 only, so the tool set available to us works on a minority of the interface. If we want more, it needs `antibody-2026-10-30` on other CDRs (IDEA-009).

### TH-011 — Refold scores cannot separate close CDRH3 designs, because predicted conformation for an 11-residue H3 loop is only accurate to roughly 2 A

- **confidence:** high
- **status:** active
- **property:** binding
- **falsified_by:** Repeated folding of the same design across seeds and models gives min_ipae spreads much smaller than the score gaps we are ranking on, or the returned SPR values track the predicted ranking.
- **evidence_for:** PMID:40814020 — **AlphaFold3**, a modern co-folding model of the same class as Boltz-2, achieves median unbound CDR-H3 RMSD **2.9 A for antibodies using twenty seeds**; PMID:37248282 — ABodyBuilder2 predicts CDR-H3 at **2.81 A** on 34 recently solved antibodies, only 0.09 A better than AlphaFold-Multimer; PMID:35099535 — ABlooper averages **2.49 A** CDR-H3 RMSD on the Rosetta Antibody Benchmark; PMID:19062174 — median backbone RMSD **2.4 A for medium (10-11 residue) loops**, and trastuzumab's CDRH3 is 11 aa; EXP-001 — our two refolded designs differed by just **0.021 min_ipae on a single Boltz-2 seed**, a gap with no established error bar
- **evidence_against:** PMID:35099535 — ABlooper's RMSD drops from 2.49 A to **2.05 A on its 75% most confident predictions**, so a model's own confidence estimate does carry signal. Confidence-based *filtering* (discard the bottom quartile) is defensible even where fine *ranking* between close designs is not. Separately, all these benchmarks measure unbound or apo loop geometry, not interface score reproducibility, which is what we actually rank on — the inference from one to the other is ours, not the papers'
- **query_terms:** CDRH3 loop conformation prediction accuracy RMSD benchmark; Boltz AlphaFold interface pAE reliability antibody antigen; folding model seed variance interface confidence score
- **history:**
  - 2026-11-14 created, confidence medium (pass-003). Found while looking for something else: the question "can we trust the ranking?" matters more than any single score, because only ~4 designs are tested and rank order is the real output. The two RMSD studies are older methods than our folding stack, hence medium not high. The cheap test is IDEA-006, which already exists: add a second Boltz-2 seed and see whether the 0.021 gap in EXP-001 survives.
  - 2026-11-14 medium → **high** (pass-004). The "older methods" caveat is gone: PMID:40814020 benchmarks **AlphaFold3**, the same class of co-folding model as Boltz-2, at 2.9 Å median CDR-H3 RMSD **with twenty seeds**. Three further independent benchmarks land in the same 2.5–2.9 Å band (PMID:37248282, PMID:35099535, PMID:19062174). Four methods spanning 2009–2025 converge, so this is a property of the problem, not of one tool. Consequence: ranking designs on a 0.021 min_ipae gap is indefensible. What survives is confidence-based *filtering*, per the ABlooper 75% result.
  - 2026-11-14 **confirmed by measurement, caveat closed** (EXP-003). pass-004 left one inference as ours rather than the literature's: that loop-RMSD error implies interface-score noise. EXP-003 measures the thing itself, on our target and our models. Boltz-2 ipTM separates real from shuffled HER2 pairs at AUROC 0.645, and 94.3% of 8,190 shuffled pairs score >=0.7. Stays `high`, now resting on primary measurement instead of an analogy. Superseded in practice by TH-013, which states the consequence.

### TH-012 — A VHH submission is better supported by our tooling than an scFv, on two independent legs: thermostability is predictable and structure prediction is more accurate

- **confidence:** medium
- **status:** active
- **property:** thermostability; binding
- **falsified_by:** Our VHH designs come back with worse measured Tm or worse SPR than our scFv designs, or the sdAb-trained predictors rank them no better than chance.
- **evidence_for:** PMID:39154093 — TEMPRO predicts **sdAb** Tm with MAE 4.03 C, RMSE 5.66 C, and there is no equivalent Tm-trained predictor for scFv (see TH-007); PMID:40814020 — AlphaFold3 median unbound CDR-H3 RMSD is **2.2 A for nanobodies versus 2.9 A for antibodies**, so the format we can predict best is also the one we can rank best; PMID:37632647 — an anti-HER2 VHH reaches sub-nanomolar affinity on trastuzumab's own epitope, so the route is not hypothetical
- **evidence_against:** PMID:37248282 — NanoBodyBuilder2 predicts nanobody CDR-H3 at 2.89 A against ABodyBuilder2's 2.81 A for antibodies, i.e. **no nanobody advantage at all** in that benchmark. So the prediction-accuracy leg rests on one study (PMID:40814020) and is contradicted by another. A VHH is also a single domain, so it forgoes the light-chain CDRs that carry 45% of trastuzumab's interface contacts (TH-010/EXP-002)
- **query_terms:** nanobody VHH structure prediction accuracy benchmark CDR3; single domain antibody developability expression yield comparison scFv; nanobody affinity maturation HER2 epitope
- **history:**
  - 2026-11-14 created, confidence medium (pass-004). A strategic rather than mechanistic theory, assembled from two findings that arrived in different passes and only mean something together: TEMPRO is sdAb-only (pass-002) and AF3 predicts nanobody H3 better (pass-004). Deliberately held at medium, not high, because the two benchmarks disagree about whether nanobody prediction is actually more accurate, and because a single domain gives up the light-chain contacts EXP-002 measured. The decision it informs is real though: IDEA-005 (submit a VHH alongside the scFv designs) is now high priority.

### TH-013 — Fold-confidence scores cannot rank HER2 CDRH3 designs, for any template-taking model tried, and the topmodel score is actively misleading

- **confidence:** high
- **status:** active
- **property:** binding
- **falsified_by:** A decoy control on our own designs shows min_ipae or ipTM separating real HER2 binders from shuffled pairs at AUROC well above 0.8, or the returned SPR values track our min_ipae ranking, or a folding model we have not yet tried separates measured binders at AUROC well above 0.65 on the same 60-variant calibration set.
- **evidence_for:** EXP-006 — **Protenix is no better**: across ten of its metrics, AUROC spans only 0.426–0.526 on the same 60 measured variants, with the chain-pair `px_iptm_ab_her2` at 0.448. EXP-004 — calibrated against **measured KD** on 60 SPR variants (30 with KD 0.023–13,863 nM, 30 non-measurable): ipTM AUROC **0.502**, pTM 0.504, pLDDT 0.356, **min_ipae 0.497**, mean_ipae 0.454, against sequence-only p_binder **0.849**; no fold metric tracks KD among the binders (all |Spearman| ≤ 0.21, p ≥ 0.28); adding fold features made the classifier *worse* (0.849 → 0.806); **ipTM ≥ 0.7 passes 100% of true non-binders**. EXP-003 — real-vs-shuffled AUROC for ipTM: boltz2 **0.645**, af3 0.754, chai1 0.568; for the topmodel score: af3 **0.397**, boltz2 0.443, chai1 0.492, i.e. at or below random **with the direction inverted** (shuffled medians higher than real). Boltz-2 ipTM real median 0.911 vs shuffled 0.884, and **94.3% of 8,190 shuffled pairs score >=0.7**. Spearman between topmodel and ipTM is -0.182 (af3) and -0.658 (chai1), so the two scores disagree with each other
- **evidence_against:** **Scope, not exception:** the negative holds for *fixed-framework CDRH3 optimisation*. For **de novo** designs, where the whole binder changes, AF3 ipTM retains signal (AUROC 0.754, EXP-003), so this theory must not be read as "structure scores never work". EXP-004 also bounds itself: n=60 rules out a strong effect (AUROC > ~0.65), not a weak one. The earlier clause here — that a structure score was still worth keeping as a coarse filter — is **retracted**: measured, ipTM ≥ 0.7 passes 100% of true non-binders, so it removes nothing on this track.
- **query_terms:** interface predicted aligned error decoy discrimination benchmark; ipTM AlphaFold screening binder hit rate validation; decoy control protein complex prediction score calibration
- **history:**
  - 2026-11-14 created at confidence **high** (EXP-003). High because it is a direct measurement on our target with our models over thousands of decoys, not an inference. It invalidates the basis of EXP-001's ranking: ipTM 0.946 and 0.965 sit inside the band 94.3% of shuffled decoys occupy. Practical consequence: stop ranking on min_ipae, keep a structure score only as a coarse filter, and never trust a new score until its decoy control has been run - which costs nothing.
  - 2026-11-14 **no change, and the absence is the finding** (pass-005). Three targeted searches for published validation of ipTM or ipAE *as a binder screen* returned nothing on the question: loop-accuracy benchmarks, a CASP complex-prediction paper, and a small-molecule docking review. So the team's own decoy control (EXP-003) is currently better evidence on this than anything the literature offered, which makes IDEA-012 (run the decoy control before trusting any score) more valuable, not less.
  - 2026-11-14 **strengthened, and one clause retracted** (EXP-004). EXP-003 used shuffled pairs, which are easy negatives. EXP-004 calibrates against *measured KD* and every fold metric lands at chance, `min_ipae` included. The compromise I had kept in `evidence_against` — that a structure score was still worth holding as a coarse filter — is measured false for this track: ipTM ≥ 0.7 passes **100%** of true non-binders. Retracted in place rather than quietly dropped. The “only 10 of 336 tokens vary” explanation does not cover `min_ipae`, whose spread is 3.0, so this is not a resolution artefact.
  - 2026-11-14 **broadened from one model to the class** (EXP-006, ours). The obvious response to EXP-004 was “try a different folder”, since `AGENTS.md` prefers template-taking models. We tested that on the teammate’s committed two-model table and **Protenix fails identically**: all ten of its metrics land between AUROC 0.426 and 0.526, including `px_iptm_ab_her2`, the antibody↔HER2 chain-pair ipTM, which is the most targeted score available and sits *below* chance. Confidence stays `high` but the scope is wider and better earned. Practical value: this closes a GPU-hungry line of enquiry for two minutes of analysis — nobody need refold with OpenFold3 or Mol-1 hoping for a screen.

### TH-014 — A one-hot sequence model ranks trastuzumab CDRH3 variants well enough to choose submissions, and should replace structure scores for the optimisation track

- **confidence:** high
- **status:** active
- **property:** binding
- **falsified_by:** The measured SPR values Adaptyv returns do not correlate with the sequence model's predicted ranking, or its advantage disappears on designs further than 4 mutations from anything in the training data.
- **evidence_for:** EXP-007 — revalidated inside the **legal >=4-substitution band**: classifier AUROC **0.913** (n=556), ranker Spearman **+0.628** (n=131), and a top-20 median true KD of **6.5 nM** against an 86 nM set median. EXP-003 — overall figures with Hamming<=2 clusters held out together: AUROC **0.927**, Spearman **+0.687**, recovering the 0.023 nM best binder in the set. Out-of-fold, keeping the top 10 gives **14 nM against 74 nM for a random pick of 10** (`NICCOLO/eda/AGENT_GUIDE.md`, commit fd1b002). Model form: logistic regression L2 `C=0.3` for the classifier, ridge `alpha=3.0` for the ranker; gradient boosting and random forests were tried and lost on both tasks, because the signal is close to additive across positions. Mutation count alone gives 0.800 AUROC and +0.387 Spearman, so the model is doing more than counting mutations
- **evidence_against:** Labels exist only for **10-residue CDRH3 variants on the trastuzumab framework**, so none of this transfers to the de novo track, where there are no labels at all. The published clonotype splits also leave **57% of test sequences 1 mutation from a training sequence**, which is why only the >=3-mutation figures should be quoted. The best-residue estimates behind it rest on n=4-5 at several positions (TH-015)
- **query_terms:** one-hot sequence model antibody affinity prediction generalisation; antibody variant effect prediction held out cluster validation; deep mutational scanning affinity model extrapolation distance
- **history:**
  - 2026-11-14 created at confidence **high** (EXP-003). The leakage controls are what earn the confidence: grouped clusters, a constraint-band report, and FACS/SPR overlap removed before fitting. Quote the in-band figures, not the raw 0.99s, not the raw 0.99 figures from the ungrouped FACS classifier. Consequence: IDEA-011 adopts this as the primary ranker for optimisation-track designs.
  - 2026-11-14 **numbers corrected, claim intact** (EXP-007). Two errors this ledger had propagated. (1) The band was wrong: we used >=3 substitutions, from `abguard`'s “>2 mutations” rule, but on a **10-residue** CDRH3 three substitutions leave exactly 70% identity and the rule demands *strictly under*, so the minimum is **4**. In-band figures become AUROC **0.913** and Spearman **+0.628**. (2) Worse, the **1.42 nM top-20 median was never achievable** — that top 20 is full of 1–2 substitution variants the identity rule forbids. Inside the legal band it is **6.5 nM**. The claim survives: 6.5 nM is still ~13x better than the 86 nM set median, and the model barely degrades in band (0.927 → 0.913), so the loss is imposed by the constraint, not the model. The cost itself is TH-018.

### TH-015 — Measured per-position SPR preferences give a data-driven CDRH3 recipe, but the best residues rest on very small samples

- **confidence:** medium
- **status:** active
- **property:** binding
- **falsified_by:** Designs built from the per-position favourites come back no tighter than trastuzumab, or the preferences fail to reproduce when n is increased.
- **evidence_for:** EXP-003 — mean log10 KD by residue among quantified variants (lower is tighter), for CDRH3 positions 1-10 = W99-D108: pos1 **F 0.92** (n5) beats wild-type W 1.25 (n104); pos4 **Q 0.15** (n4) and **E 0.18** (n4) against wild-type D absent from the top 3; pos8 **L 0.64** (n4) against wild-type A; pos9 M 0.72 (n31) is wild-type and best; pos10 **V 1.05** (n4) against D 1.50 (n120); pos6 F and pos7 Y are wild-type and best
- **evidence_against:** The standout substitutions all have **n=4-5**, against n=104-128 for the wild-type residues they beat, so each is a mean of a handful of measurements and could be noise or clone-specific. Position 11 (Y109) is **not in the dataset** at all, since SPR-768 varies a 10-mer, and Y109 is exactly the position TH-009's recipe sacrifices
- **query_terms:** deep mutational scanning CDR position substitution affinity matrix; site saturation mutagenesis antibody affinity per position measured
- **history:**
  - 2026-11-14 created, confidence medium (EXP-003). Medium, not high, purely because of sample size: the direction of each preference is measured, but the magnitudes rest on n=4-5. Still the most directly actionable data we have, because it covers 10 of our 11 CDRH3 positions with real SPR rather than a prediction. Note the awkward gap: the one position it cannot speak to is Y109.

### TH-016 — NanoMelt is the only usable Tm predictor we have, so thermostability is rankable for a VHH and out-of-domain for a VH or scFv

- **confidence:** high
- **status:** active
- **property:** thermostability
- **falsified_by:** The Tm values Adaptyv returns correlate with NanoMelt's ranking just as well for our scFv designs as for a VHH, or a Tm-trained predictor for paired-chain Fv turns up.
- **evidence_for:** EXP-005 — neither platform ships a Tm predictor (Nebius `stab-ddg` is StaB-ddG, i.e. **binding** ddG not folding stability; Foundry's AQ workflows are SandboxAQ small-molecule tools); cheap composition/charge/hydrophobicity features reach only Pearson **+0.40**, MAE 7.6 C against 8.3 C for predicting the mean; NanoMelt (mAbs 2025) has published held-out Pearson **0.853**, MAE **4.1 C** on 639 measured **nanobody** Tm, is installed on the VM, reproduces its own reference output exactly, and runs at ~45 sequences/s. Its training data also independently encodes the cysteine rule `abguard` enforces: mean Tm 41 C at 0 cysteines, 67 C at 2, 73 C at 4, 53 C at 5 (odd implies unpaired)
- **evidence_against:** Applied to trastuzumab's VH the model shows clear domain artefacts: **11 of its top 12 single mutations are ->Glu** (the only substitution with a positive mean effect across 1,577 mutants, i.e. it has learned "more acidic = more stable"), and **`V37F`, a VHH hallmark position**, is among the top hits — the model partly wants to nanobody-ify the VH. So its 4.1 C MAE is an in-domain figure and should not be quoted for a VH or scFv
- **query_terms:** NanoMelt nanobody melting temperature predictor validation; Fv scFv melting temperature predictor paired chain training data; protein stability predictor domain shift antibody format
- **history:**
  - 2026-11-14 created at confidence **high**, superseding TH-007 (EXP-005). High on the tooling facts, which are verified rather than inferred: the tool reproduces its reference output, and the alternatives were checked and are not Tm predictors at all. The part that stays uncertain is deliberately in the claim itself — "out-of-domain for a VH or scFv" — rather than hidden in a caveat. A predictor's **training domain is part of its error bar**, and that is the transferable lesson. Third leg under TH-012.

### TH-017 — Framework stabilisation is real but bounded: no single mutation clears the error bar, stacking collapses, and the best hits are artefacts

- **confidence:** medium
- **status:** active
- **property:** thermostability
- **falsified_by:** A stacked framework variant comes back with a measured Tm gain above ~4 C that survives the charge, cysteine and position-98 filters.
- **evidence_for:** EXP-005 — full non-CDR saturation scan, 1,577 single mutants: mean dTm **-0.59 C**, only **27.6% stabilising at all**, **0 above +4.1 C** (best single `R98A` at +3.74 C, still inside the error bar), 19 below -4.1 C. Combinations are strongly sub-additive then harmful: stack3 reaches +4.28 C where additivity predicts +6.23, and stack12 **collapses to +1.44 C** where additivity predicts +18.60 C. Random 6-mutation subsets beat greedy ordering, averaging +4.16 C with a best of +6.09 C, and 12 of 29 combinations clear +4.1 C
- **evidence_against:** The 12-of-29 combinations that do clear the error bar are a real signal, so "bounded" is not "impossible". But they are dominated by the ->Glu charge bias, and the scan's most fragile positions (W36, C22, C96, L20, M83, V48, D73, E89) plus position 98 at the CDRH3 junction have to be excluded by hand, so the usable search space is much smaller than 1,577
- **query_terms:** framework mutation epistasis sub-additive stability antibody; net charge mutation protein stability expression tradeoff; consensus mutation antibody VH stabilisation measured Tm
- **history:**
  - 2026-11-14 created, confidence medium (EXP-005). Medium because every number is a *prediction* from a model that is out-of-domain for this VH (TH-016), so the magnitudes are soft even though the pattern — sub-additivity, charge bias, greedy losing to random — is consistent and mechanistically plausible. Practical reading: framework stabilisation is worth one submission slot, not the whole strategy, and any candidate must pass the cysteine check and avoid position 98.

### TH-018 — The identity constraint costs roughly 4x in achievable affinity, because the tightest known binders are 1-2 substitutions from trastuzumab and are therefore illegal

- **confidence:** high
- **status:** active
- **property:** binding
- **falsified_by:** A design at >=4 substitutions comes back with a measured KD at or below the ~1.4 nM that the unconstrained top 20 achieves, or the in-band and out-of-band affinity distributions turn out not to differ.
- **evidence_for:** EXP-007 — ranking the 179 quantified SPR variants and keeping the top 20 gives a median true KD of **1.42 nM**, but that top 20 is dominated by 1-2 substitution variants `abguard` forbids; inside the legal **>=4-substitution** band the top 20 medians **6.5 nM**, a ~4.6x loss. Model quality barely moves across the same boundary (AUROC 0.927 -> 0.913, Spearman +0.687 -> +0.628), so the loss is imposed by the rule and not by a failure of ranking
- **evidence_against:** 6.5 nM is still **~13x better than the 86 nM set median**, so the constraint raises the floor we must clear rather than making the task hopeless. And the comparison that actually scores is against **measured trastuzumab**, not against the dataset's best variant — if trastuzumab's own KD sits above 6.5 nM, an in-band design can still win on binding. **EXP-008 settled this: trastuzumab's own CDRH3 measures 0.152 nM**, so the legal band's 6.5 nM top-20 median is ~43x *worse* than the bar we must clear. One legal variant of 21 does beat it (0.095 nM, 1.6x), so the task is demonstrated rather than hopeless — but narrow
- **query_terms:** affinity maturation diminishing returns sequence distance wild type; novelty constraint antibody design affinity tradeoff measured
- **history:**
  - 2026-11-14 created at confidence **high** (EXP-007). High because it is arithmetic over measured KD values, not a prediction: the in-band and out-of-band top-20 medians are both measured, and the model's own performance is near-identical across the boundary, which isolates the constraint as the cause. The strategic reading is the valuable part and belongs in the pitch: the designs we are *permitted* to submit start from a worse part of the affinity landscape than the dataset's best variants occupy. Saying that explicitly, with 6.5 nM against 1.42 nM, is a better judgment story than quoting the unconstrained figure. Immediate follow-up: get trastuzumab's own measured KD from the SPR table, since that is the number we are actually scored against (IDEA-017).

### TH-019 — Beating trastuzumab on binding inside the legal band is demonstrated but rare and narrow, so binding alone is a weak basis for an S-tier claim

- **confidence:** high
- **status:** active
- **property:** binding
- **falsified_by:** A substantial fraction of legal (>=4 substitution) designs come back with measured KD below 0.152 nM, or one of ours beats it by more than ~2x.
- **evidence_for:** EXP-009 (supersedes EXP-008, full 179-variant table) — among **131** legal variants with a quantified KD, **1 beats trastuzumab's 0.152 nM (0.8%)** and **0 beat its best replicate, 0.029 nM**; the one winner `WGGSGPYLFD` (0.095 nM, 1.6x) is itself in the provided data and so unsubmittable, as are all 131. Trastuzumab's 8 replicate rows span 0.029-0.178 nM, so the baseline must be named. Earlier figure from EXP-008, now superseded: among 21 legal variants, **exactly 1 beats it** (`WGGSGPYLFD`, 0.095 nM, 1.6x), the second best is 0.368 nM (2.4x *worse*), and the legal-band median is 127.5 nM, roughly 840x worse. The set's best binder overall sits at 3 substitutions (0.023 nM, 6.5x better) and is forbidden by the identity rule
- **evidence_against:** EXP-009 closed the sampling objection by running the full table; the rate came out **lower**, so the theory is not resting on a biased subset any more. The pre-filter also does better than chance at finding these: its in-band top 20 medians 6.5 nM against a 127.5 nM band median, so a ranked search is not the same as a random draw
- **query_terms:** affinity maturation ceiling wild type antibody improvement measured KD; novelty constrained antibody design measured affinity outcome
- **history:**
  - 2026-11-14 **revised** by EXP-009 on the full 179-variant table: the rate falls from 1 of 21 (4.8%) to **1 of 131 (0.8%)**, and to **0 of 131** against trastuzumab's best replicate. Confidence stays **high** and the claim strengthens — binding is not merely a rare win inside the legal band, it is effectively absent, and the single winner cannot be submitted. The strategic consequence recorded below is reinforced, not softened.
  - 2026-11-14 created at confidence **high** (EXP-008). High because every number is a measured KD and the comparison is arithmetic. The strategic consequence matters more than the number: ranking after testing counts **S-tier** designs that improve all three properties first, and this says binding is the hardest of the three to improve legally. A design whose case rests only on predicted binding is unlikely to be S-tier. That argues for spending the submission slots on designs with an independent second property story — framework-stabilised (IDEA-015) or VHH with in-domain Tm prediction (IDEA-005) — rather than on four variations of the same affinity bet.

### TH-020 — CDRH3 positions divide into beneficial, inert and costly, and contact frequency does not tell you which

- **confidence:** high
- **status:** active
- **property:** binding
- **falsified_by:** On measured affinity for legal variants, mutating G100 or G101 turns out neutral or harmful, or mutating M107 or F104 turns out neutral or beneficial, at n comparable to 131 and with multiple-comparison correction.
- **evidence_for:** EXP-009 — Mann-Whitney on log10 KD across all 131 legal quantified variants, Benjamini-Hochberg corrected over ten positions. Costly to mutate: Y105 (q<0.001), W99 (q<0.001), F104 (q=0.028), M107 (q=0.040), G103 (q=0.090). Beneficial: G100 (q=0.062), G101 (q=0.077). Inert: A106, D102, D108 (q>0.83). Crossed with EXP-002's contact table, **M107 never contacts HER2 in any of 8 structures and is still the costliest position to mutate**, while G100/G101 do contact (3/8, 2/8) and are the two that help
- **evidence_against:** G100 and G101 sit above the conventional 0.05 threshold, so they are the best available bet rather than a proven effect. M107's conserved arm is n=11, so its effect size is poorly determined even though two independent estimators agree on the sign. All 131 variants are in the provided data, so this is a map of a measured landscape, not evidence any new design will behave the same way
- **query_terms:** CDR position specific mutation affinity contribution paratope periphery; antibody contact residue energetic hot spot versus structural contact
- **history:**
  - 2026-11-14 created at confidence **high** (EXP-009). High because it is measured KD with correction over all ten positions, and because the direction at M107 is confirmed by two independent estimators. It carries two consequences. First, IDEA-008's recipe inverts: spend the forced substitutions on **G100 + G101 + A106 + D108**, not on M107 or F104. Second, a methodological rule the ledger should keep: **a "winners share X" claim is uninterpretable without X's frequency among non-winners** — M107 looked like a winning pattern only because it is mutated in 92% of everything. A corollary learned the same day: compare any new predictor against the existing model **out-of-fold**, never an in-sample univariate correlation against a cross-validated score
