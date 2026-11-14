# Ideas — the backlog for the de novo VHH / HER2 campaign

Same schema as [`../../../LEO/tasks/experiment-ledger/IDEAS.md`](../../../LEO/tasks/experiment-ledger/IDEAS.md).
An idea becomes an experiment by getting an `exp_id`; it dies by getting
`status: dropped` **with a reason**, which is itself a result worth keeping.

`status`: `open` | `lit-checked` | `running` | `done` | `dropped`.
`source`: `literature` | `experiment` | `team`.

**Id policy.** Imported ideas **keep their original ids** (IDEA-001 … IDEA-007)
and carry `imported_from`. Ids new to this campaign start at **IDEA-008**, the
first id Leo's ledger has not used.

---

### IDEA-001 — Fix the VHH framework to caplacizumab and spend the whole design budget on the three CDRs

- **status:** lit-checked
- **priority:** high
- **tests:** TH-005
- **cost:** no extra GPU — it is a choice of input spec, not an extra run
- **source:** literature
- **imported_from:** LEO/tasks/experiment-ledger IDEA-001
- **imported_date:** 2026-11-14
- **evidence:** PMID:40653772; `/home/hackathon/hackathon-team-3/apheris/rfd4/Data/03-vhh/pmvk_caplacizumab.cif.zst`
- **notes:** **Re-scoped on import, and this is the biggest change any import needed.** Leo's IDEA-001 was "CDRH3 diversity for binding, framework mutations for Tm, in the same molecule" and assumed `antibody-2026-10-30` could design framework residues. It cannot: every documented RFD4-Proteina antibody route keeps the framework fixed and designs only CDR loops. Framework redesign is undocumented on this platform. So the surviving, actionable half is **framework choice**: take the supplied caplacizumab VHH framework (marketed, humanised, 86 positions, three CDR sentinels at label_seq_id 23/38/76) and keep it fixed, buying stability by selection rather than by design. The stability half of Leo's idea is therefore **not** available to us and the pitch must not claim it.

### IDEA-002 — Cap arginine content in the designed CDRs and measure what it costs in predicted affinity

- **status:** open
- **priority:** medium
- **tests:** TH-003
- **cost:** post-hoc filter over an existing design set, no GPU
- **source:** literature
- **imported_from:** LEO/tasks/experiment-ledger IDEA-002
- **imported_date:** 2026-11-14
- **evidence:** PMID:28778924
- **notes:** Cheap and diagnostic. Record Arg fraction for every design *before* filtering, then check whether Arg fraction correlates with min_ipae rank. A positive correlation is a warning about the metric, not about the designs: it would mean our ranking rewards non-specific binding, which is the exact failure PMID:28778924 describes. Threshold to apply: <= 10% of CDR residues, <= 2 Arg in CDR3 (OUR JUDGEMENT — the source gives a mechanism, not a number).

### IDEA-003 — Commit to TEMPRO plus the Tm > 65 C criterion as the VHH thermostability proxy

- **status:** done
- **priority:** high
- **tests:** TH-005
- **cost:** one afternoon of setup; already paid, see notes
- **source:** team
- **imported_from:** LEO/tasks/experiment-ledger IDEA-003
- **imported_date:** 2026-11-14
- **evidence:** PMID:39154093; PMID:41564239; PMID:36683173
- **notes:** **Closed for this campaign, and closed favourably — the format question that blocked it in Leo's ledger does not apply to us.** Leo's IDEA-003 was split by format: TEMPRO and PMID:41564239 are single-domain-antibody-specific, so they applied only to a VHH submission and an scFv had no Tm-trained predictor. **We are submitting a VHH, so we get the good half.** Decision: rank on **TEMPRO** (PMID:39154093, MAE **4.03 C**, RMSE **5.66 C** on sdAb Tm — an error bar on exactly the quantity Adaptyv measures), gated by the **Tm > 65 C** screening criterion from PMID:41564239 (NBsTem_Tm, trained on 514 experimental Tm values, Pearson **0.83** on external test; the paired NBsTem_Q theoretical model has accuracy **0.84**; the full published criterion is "Tm > 65 C **and Qclass IV**"). **Obtainability verified by installing and running it, not by reading about it:** repo is `github.com/Jerome-Alvarez/TEMPRO` (first author of PMID:39154093), 43 MB, pretrained Keras weights ship in `user.zip`, no separate download, needs Python 3.9 + TF 2.10.1 + fair-esm; a real VHH sequence was pushed through end to end and returned **Tm 75.55 C**. Two caveats to carry: the repo has **no LICENSE file** (all rights reserved by default — fine for internal ranking, flag before redistributing anything derived from it), and we can run only the `ESM_650M` head, since `ESM_15B` needs ~58 GB of embeddings the author never uploaded. We can reproduce TEMPRO's Tm but **not** the `Qclass IV` half of PMID:41564239's criterion, which needs their NBsTem_Q model or its web server — so treat **Tm > 65 C as a necessary-not-sufficient gate** and say so in the pitch.

### IDEA-004 — Use the 1N8Z contact atoms as hotspots, and stop there

- **status:** lit-checked
- **priority:** medium
- **tests:** TH-008
- **cost:** none — reuse the supplied HER2 spec's `C_HOT` block
- **source:** literature
- **imported_from:** LEO/tasks/experiment-ledger IDEA-004
- **imported_date:** 2026-11-14
- **evidence:** PMID:12610629; `/home/hackathon/hackathon-team-3/apheris/rfd4/Data/12-her2/her2_cdrh3.yaml`
- **notes:** **Re-scoped on import.** Leo's IDEA-004 was "extract per-residue energetic contributions at the 1N8Z interface" in order to score epitope preservation. We are not doing that: no cached record gives per-residue energetics at this interface, and we have no PyRosetta or FoldX installed to compute them. What we take instead is the cheap, cited half — the 4 `C_HOT` atoms already in the supplied HER2 spec, anchored to "the HER2 atoms nearest CDRH3 in PDB 1N8Z". A contact is not a hot spot, so those atoms are good enough to *aim* generation and not good enough to *score* it. Epitope preservation is scored by the competition assay, full stop.

### IDEA-005 — The de novo VHH arm itself: generate VHHs against the trastuzumab epitope on HER2 domain IV

- **status:** open
- **priority:** high
- **tests:** TH-006; TH-009; TH-011
- **cost:** 100-200 designs per arm at ~24 s per VHH design on an A100, so ~40-80 min GPU per arm, plus refolding at ~25 s per design per model
- **source:** literature
- **imported_from:** LEO/tasks/experiment-ledger IDEA-005
- **imported_date:** 2026-11-14
- **evidence:** PMID:37632647; PMID:41193805
- **notes:** **This is now the campaign, not a hedge.** In Leo's ledger IDEA-005 was "submit a VHH alongside the scFv designs to hedge the track split"; here it is the whole brief. Two feasibility priors, both primary: PMID:37632647 (a real anti-HER2 VHH at sub-nM binds trastuzumab's epitope) says the target is reachable by a single domain, and PMID:41193805 (Nature 2025) says de novo VHHs against user-specified epitopes work at all. **There is no HER2 VHH example and no HER2 VHH weight on this platform** — the HER2 fine-tunes only redesign trastuzumab's CDRH3 on the Fv. So this arm needs `antibody-2026-10-30` plus a spec we author: HER2 mmCIF (PDB not accepted) + the caplacizumab VHH framework + the 4 `C_HOT` epitope atoms + three CDR length ranges. Write the hypothesis row **before** the run.

### IDEA-007 — Check our top designs against claimed patent sequences before we submit

- **status:** open
- **priority:** low
- **tests:**
- **cost:** one patentcore pass, currently unaffordable — see notes
- **source:** literature
- **imported_from:** LEO/tasks/experiment-ledger IDEA-007
- **imported_date:** 2026-11-14
- **evidence:**
- **query_terms:** patent: anti-HER2 VHH single domain antibody sequence claim; patent: nanobody HER2 domain IV epitope claim; patent: computationally designed antibody CDR freedom to operate
- **notes:** **Blocked on budget, not on method: `AMASS_API_KEY` is unset and the credits are unclaimed, so no new patentcore search can run.** The three cached patentcore files already show a crowded field (`Anti-HER2 antibody variants`, `High-affinity fully-humanized anti-Her2 monoclonal antibody`, hersintuzumab and ~35 more). **One concrete, zero-cost action is available now and should be taken:** PMID:37632647 publishes the sequence of a sub-nM anti-HER2 VHH that binds our exact epitope, so a design that lands close to it is both a freedom-to-operate risk and an `abguard` novelty risk (the de novo track requires **under 70% identity to any known antibody**, CITED). Add that published VHH to the local novelty-comparison set. Note that `abguard`'s novelty set for the *optimisation* track is 525,063 unique CDRH3 sequences; the de novo rule is the 70%-identity-to-any-known-antibody one, which is the one that binds us. "Path to value" is 20% of the pitch, so the honest answer here is "we checked identity against the one published epitope-matched VHH and we have not run a patent search, because credits were unclaimed" — not a claim of clearance.

### IDEA-008 — Mine PMID:41193805 for the de novo antibody recipe beyond its abstract

- **status:** open
- **priority:** high
- **tests:** TH-009; TH-011
- **cost:** one AMASS/fulltext fetch, currently unaffordable — or a manual read of PMC12727541
- **source:** literature
- **new_in:** this campaign
- **evidence:** PMID:41193805
- **query_terms:** RFdiffusion de novo antibody VHH CDR length hotspot designs tested success rate; de novo antibody yeast display screening hit rate epitope specific
- **notes:** **This is the single biggest gap in the brief and it is a cache limitation, not a judgement call.** The brief asked us to mine PMID:41193805 for CDR length ranges, hotspot counts, designs per target, filters and success rates. **The cache cannot answer that: biomedcore records are abstract-only.** Every field on the record is metadata or `abstract`; there is a `hasFulltext: True` flag but no fulltext body. So from the cache we get only: VHHs, scFvs and full antibodies designed de novo against user-specified epitopes; **4** disease-relevant epitopes experimentally characterised for VHHs; cryo-EM confirmation of binding pose for influenza HA and C. difficile TcdB; initial designs at **tens to hundreds of nanomolar KD**, improved to **single-digit nM** by OrthoRep affinity maturation; all **six** CDR loops atomically accurate for one TcdB scFv. **What that already settles:** expect our first-pass designs to be weak (tens-to-hundreds of nM), because the Nature paper's were, and do not promise sub-nM from generation alone. **What it does not settle:** their CDR length ranges, hotspot counts and per-target design counts, which is why TH-009 and TH-011 lean on the Apheris docs and PMID:19062174 instead. The paper is open access as **PMC12727541** (doi 10.1038/s41586-025-09721-5) — read the methods manually if anyone has 20 minutes, and supersede TH-009's length window if their numbers disagree.

### IDEA-009 — Implement the sequence-liability filter as a script, with every threshold labelled cited or judgement in the code

- **status:** open
- **priority:** high
- **tests:** TH-010; TH-003
- **cost:** no GPU; Biopython is already installed in `/home/hackathon/nick-work/nmenv/`
- **source:** literature
- **new_in:** this campaign
- **evidence:** PMID:24492299; `/home/hackathon/hackathon-team-3/apheris/challenge/ai-antibodies-hackathon-apheris.md`
- **notes:** All of TH-010's checks are computable from sequence alone: free/unpaired Cys, N-X-S/T sequon (X != P) in CDRs, NG/NS/DG counts, Arg fraction, longest hydrophobic run, total length. Run it **before** spending refold GPU, since it is free and kills designs that `abguard` would reject anyway. The code must carry the cited-vs-judgement label next to each constant, so the pitch can defend every number or admit it is ours. Ordering matters: hard-kill checks first (they are cheap and final), soft re-rank second.

### IDEA-010 — Measure our own refold pass rate on the first 100-design batch and correct the per-arm counts

- **status:** open
- **priority:** medium
- **tests:** TH-011
- **cost:** none beyond the first arm, which we are running anyway
- **source:** literature
- **new_in:** this campaign
- **evidence:** `/home/hackathon/hackathon-team-3/apheris/rfd4/16-refolding.md`
- **notes:** The 7.7-10.6% figure is upstream's, on other targets, with cut-offs calibrated on AlphaFold3. Our own rate is the only one that should set the second arm's size. Caveat to design around: for the supplied HER2 example the automated `pass` column is documented as permanently empty, because the spec keeps HER2 and the framework in relative rather than absolute coordinates, so the scorer cannot measure how far the antibody moved. If our self-authored VHH spec inherits that, we must define our own pass rule on min_ipae and say it is ours (OUR JUDGEMENT), since **no doc anywhere gives a numeric min_ipae, iptm or binder_plddt cutoff**.
