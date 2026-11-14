# Ideas — the backlog

Things we have not tested, so none are lost and none are tried twice. An idea
becomes an experiment by getting an `exp_id`; it dies by getting
`status: dropped` **with a reason**, which is itself a result worth keeping.

`status`: `open` | `lit-checked` | `running` | `done` | `dropped`.
`source`: `literature` (a research pass raised it) | `experiment` (a result
suggested it) | `team`.

---

### IDEA-001 — Split the design budget: CDRH3 diversity for binding, framework mutations for Tm, in the same molecule

- **status:** open
- **priority:** high
- **tests:** TH-005; TH-002
- **cost:** one design run on `antibody-2026-10-30` (framework-capable) plus one on `her2-novel-2026-11-09`, then refold both
- **source:** literature
- **evidence:** PMID:40653772
- **notes:** The HER2 images only redesign CDRH3, so the framework half needs `antibody-2026-10-30`. This is the only idea here that attacks two scored properties at once, which is what S-tier requires.
  **EXP-005 turned this from a hunch into the measured strategy.** CDRH3 is thermostability-**inert**: 41 variants of the identical VH differing only in CDRH3 span 4.63 °C with 0 of 40 beyond the 4.1 °C MAE, while framework point mutations in the same VH do move Tm (p = 0.002). So the split is not a compromise, it is the only way to touch Tm at all. Caveat from TH-017: the framework gains available are small (best single +3.74 °C, inside the error bar; stacking sub-additive) and candidates must avoid the →Glu charge bias, unpaired cysteines, and position 98 at the CDRH3 junction.

### IDEA-002 — Cap arginine content in designed CDRH3 and check whether it costs predicted affinity

- **status:** open
- **priority:** medium
- **tests:** TH-003
- **cost:** post-hoc filter over an existing design set, no GPU
- **source:** literature
- **evidence:** PMID:28778924
- **notes:** Cheap — it is a filter over designs we already have. If high-Arg designs dominate the top of the min_ipae ranking, that is a warning that the ranking is rewarding non-specific binding.

### IDEA-003 — Choose the thermostability predictor and commit to it

- **status:** done
- **priority:** high
- **tests:** TH-016
- **cost:** one research pass plus setup; no GPU
- **source:** team
- **evidence:** PMID:36683173; PMID:39154093; PMID:41564239
- **query_terms:** scFv Tm predictor validation benchmark antibody fragment; nanobody thermostability deep learning screening criterion
- **notes:** Named as an open decision in `AGENTS.md`. Two of the three scored properties have no built-in proxy and this is the cheaper one to fix. Blocking for any S-tier claim. **pass-002 split this in two:** for a **VHH** submission use TEMPRO (PMID:39154093, MAE 4.03 °C) cross-checked against the `Tm > 65 °C` criterion in PMID:41564239 — both sdAb-specific and both reporting on Tm directly. For an **scFv** there is still no Tm-trained predictor; the best available is the TS50-trained pair in PMID:36683173, which recovered 18 of 20 known thermostable positions on a held-out mAb. Recommendation: if we submit a VHH, this decision is effectively closed; if scFv only, flag the predictor as a known weakness in the pitch rather than overclaiming.
  **DONE by EXP-005**, and the recommendation held up. The answer: **NanoMelt** (MAE 4.1 °C on 639 measured nanobody Tm), installed on the VM and verified against its own reference output. Neither platform ships anything else — Nebius `stab-ddg` is StaB-ddG (*binding* ddG), Foundry AQ is small-molecule — and cheap features reach only Pearson +0.40. It is nanobody-trained, so it is in-domain for a VHH and out-of-domain for our scFv/VH, exactly the split pass-002 predicted. See TH-016.

### IDEA-004 — Extract per-residue energetic contributions at the 1N8Z interface, not just the contact list

- **status:** done
- **priority:** medium
- **tests:** TH-009; TH-010
- **cost:** structure analysis; `bsa_1n8z.py` in the litreview folder has started on buried surface area
- **source:** literature
- **evidence:** PMID:12610629; EXP-002
- **notes:** A contact is not a hot spot. Without this, "preserves the epitope" is an assertion. The sibling task already computes BSA, so fold that in rather than redoing it. **Done by EXP-002** (sibling branch, commits 9b45cdb and ae2cdef): per-position contact geometry across 8 structures. It retired TH-008 and produced TH-009 and TH-010. Energetics are still missing — a contact frequency is not a hot spot — so the remaining half is IDEA-010.

### IDEA-005 — Submit a VHH alongside the scFv designs to hedge the track split

- **status:** open
- **priority:** high
- **tests:** TH-006; TH-012
- **cost:** one design run plus refold
- **source:** literature
- **evidence:** PMID:37632647; PMID:39154093; PMID:40814020
- **query_terms:** anti-HER2 nanobody VHH design affinity epitope; single domain antibody cell-free expression yield
- **notes:** The two tracks have separate prize pots and a team may submit any mix. A VHH is short, so it is comfortably inside the 250 aa limit and has no VH/VL pairing to get wrong. **Promoted to high by pass-004 (TH-012):** a VHH is the only format where we can both predict thermostability (TEMPRO, sdAb-only, MAE 4.03 °C) and claim better structure-prediction accuracy (AF3 2.2 Å vs 2.9 Å). Of the three scored properties, that covers the two with no built-in proxy. Counter-consideration recorded in TH-012: a single domain forgoes the light-chain CDRs that carry 45% of trastuzumab's interface contacts.

### IDEA-006 — Establish the seed-to-seed error bar on min_ipae before ranking anything on it

- **status:** open
- **priority:** medium
- **tests:** TH-001; TH-011
- **cost:** two refold runs on the shared H100
- **source:** experiment
- **evidence:** EXP-001
- **query_terms:** folding model seed variance interface score reproducibility; antibody antigen co-folding template confidence metric
- **notes:** EXP-001 ranked on a single seed and refolded only 2 of 4 designs. The gap was 0.021 min_ipae (0.310 vs 0.331) — inside the noise one seed can tell us about. **pass-003 raised the stakes** (TH-011: predicted H3 conformation for an 11-residue loop is good to only ~2 Å), and **pass-004 rescoped this idea**: PMID:40814020 needed **twenty** seeds to establish AlphaFold3's median CDR-H3 accuracy. Two seeds would not give an error bar either. So this is no longer "add a second seed" — it is: refold one design across as many Boltz-2 seeds as the shared H100 allows (5 at minimum, more if cheap), report the **spread** of min_ipae, and only then decide whether any gap between designs is real. Rank order is the submission output, so a metric without an error bar is the single largest hole in the pipeline. Runs on one design, so it is cheap.
  **Downgraded high → medium by EXP-003**, and worth saying why: this idea was going to measure *how noisy* min_ipae is. EXP-003 answered the prior question instead — whether min_ipae discriminates at all — and the answer is barely (Boltz-2 ipTM AUROC 0.645, 94.3% of shuffled decoys scoring ≥0.7). Establishing the error bar on a metric we are no longer ranking on is a lower priority than switching to one that works (IDEA-011). Keep it, because a seed spread is still needed if we use a structure score as a coarse filter, but it is no longer the top of the queue.

### IDEA-007 — Check whether any of our designs sit near a claimed patent sequence

- **status:** open
- **priority:** low
- **tests:**
- **cost:** one patentcore pass
- **source:** literature
- **evidence:**
- **query_terms:** patent: anti-HER2 antibody CDR sequence claim variant; patent: trastuzumab biosimilar CDRH3 sequence claim; patent: computationally designed antibody CDR freedom to operate
- **notes:** The digest already shows a crowded field — `Anti-HER2 antibody variants`, `High-affinity fully-humanized anti-Her2 monoclonal antibody` and ~36 more. Not a hard filter for the hackathon, but "path to value" is 20% of the pitch and freedom-to-operate is part of that answer.

### IDEA-008 — Design the minimal-risk CDRH3 variant: substitute A106, M107, D108, Y109 and nothing else

- **status:** open
- **priority:** high
- **tests:** TH-009
- **cost:** one `her2-novel-2026-11-09` design run with CDRH3 positions constrained, then refold and score
- **source:** experiment
- **evidence:** EXP-002
- **query_terms:** CDRH3 position restricted library design epitope retention; conservative CDR substitution affinity retention antibody
- **notes:** The one design in the backlog whose track eligibility is arithmetically guaranteed before any model runs: 4 substitutions on an 11-residue CDRH3 gives 63.6% identity, clearing the `<70%` gate, while sacrificing only Y109, which contacts HER2 in 1 of 8 structures. Every other contacting position is spared. Lowest constructible epitope risk, defensible position by position in the pitch.
  **EXP-003 makes the substitutions specific instead of arbitrary** (measured mean log10 KD, lower is tighter): A106 → **L** (0.64 vs wild-type A, which does not reach the top 3), D108 → **V** (1.05 vs D 1.50), M107 → **H** (0.96; note wild-type M at 0.72 is better, so this position *costs* affinity and is being spent only to reach the fourth substitution), Y109 → no data, since SPR-768 varies a 10-mer and stops at D108.
  **Open choice worth deciding deliberately:** the largest measured gain anywhere in the loop is D102 → **Q (0.15)** or **E (0.18)**, roughly a 10-fold improvement. D102 contacts HER2 in 6 of 8 structures, so TH-009 calls it risky — but a substitution that *improves* measured binding is not an epitope loss, it is a better contact. Substituting D102 instead of M107 would raise expected affinity and keep the 4-substitution count. The counter-argument is that epitope preservation is scored by competition with trastuzumab, not by affinity, so a changed contact geometry is a risk the SPR number cannot rule out. Flagging rather than deciding: this is a judgment call for the pitch, and both variants are cheap to submit.
  **REVISED by EXP-008, which replaces the guesswork with a measured pattern.** Across the 8 tightest *legal* variants in the SPR data, **M107 is mutated in 8 of 8, D102 in 8 of 8, A106 in 7 of 8**, F104 in 3 of 8 and D108 in 3 of 8. The one legal variant that beats trastuzumab is `WGGSGPYLFD` = **D102S F104P A106L M107F** at 0.095 nM against trastuzumab's 0.152 nM, and its A106L is also the best residue at that position in TH-015's table. So the recipe becomes **D102 + A106 + M107 + one of F104/D108**, not the A106/M107/D108/Y109 we first proposed: Y109 is invisible to the SPR table and D102 is mutated in every tight legal variant. Two hard constraints on using this: `WGGSGPYLFD` is in the HER2 dataset so it **fails the novelty filter** and can never be submitted — it is a template, not a candidate; and per TH-019 only 1 of 21 legal variants beat trastuzumab at all, so this design needs a second property story (Tm or epitope) to be S-tier.
  **REVISED AGAIN by EXP-009, which overturns the EXP-008 pattern on the full table.** The convergence was a base-rate artefact: M107 is mutated in **92% of all 131 legal variants** and only 75% of the 8 tightest, so "8 of 8" carried no information. With Benjamini-Hochberg correction across all ten positions, mutating **M107 (q=0.040) and F104 (q=0.028) measurably hurts** affinity — as do W99 and Y105 (both q<0.001) and G103 (q=0.090) — while **G100 (q=0.062) and G101 (q=0.077) help** and A106, D102 and D108 are inert (q>0.83). Crossed with EXP-002's contacts, **M107 never touches HER2 in 8 structures and is still the costliest position to mutate**, so non-contact means tolerated, not free.
  **The recipe becomes `G100 + G101 + A106 + D108`**: 4 substitutions, 7/11 = 63.6% identity, gate cleared. G100/G101 improve measured affinity and are the lowest-contact of the contacting positions (3/8, 2/8), so epitope risk stays low; A106/D108 never contact and are affinity-inert, so they are free filler to reach four. It avoids every position that either contacts heavily or costs affinity, and it drops the dependence on Y109, which SPR-768 cannot see — so IDEA-010's Y109 energetics question is no longer blocking. Caveats: G100/G101 sit above q=0.05, and per the revised TH-019 only 1 of 131 legal variants beat trastuzumab at all, so the second-property story is now *more* necessary, not less.

### IDEA-009 — Redesign H2 or L3 with `antibody-2026-10-30` instead of chasing CDRH3

- **status:** open
- **priority:** medium
- **tests:** TH-010
- **cost:** one design run on a weight we have not used yet, plus refold; higher risk of wasted GPU
- **source:** experiment
- **evidence:** EXP-002
- **query_terms:** CDR H2 redesign affinity antibody engineering; light chain CDR L3 contribution antigen binding affinity maturation
- **notes:** H2 carries 14-47% of interface contacts and L3 11-29%, comparable to or exceeding H3's 12-34%. All four HER2 weights redesign CDRH3 only, so the obvious tool is aimed at a minority of the interface. `antibody-2026-10-30` is the only weight that can touch other CDRs. Risk: it is not HER2-tuned, so novelty and interface quality are unproven for this target, and the `<70% CDRH3 identity` gate still has to be met by CDRH3 changes regardless. Treat as the exploratory slot, not a core submission.

### IDEA-010 — Get per-residue energetics for the CDRH3 contacts, not just contact frequency

- **status:** open
- **priority:** high
- **tests:** TH-009
- **cost:** computational alanine scan or MM-GBSA style decomposition on 1N8Z; no GPU needed for the cheap version
- **source:** experiment
- **evidence:** EXP-002
- **query_terms:** computational alanine scanning antibody antigen interface energy decomposition; MM-GBSA per residue binding free energy antibody
- **notes:** TH-009's weak point, stated in its own `evidence_against`: contact frequency across crystal structures is geometry, not energy. Y109 contacts in only 1 of 8 structures, but a rare contact can still be energetically load-bearing, and Y109 is the one position the IDEA-008 recipe sacrifices. If an energetic estimate says Y109 is a hot spot, IDEA-008 needs a different fourth substitution — and there is no safe alternative, which would make the tension in TH-009 sharper rather than looser.

### IDEA-011 — Adopt the sequence pre-filter as the primary ranker for optimisation-track designs

- **status:** open
- **priority:** high
- **tests:** TH-014; TH-013
- **cost:** none on the GPU. The code is `NICCOLO/eda/prefilter.py` on branch `origin/nick`; it needs its inputs pulled (SPR-768, FACS 524k, clonotype splits) and runs on a laptop at ~1.29M sequences/s
- **source:** experiment
- **evidence:** EXP-003
- **query_terms:** sequence based antibody affinity ranker deployment validation
- **notes:** **A guide now exists:** `NICCOLO/eda/AGENT_GUIDE.md` (commit fd1b002) documents the two models, the setup (`scp` the SPR-768 table and the FACS splits into `NICCOLO/eda/data/`), and the mistakes to avoid — read it before ranking anything. Model form: logistic L2 `C=0.3` classifier, ridge `alpha=3.0` ranker, one-hot CDRH3 over 10 positions. Out-of-fold the top 10 average 14 nM against 74 nM for a random 10. One stale path to watch: the guide points base validation at `scripts/check_*.py`, which moved to `NICCOLO/scripts/` in commit f4d8da8. The highest-value change available and it costs no GPU time. We currently rank on min_ipae, which EXP-003 shows separates real HER2 pairs from shuffled ones at AUROC 0.645 while a one-hot sequence model reaches Spearman +0.622 against measured log10 KD in the >=3-mutation band. Use the two stages as intended: classifier for "will this bind at all", ranker for "how tightly", and keep a structure score only as a coarse filter that discards the clearly bad. Coordination note: this is a teammate's code in `NICCOLO/`, which we must not write into — use it, cite it, and ask before changing it. Scope limit from TH-014: optimisation track only, since no labels exist for de novo designs.

### IDEA-012 — Run a decoy control before trusting any score we rank on

- **status:** open
- **priority:** high
- **tests:** TH-013
- **cost:** near zero: shuffle design-target pairings, refold a sample, compare score distributions
- **source:** experiment
- **evidence:** EXP-003
- **query_terms:** decoy control scoring function validation protein design benchmark
- **notes:** The transferable lesson from EXP-003, worth making a standing rule rather than a one-off. EXP-001 reported a ranking from a 0.021 min_ipae gap without anyone having established that min_ipae discriminates binders on this target; it does not. The control is to shuffle which design is paired with which target, score both sets, and check the distributions actually separate. It costs almost nothing and it is the difference between a number and a measurement. Applies to any new metric we adopt, including the thermostability predictor once chosen (IDEA-003) and the sequence ranker itself on de novo designs, where TH-014 says it has no support.

### IDEA-013 — Check ESM-IF1 against ProteinMPNN for the inverse-folding step

- **status:** open
- **priority:** medium
- **tests:**
- **cost:** low: both are sequence-design tools run on an existing backbone, no new design run needed
- **source:** literature
- **evidence:** PMID:37880344
- **query_terms:** ESM-IF1 inverse folding binder design success rate; inverse folding model comparison protein binder design benchmark
- **notes:** `AGENTS.md` names ProteinMPNN with soluble weights as a developability filter, which is a different job from designing a binder. PMID:37880344 reports ESM-IF1 producing successful binder designs toward 185 heteromeric interfaces (6.5%) against ProteinMPNN's 18 (1.5%) from 100 samples each — a ~4x difference at exactly the step where we would pick sequences for a designed backbone. **Caveat that keeps this at medium:** "successful" there means *predicted* interface RMSD <= 2 A, not measured binding, so it is one predicted metric preferring another. Per IDEA-012 this is precisely the kind of number that needs a control before it drives a decision. Worth a cheap comparison, not worth reorganising the pipeline around.

### IDEA-014 — Score every submission candidate with NanoMelt, and say which format the number is in-domain for

- **status:** open
- **priority:** high
- **tests:** TH-016; TH-012
- **cost:** negligible: ~45 sequences/s on the H100, already installed and verified on the VM
- **source:** experiment
- **evidence:** EXP-005
- **query_terms:** nanobody thermostability prediction applicability domain
- **notes:** Thermostability is one of three judged properties and the only one with no labelled data, so a predicted Tm with an honest domain caveat beats no number at all. Run NanoMelt over the final ranked list. For a **VHH** the 4.1 °C MAE applies and the ranking is defensible; for an **scFv/VH** the number is out-of-domain and must be reported as such — the ->Glu bias and `V37F` hits show the model reaching for its training domain. Pair this with the cysteine check: NanoMelt's own training data reproduces the odd/even cysteine effect (41 °C at 0, 67 °C at 2, 73 °C at 4, 53 °C at 5), which is the same rule `abguard` enforces, so a design that fails one will likely look bad in the other.

### IDEA-015 — Spend one submission slot on a framework-stabilised variant, filtered for the three known traps

- **status:** open
- **priority:** medium
- **tests:** TH-017; TH-005
- **cost:** sequence-level work plus one refold as a plausibility check; no design run needed
- **source:** experiment
- **evidence:** EXP-005
- **query_terms:** framework stabilising mutation combination epistasis antibody VH
- **notes:** 12 of 29 tested combinations clear the 4.1 °C error bar and random 6-mutation subsets reach +6.09 °C, so there is something to win — but TH-017 says it is bounded and contaminated. Any candidate must: exclude the fragile positions (W36, C22, C96, L20, M83, V48, D73, E89), exclude position 98 at the `...YYC-S-R|WGGDGFYAMD` junction since it dominates the top hits and is nearest the paratope, keep the cysteine count unchanged so `abguard` does not reject it, and not be a pure net-charge shift (11 of the top 12 singles are ->Glu, with unknown cost to expression). One slot, because the Tm gain is predicted by a model that is out-of-domain for a VH.

### IDEA-016 — Score the final list for expression and self-association, not just the three judged properties

- **status:** open
- **priority:** medium
- **tests:**
- **cost:** low: all four tools are sequence/structure-level, no GPU design run
- **source:** literature
- **evidence:** PMID:31181988; PMID:23995618
- **query_terms:** scFv self association viscosity concentration dependent developability; cell-free expression failure sequence determinants antibody fragment
- **notes:** A design that fails to express is not scored badly, it is **unmeasurable** — and Adaptyv reports expression success improved when participants used Aggrescan3D, SoluProt, NetSolP and ProteinMPNN soluble weights (root `AGENTS.md`). Nothing in this backlog covered that until now, which is a gap given the whole submission depends on the assays returning numbers at all. Also worth noting that a high Tm does not clear developability on its own: tibulizumab's anti-IL-17 scFv had Tm 73.1 °C by DSC yet showed concentration-dependent, initially reversible self-association (PMID:31181988), while CDR-restricted scFv engineering reached 100 mg/ml with <2% HMW species after 7 weeks (PMID:23995618). Cheap insurance on the four designs that actually get tested.

### IDEA-017 — Get the full SPR-768 table into the analysis, not just the 60-variant calibration subset

- **status:** done
- **priority:** high
- **tests:** TH-019; TH-018; TH-015
- **cost:** one `scp` from the VM, then re-run two scripts we already have
- **source:** experiment
- **evidence:** EXP-008; EXP-007
- **query_terms:**
- **notes:** Every affinity number in this ledger comes from the **60-variant fold-calibration subset**, the only SPR table committed to the repo. That subset was deliberately built as 30 binders spanning the KD range plus 30 non-binders, so it cannot give a base rate: TH-019's "1 of 21 legal variants beats trastuzumab" is an existence proof, not a frequency. The full table has **179 quantified** variants of 620 outcomes. Re-running `analysis/legal_band_affinity.py` against it would turn three soft claims into real ones — the legal-band hit rate against 0.152 nM, the per-position substitution frequencies behind IDEA-008's recipe, and the n backing TH-015's per-position preferences (currently 4-5 at the standout residues). The path is in `NICCOLO/eda/AGENT_GUIDE.md`: `scp` `trastuzumab_her2_spr_all_768.csv` from `hackathon-data/her2-trastuzumab/` on the VM. Cheapest high-value item on the board.
  **DONE by EXP-009** — no `scp` was needed, the full table was already on the `nick` session's laptop. It cost this ledger two claims: the legal-band hit rate fell from 1 of 21 to **1 of 131** (0 of 131 against trastuzumab's best replicate), and the per-position frequencies behind IDEA-008 turned out to be base-rate artefacts. TH-015's thin n=4-5 residue preferences are superseded at the position level by n=131.
