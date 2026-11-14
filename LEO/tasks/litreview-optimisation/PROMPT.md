# Task: literature review — trastuzumab optimisation

You have AMASS access (papers, clinical trials, patents, regulatory records,
molecules, genes). Use it as the primary source and cite specific records.

## Goal

Gather what has been published on engineering **trastuzumab** (anti-HER2 IgG1,
PDB 1N8Z), and convert it into **actionable, position-level guidance** for
designing improved variants. I am not looking for a summary of the field — I
need specific mutations, positions and mechanisms I can act on, each traceable
to a source.

Scope this pass to the **optimisation** of trastuzumab (engineering the existing
antibody). A second pass will cover de novo anti-HER2 binders; do not cover that
here, but flag anything that transfers.

## The three properties being measured

Designs are scored against measured trastuzumab baseline on exactly three
properties. Organise findings around them:

1. **Binding affinity** — measured by SPR. Report KD values and the assay
   conditions they came from, since cross-paper comparisons are often invalid.
2. **Epitope preservation** — measured by competitive binding against
   trastuzumab, i.e. the variant must still engage HER2 domain IV at the same
   site. Mutations that improve affinity by shifting the epitope are failures,
   not wins. Identify which contact residues are load-bearing for the epitope.
3. **Thermostability** — directly measured. Report Tm/Tagg values and method
   (DSF, DSC, nanoDSF) where given.

Note explicitly where these three **trade off** against each other — affinity
gains that cost stability are the main known failure mode and I need to know
where that boundary sits.

## Constraints that decide whether a finding is usable

Treat these as filters on what you report, and say when a promising finding
fails one:

- Variants must differ from trastuzumab in **CDRH3** by enough to bring CDRH3
  identity **below 70%** — several substitutions, not one or two. So classic
  single-point affinity-maturation results are informative about mechanism but
  not directly transplantable. Prioritise work on **CDRH3 loop redesign,
  tolerated loop length/conformation changes, and which CDRH3 positions tolerate
  wholesale change while preserving the HER2 domain IV epitope**.
- Designs must be **at least 3 point mutations** from any published/measured
  variant, so exact known variants cannot be reused — I need the underlying
  rules, not a list to copy.
- Total length ≤250 aa (scFv, VHH, or VH+VL combined); no unpaired cysteines.
  Flag any proposed mutation that introduces a free cysteine.
- Framework and non-CDRH3 CDR changes are relevant mainly for stability and
  expression.

## Also needed

- **Expression/developability:** which sequence features predicted
  cell-free expression failure for antibodies, and what the literature says
  about Aggrescan3D, SoluProt, NetSolP and ProteinMPNN soluble weights as
  filters.
- **Thermostability prediction:** which computational predictors are actually
  validated on antibody Fv/scFv domains. I have not chosen one yet and need a
  defensible recommendation.
- **Patents:** which trastuzumab variants and HER2-binding CDR sequences are
  already claimed, so I can see what is crowded prior art.

## Deliverable

1. **Findings table:** `property | position(s)/region | observed effect (with
   magnitude) | mechanism | source record | confidence | usable under the
   constraints above (yes/no/why)`.
2. **Ranked recommendations:** the handful of design strategies best supported
   by evidence, each with its reasoning and the evidence against it. Only ~4
   designs will be experimentally tested, so precision beats coverage —
   a short, well-argued list is worth more than a long one.
3. **Contradictions and gaps:** where the literature disagrees, and what is
   simply unknown.
4. **Sources:** every claim linked to its AMASS record. Distinguish measured
   results from computational predictions.

## Working rules

- Prefer primary measured data over reviews; name the assay behind every number.
- State uncertainty plainly. "Unknown" is a useful answer; a confident guess is not.
- Credit budget is limited — cache results, do not re-run near-identical
  searches, and batch related queries.
