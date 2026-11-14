# Rule verification against the primary source

[AGENTS.md](../../../AGENTS.md) instructs agents to verify rules against the challenge documents
rather than trust any summary, including itself. Checked against
`apheris/challenge/ai-antibodies-hackathon-apheris.md` and
`apheris/challenge/2026-11_Oslo_hackathon_Apheris_challenge.md`.

## 1. The `≤70%` vs `<70%` ambiguity does not affect us

The two documents disagree:

| Source | Wording |
|---|---|
| Deck, line 82 | "CDRH3 identity to trastuzumab **≤70%**" |
| Task description, line 246 | "**under** 70% CDRH3 sequence identity to trastuzumab" |
| Flowchart, line 49 | "CDRH3 identity to trastuzumab **under** 70%" |

**It is moot for an 11-residue CDRH3.** Identity after k substitutions is (11−k)/11:

| k | Identity | <70%? | ≤70%? |
|---|---|---|---|
| 2 | 81.8% | no | no |
| 3 | **72.7%** | no | no |
| 4 | **63.6%** | **yes** | **yes** |

No substitution count lands on exactly 70%, so **4 substitutions is the minimum under either
reading.** The [structural finding](06-structural-paratope-measured.md) stands unchanged.

## 2. The mutation-count off-by-one is confirmed, and `3` is the safe reading

| Source | Wording |
|---|---|
| Task description, line 36 | "**at least 2** mutations from any provided antibody" |
| Task description, line 246 | "**strictly more than 2** point mutations from any data point provided" |
| Flowchart, line 49 | "**More than 2** point mutations from any provided antibody" |

Two of three sources say *more than* 2. **AGENTS.md's advice to use 3 as the minimum is correct and
now verified against the primary source.**

## 3. Track assignment confirmed — we are in the optimisation track

The flowchart assigns the track on "**Identity to nearest known antibody** under 70%?" — not identity
to trastuzumab. A design that keeps trastuzumab's framework, light chain and five of six CDRs will be
~95%+ identical to trastuzumab, which is a known antibody. **So it routes to optimisation, as
intended**, and then faces the two extra checks: CDRH3 identity <70% ✓ and >2 point mutations from
any provided data point.

## 4. ⚠ A RISK MY OWN FINDING CREATES — read this before designing

The provided data includes (task description, line 95):

> **524,346 variants, up to 10 mutated positions in CDRH3**

The optimisation rule requires **more than 2 point mutations from *any* provided data point**. That is
a Hamming-distance-≥3 requirement against every one of those variants, measured over the whole VH+VL.

**My "mutate only the 6 contact-free positions" recommendation shrinks the design space, and a
smaller space is more likely to collide with a 524,346-variant library that already covers up to 10
CDRH3 positions.**

Space available if we restrict to the 6 robust contact-free positions (`G100, G101, A106, M107, D108,
Y109`):

| Substitutions | Position sets | Sequences |
|---|---|---|
| exactly 4 | C(6,4) = 15 | 15 × 19⁴ ≈ **1.95 million** |
| exactly 5 | C(6,5) = 6 | 6 × 19⁵ ≈ 14.9 million |
| all 6 | 1 | 19⁶ ≈ 47 million |

So the space is large in the abstract (~64 M), **but we will not sample it randomly** — we will pick
chemically sensible substitutions, and so did whoever built a FACS library of HER2 binders. The
plausible, binding-compatible corner of that space is exactly where the library is dense. **A design
that is 1 or 2 mutations from some library variant is rejected outright, not reclassified.**

### Consequences

1. **The novelty check is not optional and cannot be done on this machine.**
   `apheris/hackathon-data/her2-trastuzumab/` is absent here (git-ignored, lives on the VM). Every
   candidate must be run through `scripts/check_optimisation.py` with
   `--provided trastuzumab_her2_spr_all_768.csv` **and** against the full CDRH3 novelty set on the VM
   before it is ranked. A locally "clean" candidate means nothing.
2. **Prefer substitutions the library is less likely to contain.** A FACS screen selected for
   *binders*; the library's dense regions are the conservative, obviously-binding substitutions.
   Mutating the contact-free positions is in our favour here — those positions matter less for
   binding, so a binder-selected library had less reason to explore them systematically. *This is an
   argument, not a measurement: it needs checking against the actual tables.*
3. **Over-shoot the thresholds.** Aim for **5 substitutions, not 4**, and for Hamming distance ≥4
   from the nearest provided variant rather than exactly 3. Both thresholds are cliff edges with
   rejection on the far side, and the cost of one extra substitution at a contact-free position is low.
4. **Generate more candidates than we submit.** With a hard novelty filter we cannot evaluate locally,
   some candidates will die at validation. Build a surplus.

## 5. Also noted

- **The HER2 weights only redesign CDRH3** (line 194). To change other CDRs we would need
  `antibody-2026-10-30`. Since our strategy deliberately confines change to CDRH3, the HER2 weights
  are the right tool.
- **"A known binder fails the novelty filter, so filter designs before ranking them"** (line 194) —
  the documents say this explicitly, matching AGENTS.md's warning that the conservative weights' best
  scorers are exactly the ones that get rejected.
- The **524,346** figure here is the FACS variant count; the novelty set is **525,063** CDRH3s
  (AGENTS.md). Two different numbers for two different purposes — do not conflate them.
