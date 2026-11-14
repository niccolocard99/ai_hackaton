# De novo VHH on the trastuzumab epitope: reasoning

Goal: new VHH binders to HER2 domain IV that sit on trastuzumab's epitope, so they compete with
trastuzumab in the epitope-preservation assay. This is the de novo brief in the task description
("a new antibody that keeps the original epitope"). Everything below is either measured from a
structure in this repo (and says so) or marked as an assumption.

**Status: nothing has been run on a GPU.** The machine this was built on has no Docker or GPU. The
specs are validated structurally by `check_spec.py` only. Whether RFD4 accepts them, and what it
produces, is the first thing to learn on the VM (see STATE.md, Next).

## 0. Which structure, and why not 5MY6

The model is given **PDB 1N8Z, cropped by the organisers to HER2 domain IV** (`her2_fv_1n8z.cif`):
chain C 489–580 and chain D 591–607, 109 residues, in 1N8Z coordinates. Trastuzumab is not in the
model's input; it is used offline to choose and check the hotspots.

5MY6 was proposed as the VHH–HER2 structure to take hotspots from, and was checked first
(`compare_5my6.py`, output in `compare_5my6.txt`). It is HER2 bound to nanobody **2Rs15d** (2.25 Å,
Clin Cancer Res 2017). Its epitope is nowhere near trastuzumab's:

| Measure | Result |
|---|---|
| 2Rs15d epitope (HER2 residues within 4.5 Å) | 18 residues, mature numbering 87–195 (domain I) |
| Trastuzumab epitope in 1N8Z | 19 residues, 557–603 (domain IV) |
| Shared residues | **0** |
| Residues of trastuzumab's epitope even modelled in 5MY6 | 5 of 19 (its HER2 stops at mature 575) |
| Closest 2Rs15d atom to trastuzumab's Fv, after superposing on domain IV (80 Cα, RMSD 0.73 Å) | **43.7 Å** |
| Distance between the two epitope centroids in that frame | **85.8 Å** |

(5MY6 numbers HER2 22 higher than 1N8Z; 80 of 80 overlapping domain IV residues are identical at that
offset.) Hotspots taken from 5MY6 would put a VHH on domain I, where it cannot compete with trastuzumab,
so it would score nothing on the epitope property.

## 1. The epitope is too wide for one VHH, so we target a patch

Measured from `her2_fv_1n8z.cif` (HER2 residues within 4.5 Å of any trastuzumab heavy atom, `epitope_patch.py`,
output in `epitope_patch.txt`):

- 19 residues, spanning 32.4 Å between the two farthest residue centroids. A VHH paratope spans
  roughly 18–22 Å, so one VHH cannot cover the whole footprint.
- "Contact weight" is the number of distinct trastuzumab residues each epitope residue touches
  (42 in total). A patch's coverage is the share of that weight inside a sphere around it.

| Patch radius | Centre | Epitope residues covered | Contact weight covered | Heavy / light chain |
|---|---|---|---|---|
| 9 Å | F573 | 7 / 19 | 26 / 42 (62%) | 74% / 47% |
| 11 Å | D570 | 13 / 19 | 30 / 42 (71%) | 65% / 79% |
| 15 Å | P571 | 17 / 19 | 40 / 42 (95%) | 96% / 95% |

All four hotspot residues below (P572, F573, I591, K593) lie inside the 11 Å patch around D570.
Competition needs steric overlap, not full coverage, so a VHH on this patch should displace a Fab.
Caveat inherited from `litreview-optimisation`: a contact count is geometry, not energy.

## 2. Hotspots: the organisers' four atoms, checked against the guide's definition

The RFD4-Proteina annotations reference defines a hotspot by how the model was trained: *a target-chain
atom within 4.5 Å of any binder-chain atom. The count is sampled geometric(p=0.2), capped at 10% of
eligible atoms, and 5% of the time none are given. Trained alongside a target conditioned on both
coordinates and sequence, with the binder generated.* So: a handful of target interface atoms (mean
about 5), with the target fixed in coordinates and sequence. Four atoms on a target held by `C_CRD` and
`C_SEQ` is inside that range.

The atoms are those the organisers use in `Data/12-her2/her2_cdrh3.yaml` ("the HER2 atoms nearest CDRH3
in PDB 1N8Z"). `check_spec.py` re-derives each distance from the structure:

| Atom | Nearest trastuzumab CDRH3 residue | Distance |
|---|---|---|
| D / K593 NZ | G103 | 2.42 Å |
| C / F573 CD1 | Y105 | 3.44 Å |
| D / I591 CG2 | G103 | 3.74 Å |
| C / P572 CD | Y105 | 4.00 Å |

- All four are side-chain atoms with full occupancy, and the largest separation between any two is
  12.4 Å, inside one paratope. A VHH's CDR3 is the loop most likely to reach them.
- A parallel check pasted by Leo found all four present with complete side chains and reported them in the
  VH-reachable set. I re-verified presence and the contact distances above; I did not re-derive
  "VH-reachable".
- Dropped from the first draft: E558 OE1 and D560 OD1 (the acidic pair that trastuzumab's CDR-H2 binds
  with R50 and R59, not its CDR-H3). They can be added back as a second design if a run suggests the
  CDR3-centred set is too narrow.
- `C_HOT` is a trained conditioning signal, not a geometric constraint. Nothing in the spec stops a design
  binding next to the patch, so every generated complex must be checked for contacts at these atoms.
- The compactness bound in `check_spec.py` was first 12 Å, a number I chose without a source. The set
  measures 12.4 Å, so I replaced it with the paratope span (20 Å), which is the principled limit.

## 3. Framework: caplacizumab's, from the organisers' own VHH example

`Data/03-vhh/pmvk_caplacizumab.cif.zst` is the only VHH structure the organisers supply. It is a
humanised VHH framework, so it carries the hallmark residues that let it fold without a VL partner,
and the Cys22–Cys75 disulfide. We keep its backbone (`C_CRDGRP`) and let RFD4 place it rigidly.

CDR lengths are the ones the organisers encode in their example's placeholder residues:
CDR1 13, CDR2 9–10, CDR3 12–21. They are not tuned to this epitope. Total length 117–127 aa, well
inside the 250 aa limit.

## 4. The 70% identity line decides which variant is worth running

The de novo track needs under 70% identity to the nearest known antibody (VHH vs known VHHs and
heavy chains). With the framework fully fixed, the identity floor to caplacizumab itself is
fixed residues / total length, before any chance matches in the CDRs:

| Variant | Fixed framework | Floor, longest CDRs | Floor, shortest CDRs |
|---|---|---|---|
| `fixedfw` | 83 of 83 | 65.4% | 70.9% |
| `freefw` (18 positions redesigned) | 65 of 83 | 51.2% | 55.6% |

So `fixedfw` sits on the line. The floor is above 70% for short CDR3s, and chance matches in the CDRs
push it up further. It could be assigned to the optimisation track by identity alone.
`freefw` redesigns the 18 most exposed framework positions (rule in `build_inputs.py`: at most
12 CB neighbours within 10 Å, never Cys/Gly/Pro, not within 3 residues of a CDR, not the termini or
the J-region). The smoke test measured recovery of the source residue at 83% (15 of 18 freed positions
unchanged in each design), not the 30% first assumed. That puts the framework-only identity at about 65%
for a 123-aa design, against 67.5% for `fixedfw` (83/123): a small margin (section 5, item 2).

Why not just run `freefw`: a model tuned for CDRs redesigning framework surface is the less
documented path, and thermostability is one of the three scored properties. We run both at small N
and let `abguard` decide, rather than pick from this arithmetic.

What this arithmetic cannot settle: the local track checker compares only against trastuzumab's
variable domains, while the official check uses Adaptyv's antibody database, which may hold a closer
VHH than caplacizumab. A local PASS is necessary, not sufficient (AGENTS.md).

## 5. Assumptions that could be wrong

1. **The spec is untested.** It follows the organisers' Fv spec (`Data/03-fv/design.json`) closely
   and every selector resolves, but there is no organiser VHH spec to copy, so this is the main
   first-run risk. Smoke test with N=2 before anything else.
2. **`freefw` assumes omitting `C_SEQ` redesigns a position. Smoke result: it does, weakly.** In the two
   smoke designs all 65 fixed positions were unchanged, but only 3 of the 18 freed positions changed
   (E44K, R56N, A67P/V, the same three in both): 83% recovery of the source residue, not the 30% assumed
   in section 4. The abguard-based checker put identity at 68% and 66% for `freefw` against 66% and 65%
   for `fixedfw`, so on this evidence `freefw` buys no margin and carries the framework-stability risk.
   Whether the recovery comes from the side-chain coordinates in the input or from the model's own
   preference for native residues is not separable from this data. Two designs per variant is a small sample.
3. **HER2 is held with `C_CRD` (absolute), the upstream how-to prefers `C_CRDGRP`.** I followed the
   organisers' Fv example and the hotspot training description ("target conditioned on both coordinates
   and sequence"). The upstream RFD4-Proteina how-to says that for protein-protein binders it prefers
   backbone-only `C_CRDGRP` with `C_SEQ` on the target, and no `C_CRD`, because that reduces target-like
   outputs. **Measured on the smoke designs: the output does not keep the 1N8Z frame** (HER2 CA RMSD
   about 150 A to the input), but HER2 stays rigid (0.05 A after superposition). `epitope_check.py`
   therefore superposes HER2 onto 1N8Z before comparing with trastuzumab, which is valid. If designs come
   back target-like or badly docked, try one `C_CRDGRP` group spanning both HER2 chains.
4. **Two-chain target** (C and D, split at the unresolved 581–590 loop, as the organisers do) is
   untested for a design from scratch.
5. **Hit rate.** The antibody docs quote upstream binder pass rates around 8–11%, so plan for
   100–200 designs per variant, not a handful.
6. **Pre-flight.** The upstream FAQ recommends materialising the condition CIFs
   (`scripts/write_condition_spec_cifs.py`, from the upstream clone, not shipped in the images) or
   `prevalidate_inputs=true`. Whether the Apheris images expose either is unknown, so the N=2 smoke test
   is our check.

## 6. What a result means

Passing the refold score and base validation says a design is plausible and eligible. It says nothing
about SPR affinity versus trastuzumab, competition, or thermostability. An S-tier VHH (better than
trastuzumab on all three) is a long shot; an A-tier one (epitope kept plus one more) is a more
realistic target. With only 4 guaranteed tests, rank by conviction, not volume.

Fold confidence is a weak guide: `NICCOLO/eda/AGENT_GUIDE.md` section 6 measured Boltz-2 and Protenix
ipTM, `min_ipae` and pLDDT at chance for separating binders from non-binders (60 to 180 trastuzumab
CDRH3 variants). That was the optimisation setting and there is no labelled de novo data, so it may not
transfer, but it is the only measurement we have. Use refolding to check the designed pose and the
hotspot contacts, not to rank.

## 7. Run protocol

1. `python3 build_inputs.py && python3 check_spec.py spec/design_*.json` (no GPU needed).
2. On the VM, `./run_design.sh fixedfw 2 1`, then `freefw 2 1`. The script runs `check_spec.py`, the
   generation-ledger `check-run` / `add-run` / `set-status` (mandatory, AGENTS.md), and refuses to start
   while another `oslo-hackathon` container holds the shared GPU. Read the logs and outputs.
3. **Base validation (mandatory, AGENTS.md)** on every design before ranking: build the table
   (`sequence_id, kind=vhh, heavy_sequence, light_sequence, design_description`; the VHH goes whole in
   `heavy_sequence`), then run `NICCOLO/scripts/check_denovo.py` (moved there from `scripts/`) and
   `NICCOLO/scripts/check_optimisation.py`. `abguard` and `immunum` are on `PATH` on the VM. A design is
   eligible only on a `PASS` in the script for its assigned track. Then log it with the ledger (`add`),
   use `check --new-only` before refolding or ranking, and `add` again with the `base_validation`,
   `track_assigned` and `cdrh3` columns. The table builder is not written yet because the per-design CSV
   schema is only known from a real run.
4. Scale the surviving variant(s): `./run_design.sh <variant> 100 <seed>`.
5. `./refold.sh <run> 20 boltz2`, then a second model, as a sanity check on pose and hotspot contacts
   (section 6), and as provenance for the designs we submit (`NICCOLO/PIPELINE.md`, stage 6).
6. Shortlist.
   - **Thermostability:** NanoMelt is installed on the VM (`NICCOLO/thermo/`, `NICCOLO/PIPELINE.md` stage 4;
     held-out Pearson 0.853, MAE 4.1 °C on 639 measured nanobody Tm). It is trained on nanobodies, so it
     applies to VHHs; the team has so far used it only on trastuzumab VH variants. Use it to compare
     `fixedfw` with `freefw` and to shortlist, and treat differences under its 4.1 °C MAE as noise.
   - **Epitope preservation:** the team pipeline has no screen for it (`PIPELINE.md` stage 5: not built). The
     frame-based check against trastuzumab's Fv (section 5, item 3) would be one; it is not written.
   - No developability scorer (Aggrescan3D, SoluProt, NetSolP, ProteinMPNN) is installed on the VM,
     according to a teammate's note on 14 Nov.

## 8. Results so far

**Smoke test (N=2 per variant, seed 1, H100, 14 Nov).** The spec was accepted by RFD4-Proteina and every
condition resolved as designed (192 sequence residues, 436 coordinate atoms, 4 hotspot atoms, 332
internal-coordinate atoms). All four designs have `hotspot_obedience` 1.0, `sequence_adherence` 1.0, zero
clashes, and 3 to 6 interchain H-bonds. All four pass `NICCOLO/scripts/check_denovo.py` with
`assigned=de_novo`; the identity it reports is 65% and 66% for `fixedfw`, 68% and 66% for `freefw`. That is
2 to 5 points under the 70% line, so some of a larger batch will cross it and be filtered out.

`epitope_check.py` on the four designs (HER2 superposed onto 1N8Z; trastuzumab's epitope is 19 residues):

| Design | Hotspots within 4.5 Å | VHH contacts on HER2 | Share of epitope touched | Contacts inside epitope | VHH/Fv atom pairs under 3 Å |
|---|---|---|---|---|---|
| fixedfw_0 | 4 of 4 | 11 | 53% | 91% | 2969 |
| fixedfw_1 | 3 of 4 (D593 at 5.5 Å) | 12 | 42% | 67% | 1926 |
| freefw_0 | 4 of 4 | 13 | 63% | 92% | 5042 |
| freefw_1 | 4 of 4 | 12 | 58% | 92% | 2606 |

So the designs do land on trastuzumab's epitope, and each overlaps trastuzumab's Fv heavily, which is the
expected signature of competing for the same site. That is a geometric proxy for the competition assay,
not a measurement of it. Four designs say nothing about the rate across a batch.
<<<<<<< HEAD

**Batch 1 (`fixedfw`, 100 designs, seed 2) plus the 4 smoke designs, scored 14 Nov.**

- **Base validation:** 99 of 104 pass the de novo track. The 5 failures all have a 12-residue CDR3 and 70 to
  71% identity to trastuzumab, were assigned to the optimisation track, and are not reclassified.
- **CDR3 length and track risk:** every design with a CDR3 of 15 residues or fewer has an estimated framework
  identity of 70% or more to its caplacizumab scaffold (`score_designs.py`, arithmetic, not abguard's number),
  so only 17-residue and longer loops have a comfortable margin.
- **RFD4's own hotspot metric:** among the 99 passing designs, 72 have `hotspot_obedience` 1.0, 22 have 0.75
  and 5 have 0.5.
- **Epitope (the axis that carries signal):** 8 of 99 reach the team's bar (`epitope_cov` >= 0.6 and
  `contacts_in_ep` >= 0.7). Mean coverage is flat, 0.51 to 0.53, across CDR3-length buckets, so length is not
  what drives it.
- **Epitope-first ranking** puts `fixedfw_48`, `fixedfw_89`, `freefw_0` and `fixedfw_42` first, the same four
  that Nick's scoring found (`NICCOLO/SCORED_DENOVO_VHH.md`, EXP-011).
- **Boltz-2 refold of the 12 (kept, none further):** the VHH fold is confident (median pLDDT 93) but only 2 of
  12 reproduce the designed pose within 5 A (`fixedfw_6` 2.2 A, `freefw_0` 1.7 A, both ipTM about 0.96); the
  other ten land 7 to 49 A away. EXP-011 shows why this cannot be read as endorsement: Boltz-2 gives ipTM >=
  0.7 to 94% of known non-binding VHH-antigen pairs (real-versus-shuffled AUROC 0.64, against AF3 0.75).
- **Thermostability:** NanoMelt predicted Tm for the 12 is 65.6 to 73.5 C (Nick's run; mine agrees, for example
  69.7 C for `fixedfw_89`), against a 4.1 C error. It has no baseline for a VHH against trastuzumab, so it is
  a tie-break among our own designs and must not be presented as a thermostability claim.

## 9. Round 2: what changed after batch 1, and why

Information that arrived after batch 1 (`origin/nick`, PR #17, not yet on main: EXP-010, EXP-011, IDEAs 018 to
020) and what it changed:

| Finding | Change |
|---|---|
| Epitope coverage is the only axis with a usable signal; binding is unpredicted; thermostability has no baseline | Ranking is epitope-first (`score_designs.py`); the team's thresholds are the definition of a strong epitope |
| Only 8% of valid designs reach that bar | Refine the designs that already reach it (below) instead of re-rolling at 8% |
| Boltz-2 cannot separate binders from non-binders (and the team decided on 14 Nov: no Boltz) | Boltz option removed from `refold.sh`. Leo's directive is RFD4-Proteina only, so the Protenix v2 refold that IDEA-019 recommends is not run either |
| Short CDR3s carry HIGH estimated track risk at the same epitope quality | The fresh batch constrains CDR3 to 17 to 21 residues |

`queue_round2.sh` runs two steps, RFD4-Proteina only, strictly one GPU job at a time, each waiting for the card:

1. **Refine** the 5 parents that are valid, reach the epitope bar and are not HIGH track risk (`fixedfw_48`,
   `fixedfw_89`, `freefw_0`, `fixedfw_42`, `fixedfw_18`): RFD4 partial diffusion, fraction 0.4, 16 children
   each (80 designs), antibody weights, hotspot resampling off (`make_refine_spec.py`, `PARTIAL=0.4
   ./run_design.sh refine_strong 16 4`). Partial diffusion keeps the docked pose, so the epitope stays fixed
   by construction: the parent's annotations (HER2 fixed, the four hotspots, the framework) travel with its CIF.
2. **Fresh batch**, `fixedfw_long` (CDR3 17 to 21), 100 designs, seed 3. Same hotspots, same framework.
Not queued: the Protenix v2 refold of the 12 batch-1 designs that EXP-011 recommends. No model other than
RFD4-Proteina runs. `refold.sh` and `refold_each.py` stay in the folder, unused.

What is unchanged on purpose: the epitope and the four hotspots, the framework, the antibody weights.

How to read the result: the number that matters is the share of children (and of the fresh batch) that reach
the strong-epitope bar, against batch 1's 8%. A child that scores well is a design-pose measurement, not
evidence of binding. Untested: a multi-parent refine spec (`input_as_annotated_cif` with a partial-trajectory
fraction); if RFD4 rejects it, fall back to one parent per run with `--input parent.cif`.
=======
>>>>>>> origin/main
