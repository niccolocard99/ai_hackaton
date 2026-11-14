# Scaffold library for de novo anti-HER2 loop design

30 scaffolds: 20 VHH frameworks + 2 anti-HER2 VHH positive controls + 8 paired Fv on
non-VH3 / non-Vkappa1 germlines, plus the two RFD4-bundle references
(caplacizumab VHH, adalimumab Fv). Structures in `cif/`, one file per scaffold,
Fv-only, hydrogens stripped, renumbered 1..N, chains renamed **H** (VHH / VH) and
**L** (VL) so nothing collides with HER2 antigen chains C and D.
Per-scaffold numbers in `scaffolds.csv`.

## How the numbers were measured

* CDR/framework boundaries come from **immunum IMGT segmentation** — the same
  numbering `abguard` uses inside the gate, so the split is the gate's own.
* For every scaffold, 100 random-CDR fills were generated with the framework held
  fixed and the three CDRs re-drawn at the **shipped CDR lengths** from a CDR-like
  amino-acid composition (Cys excluded, so an odd cysteine count can only come from
  the framework). A second set of 50 fills per scaffold widens each CDR span to the
  **RFD4-bundle contig convention** (CDR1 +3/+2, CDR2 +1/+1, CDR3 unchanged;
  verified to reproduce `Data/03-fv/design.json` exactly). The `_wideCDR` columns are
  the ones that describe what RFD4 will actually emit.
* All 4530 fills were scored with the real gate,
  `/home/hackathon/.venvs/abguard/bin/python NICCOLO/scripts/check_denovo.py`.
* `id_own_parent` uses abguard's own identity formula
  (`1 - (d_H + d_L) / (max|VH| + max|VL|)`, Levenshtein on the numbered Fv) against
  the scaffold's **own parent** variable domain, which abguard does *not* check.

## The finding that should change the plan

**Identity to trastuzumab is a solved problem. Identity to the scaffold's own parent
is not, and it is the binding constraint for every kept-framework design.**

| gate | best VHH | worst VHH | best paired Fv |
|---|---|---|---|
| vs trastuzumab (<70%) | 52.0% | 65.8% | 44.2% |
| vs own parent (<70%) | 68.8% | 78.8% | 71.5% |

Every scaffold in the library clears the trastuzumab gate with 0/100 failures except
the bundle's adalimumab Fv. **No scaffold reliably clears a <70%-to-own-parent bar at
its shipped CDR lengths**: the best, `vhh_6XYF_A`, still fails 30% of the time, and
the paired Fv fail 92-100% of the time. This is structural, not scaffold-specific: keeping a natural framework means identity
to that framework's parent is roughly `1 - k * n_randomised / L_Fv`, and the measured
`k` is about 0.78, not the 0.93 a naive "random residues almost never match" argument
gives - a random draw sometimes hits the parent residue, and Levenshtein can find a
cheaper alignment than position-by-position. So clearing 70% needs **~50 randomised
residues in a ~120-residue VHH**, where the shipped contigs randomise only 34-43.

Two consequences:

1. **The lever is CDR3 length, not scaffold choice — and it is bigger than it looks.**
   Measured by sweeping CDR3 length on the five best VHHs (40 fills per point,
   1200 designs, all 1200 passing the real gate):

   | CDR3 extension | own-parent median | own-parent max | own-parent FAIL frac | trastuzumab median | trastuzumab FAIL frac |
   |---|---|---|---|---|---|
   | +0 (shipped) | 0.692-0.704 | 0.717-0.750 | 0.33-0.72 | 0.565 | 0.00 |
   | +4 | 0.682-0.695 | 0.710-0.734 | 0.12-0.45 | 0.550 | 0.00 |
   | +8 | 0.669-0.686 | 0.707-0.727 | 0.03-0.12 | 0.541 | 0.00 |
   | **+12** | 0.657-0.676 | 0.691-0.706 | **0.00-0.03** | 0.526 | 0.00 |
   | **+16** | 0.643-0.660 | 0.678-0.700 | **0.00-0.03** | 0.517 | 0.00 |
   | +20 | 0.633-0.648 | 0.653-0.676 | **0.00** | 0.507 | 0.00 |

   **CDR3 must be ~12-16 residues longer than shipped** (so ~29-35 residues total for
   these VHHs, long but well inside the camelid range) before a kept-framework design
   is reliably under 70% against its own parent. Trastuzumab identity only improves as
   CDR3 grows, so there is no trade-off between the two gates - only a trade-off
   against designability of a very long loop. `scaffolds.csv` carries
   `recommended_cdr3_len` and a ready-made `rfd4_contig_cdr3_extended` per VHH.
   For paired Fv the same calculation needs ~30 extra randomised residues across both
   chains, which is framework redesign rather than loop design.
2. **Prefer a scaffold whose parent is obscure.** If the organisers read "any known
   antibody" as "anything in OAS / the PDB / abYsis", a deposited research nanobody is
   as findable as a drug. If they read it as "any clinical or named antibody", the
   difference between an unpublished phasing scaffold and caplacizumab is everything.
   Both readings are served by the same choice, so pick the obscure parent.

## Ranked recommendation for round 1

Ranked on own-parent identity first (the binding gate), then trastuzumab identity,
then resolution and parent obscurity.

| # | scaffold | why | framework SEQUENCE kept? |
|---|---|---|---|
| 1 | `vhh_6XYF_A` (PDB 6XYF chain A, 1.11 A) | Lowest own-parent identity in the library (68.8%, 33% fail at shipped length; 0.03 at CDR3 +12). Parent is "Nanobody 22", an engineered lanthanide-binding **crystallographic phasing scaffold with no primary publication** - about as far from "a known antibody" as a real VHH gets. CDR3 19 already. Trastuzumab identity 56.8%, 0/100 fail. | **Yes, with CDR3 >= 31** (shipped 19, +12). Safest kept-framework option. |
| 2 | `vhh_8E2N_B` (8E2N chain B, 1.10 A) | Own-parent 69.3% (44-47% fail at shipped length; 0.00 at CDR3 +12). Divergent FR1 (`...GGGLVPAGGSLRLSCTTS`), unusually long CDR2 (11) *and* CDR3 (17) = 36 shipped CDR residues, the most of any VHH here. Parent VHH113 (anti-Pseudomonas Cif), no primary publication. Trastuzumab 57.5%. | **Yes, with CDR3 >= 29** (shipped 17, +12). |
| 3 | `vhh_6JB9_A` (6JB9 chain A, 1.15 A) | Longest shipped CDR3 (20) -> own-parent 70.1%. **Apo structure**: the loops are not moulded onto an antigen, which matters when you are about to replace them. Trastuzumab 54.3%, the second lowest of any VHH. | **Yes, with CDR3 >= 32** (shipped 20, +12). Parent (Nb D3-L11, anti-lysozyme) is published. |
| 4 | `vhh_8CDS_A` (8CDS chain A, 1.53 A) | **Lowest trastuzumab identity of any VHH (52.0%, max 55.2%)** - the biggest safety margin on the gate that is actually enforced. Most germline-divergent framework in the set (`QVQLVESGGGVVQTGDSLMLSCTVS`, FR2 `WFRQT/QREFVA`). Own-parent 70.4% (72% fail at shipped length, 0.00 at +12). | **Yes, with CDR3 >= 30** (shipped 18, +12). |
| 5 | `vhh_5VLV_A` (5VLV chain A, 1.35 A) | **Apo** hapten-binder framework (VHH T9), FR1 `LVQTGDSLRLSC`, CDR3 17. Own-parent 70.2%, trastuzumab 57.3%. Good structural diversity partner to #1-#4 (different FR1 and FR3). | **Yes, with CDR3 >= 29** (shipped 17, +12). |
| 6 | `fv_2XZC_HL` (2XZC chains H/L, 1.36 A) | The one paired Fv worth a round-1 slot. IGHV4-like + IGLV1-like - **two germline families away from trastuzumab's VH3/Vkappa1 on both chains**: trastuzumab identity 49.8% (max 52.0%). Parent is the A.17 "reactibody", a laboratory catalytic antibody, **not a therapeutic** - the lowest own-parent reputational risk of any Fv here even though its own-parent *identity* (73.3%) is worse than any VHH. | **Only with framework diversification.** Even CDR3H +16 leaves it near 70% against its own parent. Do not ship it with the framework sequence untouched. |

Runner-up worth keeping warm: `fv_7SD5_HL` (1.53 A) has the **lowest trastuzumab
identity in the whole library (44.2%)** and the longest CDRH3 (22), so it is the best
paired Fv on the enforced gate; its parent (mAb 10-40) is a clinical-stage
SARS-CoV-2 antibody, which is the own-parent risk.

## Safe to use with the framework SEQUENCE KEPT vs needs diversification

**Keep the framework (VHH only, and only with CDR3 at or above
`recommended_cdr3_len` in `scaffolds.csv` = shipped + 14):** `vhh_6XYF_A`,
`vhh_8E2N_B`, `vhh_6JB9_A`, `vhh_8CDS_A`, `vhh_5VLV_A`, `vhh_7AQZ_C`, `vhh_6HER_B`.
All seven sit at 68.8-71.0% own-parent at shipped CDR lengths - i.e. **on the wrong
side of the line a third to two thirds of the time** - and only drop clear of it with
a +12 to +16 CDR3 extension, while staying 12-19 points under the trastuzumab gate
throughout. **At shipped CDR lengths, none of them is safe against an own-parent
rule**; they are the scaffolds for which the fix is a longer loop rather than a new
framework.

**Needs framework diversification (do not ship with the sequence kept):**

* **All 8 paired Fv** and **both bundle references**. Own-parent identity is
  71.5-76.4% at the bundle contig spans; getting under 70% needs roughly 30 extra
  randomised residues spread over both chains, which is framework redesign, not
  loop design.
* `ref_fv_adalimumab` additionally **fails the trastuzumab gate itself**: 73.5%
  median (100/100 fail) at IMGT CDR spans, 69.3% median with 14/50 fail at the
  bundle's own wider contig spans. It is not usable as shipped.
* `ref_vhh_caplacizumab` passes the trastuzumab gate 100/100 (median 65.0%) but its
  parent is an **approved drug**, and own-parent identity is 72.5%. Use it as a
  reference/control, not as a submission framework.
* `vhh_8CDT_A` and `vhh_7UTG_A`, `vhh_4B41_A`, `vhh_7R20_B`: shortest CDR3s (11-16)
  and the highest VHH trastuzumab identities (62-66%; `vhh_8CDT_A` has a *classic
  VH3* FR2 `WVRQAPGKGLEWVS` rather than a camelid hallmark, which is exactly why).
  Usable, but with the least margin on both gates.

## Hard exclusions

* **`vhh_5IVN_A` (BC2 nanobody, 1.00 A): EXCLUDE.** 100/100 fills fail the gate's
  hard check with `unpaired_cysteine` - the framework/CDR scaffold carries extra
  cysteines (FR2 `EREGVSC`, plus one in CDR3). The only scaffold in the library with
  this problem. It is otherwise the highest-resolution structure here, so if it is
  wanted, mutate the extra Cys to Ser first and re-measure.
* **`fv_9GWZ_HL`: EXCLUDE as a framework.** Clean at shipped IMGT spans, but 50/50
  of the wide-CDR fills fail `unpaired_cysteine`: the bundle-style CDR1 widening eats
  into a cysteine-bearing position. Its parent (23ME-00610) is also a
  first-in-human clinical antibody.
* **`ctrl_vhh_9MTE_A` / `ctrl_vhh_9MTE_B`: positive controls only, never frameworks.**
  These are the anti-HER2 VHHs 27A05 and 47D05 from PDB 9MTE (2.66 A cryo-EM, HER2
  ECD in chain C of the original entry; PMID 41504346, the TPP-45142 HER2 T-cell
  engager). They are the useful thing the brief asked to note: a VHH-HER2 complex to
  check an epitope and a pipeline against. 27A05 has only 109 of its residues
  modelled (CDR3 7), so 47D05 (`ctrl_vhh_9MTE_B`, 117 aa, CDR3 13) is the better
  control. Both would be their own identity risk as frameworks (own-parent 71.8% and
  76.2%).

## Numbering confidence and cysteines

Every one of the 4530 fills numbered with immunum confidence **>= 0.50** (per-scaffold
minimum 0.50-0.68), far above the 0.1 floor the hard check applies, and no fill failed
on `numbering_failed`. Cysteine parity is even in every domain of every scaffold
except the two exclusions above (`vhh_5IVN_A`, and `fv_9GWZ_HL` at widened spans).

## Premises in the brief that the measurements corrected

1. **"Caplacizumab VHH: median 61.4% (55.9-65.1), 200/200 pass."** Direction
   confirmed, margin smaller: measured **median 65.0%, range 62.5-68.3%**, 0/100 fail.
   The maximum sits 1.7 points under the gate, not 4.9. Still safe, but it is not the
   comfortable scaffold the brief implies, and six VHHs in this library beat it by
   5-13 points.
2. **"Adalimumab Fv: median 68.5% (65.5-71.3), 16/200 FAIL -> sits ON the cliff."**
   It is worse than on the cliff: **median 69.3% with 14/50 (28%) failing** at the
   bundle's own contig spans, and **median 73.5% with 100/100 failing** when only the
   IMGT CDRs are randomised. The adalimumab parent itself is **77.2%** identical to
   trastuzumab. Unusable as shipped, as the brief concluded - just more decisively.
3. **"83 kept framework residues out of ~117-127 total" for caplacizumab.** Correct
   for the shipped file (the bundle keeps 83 real residues + 3 UNK placeholders = 86
   modelled positions, over an Fv of ~120). Worth stating precisely because the ratio
   is the whole argument: the bundle's framework/CDR split keeps **7 fewer** framework
   residues per heavy chain than IMGT does, and that 7-residue difference is what
   moves adalimumab from 73.5% to 69.3%.
4. **"VHH first; paired Fv only on non-VH3/non-Vkappa1 germlines."** Confirmed and
   strengthened: non-VH3/non-Vkappa1 paired Fv land at **44-57%** identity to
   trastuzumab, better than any VHH. But they are *worse* on the own-parent gate
   (71.5-76.4% vs 68.8-71.0%) purely because an Fv is twice as long for roughly the
   same CDR content. **If the own-parent rule is enforced, VHH wins for the opposite
   reason to the one in the brief** - not because of germline distance, but because a
   120-residue domain has a higher randomised-fraction than a 230-residue one.
5. **"Extending the loops a little is enough" (my own first estimate) is wrong,
   and it is worth recording because the arithmetic is seductive.** A naive model
   (random CDR residues differ from the parent ~93% of the time) says a 0-5 residue
   CDR3 extension clears the own-parent gate. Measured, the effective rate is ~78%
   and the required extension is **+12 to +16**, three times larger. Anyone
   re-deriving this should measure, not estimate.
6. **The identity formula is not symmetric in the way "65-75% by construction"
   suggests.** abguard pools edit distance over the chains the design *has*, so a VHH
   is scored on its heavy chain alone. That is why VHH own-parent identity (74-79% at
   IMGT spans) is *lower* than paired-Fv own-parent identity (78-82%) despite VHHs
   keeping a larger fraction of a shorter framework.

## Full table: VHH (sorted by own-parent identity, best first)

| scaffold | res A | CDR1/2/3 H | id→trast med | id→trast max | trast FAIL frac | id→own parent med | own-parent FAIL frac | CDR3 ext needed | numbering conf | cysteine |
|---|---|---|---|---|---|---|---|---|---|---|
| vhh_6XYF_A | 1.110 | 8/8/19 | 0.592 | 0.632 | 0.000 | 0.688 | 0.300 | 0 | 0.579 | ok (even per domain) |
| vhh_8E2N_B | 1.100 | 8/11/17 | 0.575 | 0.606 | 0.000 | 0.693 | 0.440 | 0 | 0.520 | ok (even per domain) |
| vhh_6JB9_A | 1.150 | 8/8/20 | 0.567 | 0.591 | 0.000 | 0.701 | 0.560 | 0 | 0.527 | ok (even per domain) |
| vhh_5VLV_A | 1.350 | 8/8/17 | 0.589 | 0.637 | 0.000 | 0.702 | 0.700 | 1 | 0.548 | ok (even per domain) |
| vhh_8CDS_A | 1.530 | 8/8/18 | 0.520 | 0.552 | 0.000 | 0.704 | 0.540 | 0 | 0.528 | ok (even per domain) |
| vhh_7AQZ_C | 1.300 | 8/8/17 | 0.557 | 0.590 | 0.000 | 0.705 | 0.640 | 0 | 0.500 | ok (even per domain) |
| vhh_6HER_B | 1.200 | 8/8/17 | 0.613 | 0.637 | 0.000 | 0.710 | 0.720 | 1 | 0.626 | ok (even per domain) |
| vhh_8CDT_A | 1.410 | 8/7/16 | 0.620 | 0.645 | 0.000 | 0.711 | 0.860 | 2 | 0.617 | ok (even per domain) |
| vhh_5J1S_C | 1.400 | 8/7/16 | 0.598 | 0.631 | 0.000 | 0.713 | 0.820 | 2 | 0.565 | ok (even per domain) |
| vhh_4B50_A | 1.300 | 8/7/16 | 0.607 | 0.631 | 0.000 | 0.713 | 0.860 | 2 | 0.572 | ok (even per domain) |
| vhh_2X1O_A | 1.340 | 8/8/15 | 0.592 | 0.633 | 0.000 | 0.717 | 0.820 | 1 | 0.593 | ok (even per domain) |
| ctrl_vhh_9MTE_B | 2.660 | 8/7/13 | 0.617 | 0.650 | 0.000 | 0.718 | 1.000 | 3 | 0.557 | ok (even per domain) |
| vhh_2P45_B | 1.100 | 8/8/14 | 0.603 | 0.628 | 0.000 | 0.719 | 0.980 | 3 | 0.565 | ok (even per domain) |
| vhh_5IVN_A | 1.000 | 8/8/15 | 0.590 | 0.631 | 0.000 | 0.721 | 0.920 | 2 | 0.570 | UNPAIRED-CYS in 100% of fills |
| vhh_6APO_A | 1.170 | 9/7/11 | 0.542 | 0.575 | 0.000 | 0.724 | 1.000 | 4 | 0.500 | ok (even per domain) |
| ref_vhh_caplacizumab |  | 8/7/14 | 0.650 | 0.683 | 0.000 | 0.725 | 1.000 | 3 | 0.598 | ok (even per domain) |
| vhh_4B41_A | 1.190 | 8/7/12 | 0.658 | 0.683 | 0.000 | 0.727 | 1.000 | 4 | 0.606 | ok (even per domain) |
| vhh_5M2W_B | 1.500 | 8/7/13 | 0.608 | 0.642 | 0.000 | 0.731 | 1.000 | 4 | 0.579 | ok (even per domain) |
| vhh_7UTG_A | 1.250 | 8/8/11 | 0.642 | 0.675 | 0.000 | 0.735 | 1.000 | 4 | 0.597 | ok (even per domain) |
| vhh_7R20_B | 1.420 | 8/7/12 | 0.633 | 0.658 | 0.000 | 0.737 | 1.000 | 5 | 0.580 | ok (even per domain) |
| ctrl_vhh_9MTE_A | 2.660 | 8/6/7 | 0.600 | 0.625 | 0.000 | 0.761 | 1.000 | 8 | 0.569 | ok (even per domain) |

## Full table: paired Fv (sorted by own-parent identity, best first)

| scaffold | res A | CDR1/2/3 H | id→trast med | id→trast max | trast FAIL frac | id→own parent med | own-parent FAIL frac | CDR3 ext needed | numbering conf | cysteine |
|---|---|---|---|---|---|---|---|---|---|---|
| fv_7SD5_HL | 1.530 | 10/7/22 | 0.454 | 0.479 | 0.000 | 0.715 | 0.920 | 5 | 0.597 | ok (even per domain) |
| fv_8F6L_AB | 1.850 | 8/7/15 | 0.567 | 0.593 | 0.000 | 0.727 | 1.000 | 10 | 0.655 | ok (even per domain) |
| fv_2XZC_HL | 1.360 | 9/7/12 | 0.498 | 0.520 | 0.000 | 0.733 | 1.000 | 10 | 0.599 | ok (even per domain) |
| fv_7VYT_HL | 1.530 | 8/8/14 | 0.550 | 0.576 | 0.000 | 0.734 | 1.000 | 10 | 0.672 | ok (even per domain) |
| fv_2VXV_HL | 1.490 | 8/8/14 | 0.546 | 0.576 | 0.000 | 0.734 | 1.000 | 11 | 0.602 | ok (even per domain) |
| fv_7SEG_HL | 2.160 | 8/8/13 | 0.485 | 0.511 | 0.000 | 0.739 | 1.000 | 12 | 0.666 | ok (even per domain) |
| ref_fv_adalimumab |  | 8/8/14 | 0.735 | 0.754 | 1.000 | 0.741 | 1.000 | 12 | 0.678 | ok (even per domain) |
| fv_6OC7_HL | 1.300 | 8/8/10 | 0.478 | 0.493 | 0.000 | 0.748 | 1.000 | 13 | 0.647 | ok (even per domain) |
| fv_9GWZ_HL | 2.120 | 8/7/12 | 0.558 | 0.580 | 0.000 | 0.764 | 1.000 | 11 | 0.614 | ok (even per domain) |

## Sequences, germlines and ready-made RFD4 contigs

CDRs are lowercase. Contigs are in the numbering of the files in `cif/` (1..N per
chain) and use the bundle's CDR-span convention, verified to reproduce
`Data/03-fv/design.json` exactly, so they are drop-in replacements for the
`contig` field of a `design.json`. Kept blocks carry input-structure numbering;
only the loop lengths change. `scaffolds.csv` also carries
`rfd4_contig_cdr3_extended` with the longer CDR3 that clears the own-parent gate.

| scaffold | source | germline | VH / VHH (CDRs lowercase) | VL (CDRs lowercase) | RFD4 contig (shipped lengths, output numbering) |
|---|---|---|---|---|---|
| fv_2VXV_HL | PDB 2VXV chains H(H)/L(L) | IGHV5-51-like (GTEVKKPGESLKISCKGS) + IGKV3-15-like (EIVMTQSPATLSVSPGERATLSC) | EVQLVQSGTEVKKPGESLKISCKGSgytvtsywIGWVRQMPGKGLEWMGFiypgdsetRYSPTFQGQVTISADKSFNTAFLQWSSLKASDTAMYYCarvgsgwypytfdiWGQGTMVTVSS | EIVMTQSPATLSVSPGERATLSCRASesissnLAWYQQKPGQAPRLFIYtasTRATDIPARFSGSGSGTEFTLTISSLQSEDFAVYYCqqynnwpsitFGQGTRLEIK | H1-22, 13-13, H36-49, 10-10, H60-96, 14-14, H111-121, /0, L1-23, 11-11, L35-48, 8-8, L57-88, 10-10, L99-108 |
| fv_2XZC_HL | PDB 2XZC chains H(H)/L(L) | IGHV4-like (ESGPGLVKPSETLSLTC) + IGLV1-like (SVLTQPPSVSAAPGQKVTISC) | VQLQESGPGLVKPSETLSLTCAVSgysissgyyWGWIRQPPGKGLEWIGSiyhsgstYYNPSLKSRVTISVDTSKNQFSLKLSSVTAADTAVYYCagltqsshndanWGQGTLTTVSS | SVLTQPPSVSAAPGQKVTISCSGSssnignnyVSWYQQLPGTAPKLLIYdnnKRPSGIPDRFSGSKSGTSATLGITGLQTGDEADYYCgtwdsslnpvFGGGTKLEI | H1-21, 14-14, H36-49, 9-9, H59-95, 12-12, H108-118, /0, L1-21, 13-13, L35-48, 8-8, L57-88, 10-10, L99-107 |
| fv_6OC7_HL | PDB 6OC7 chains H(H)/L(L) | IGHV1-46-like + IGLV3-like (YELTQPPSVSVSPGQTARITCSGD) | VQLVQSGAEVKKPGASVKVSCKASgytftgyyMHWVRQAPGQGLEWMGWinpnsggtNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCardrgyymdvWGKGTTVTVSS | YELTQPPSVSVSPGQTARITCSGDalpkqyAYWYQQKPGQAPVLVIYkdsERPSGIPERFSGSSSGTTVTLTISGVQAEDEADYYCqsadsgtyvvFGGGTKLTVL | H1-21, 13-13, H35-48, 10-10, H59-95, 10-10, H106-116, /0, L1-21, 11-11, L33-46, 8-8, L55-86, 10-10, L97-106 |
| fv_7SD5_HL | PDB 7SD5 chains H(H)/L(L) | IGHV4-like + IGLV8-like (NFMLTQPHSMSESPGKTVTISC) | QLQLQESGPGLVKPSETLSLTCTVSggsisssnfyWGWIRQPPGKGLEWIASitysgrtFYNPSLKSRVAISVDTSKNQFSLKLSSVTAADTAVYYCartfpsyydrsgyhylnygmdvWGQGTTVTVSS | NFMLTQPHSMSESPGKTVTISCTRSsgsiasnyVQWYQQRPGSSPTTVIYednQRPSGVPDRFSGSIDSSSNSASLTISGLKTEDEADYYCqsydssswvFGGGTKLTVL | H1-22, 15-15, H38-51, 9-9, H61-97, 22-22, H120-130, /0, L1-22, 13-13, L36-49, 8-8, L58-91, 9-9, L101-110 |
| fv_7SEG_HL | PDB 7SEG chains H(H)/L(L) | IGHV5-51-like (GESLKVSCKAS) + IGLV6-57-like (SYVLTQPSSVSVAPGQTATISC) | EVQLVQSGAEVKKPGESLKVSCKASgytftsyyMHWVRQAPGQGLEWMGAiepmygstSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCargsayyydfadyWGQGTLVTVSS | SYVLTQPSSVSVAPGQTATISCGGHnigsknVHWYQQRPGQSPVLVIYqdnKRPSGIPERFSGSNSGNTATLTISGTQAMDEADYYCqvwdnysvlFGGGTKLTVL | H1-22, 13-13, H36-49, 10-10, H60-96, 13-13, H110-120, /0, L1-22, 11-11, L34-47, 8-8, L56-87, 9-9, L97-106 |
| fv_7VYT_HL | PDB 7VYT chains H(H)/L(L) | IGHV1-46-like (GAEVKKPGASVKVSCKAS) + IGKV3-20-like (EIVLTQSPGTLSLSPGERATLSC) | EVQLVQSGAEVKKPGASVKVSCKASgytftsyyMHWVRQAPGQGLEWMGIinpsggstSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCasrsgsgwfgaldyWGQGTLVTVS | EIVLTQSPGTLSLSPGERATLSCRASqsvsssyLAWYQQKPGQAPRLLIYgasSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCqqygsspggtFGQGTKVEIK | H1-22, 13-13, H36-49, 10-10, H60-96, 14-14, H111-120, /0, L1-23, 12-12, L36-49, 8-8, L58-89, 10-10, L100-109 |
| fv_8F6L_AB | PDB 8F6L chains A(H)/B(L) | IGHV1-46-like + IGKV4-1-like (DIVMTQSPDSLAVSLGERATINC) | QVQLVQSGAEVKKPGASVKVSCKASgfslwtyaMNWVRQAPGQGLEWMGIisddgttYYATWAKGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCardagaggvqdyltlWGQGTLVTVSS | DIVMTQSPDSLAVSLGERATINCQASrniynfLAWYQQKPGQPPKLLIYsasTLASGVPDRFSGSGSGTDFTLTISSLQAEDVAVYYCqqgssnsnidnpFGGGTKVEIK | H1-22, 13-13, H36-49, 9-9, H59-95, 15-15, H111-121, /0, L1-23, 11-11, L35-48, 8-8, L57-88, 12-12, L101-110 |
| fv_9GWZ_HL | PDB 9GWZ chains HHH(H)/LLL(L) | IGHV4-59-like + IGKV4-1-like (DIVLTQSPDSLAVSLGERATINC) | EVQLQESGPGLVKPSETLSLTCTVSgfsltnyaVSWVRQPPGKGLEWLGVmwagggtNYNSVFKSRLTISKDNSKNQVSLKLSSVTAADTAVYYCarerpltgvmdyWGQGTLVTVSS | DIVLTQSPDSLAVSLGERATINCRASesvdysgnsfMHWFQQKPGQPPKLLIYrasNLESGIPDRFSGSGSRTDFTLTISSLQAEDVAVYYChqsnedpptFGGGTKVEIK | H1-22, 13-13, H36-49, 9-9, H59-95, 12-12, H108-118, /0, L1-23, 15-15, L39-52, 8-8, L61-92, 9-9, L102-111 |
| ref_fv_adalimumab | apheris/rfd4/Data/03-fv/framework_adalimumab.cif | IGHV3-9-like (VH3) + IGKV1-like -- same families as trastuzumab | EVQLVESGGGLVQPGRSLRLSCAASgftfddyaMHWVRQAPGKGLEWVSAitwnsghiDYADSVEGRFTISRDNAKNSLYLDMNSLRAEDTAVYYCakvsylstassldyWGQGTLVTVSS | DIQMTQSPSSLSASVGDRVTITCRASqgirnyLAWYQQKPGKAPKLLIYaasTLQSGVPSRFSGSGSGTDFTLTISSLQPEDVATYYCqrynrapytFGQGTKVEIK | H1-22, 13-13, H36-49, 10-10, H60-96, 14-14, H111-121, /0, L1-23, 11-11, L35-48, 8-8, L57-88, 9-9, L98-107 |
| ctrl_vhh_9MTE_A | PDB 9MTE chain A | VHH IGHV3-like, FR1 GGGLVQAGG | VQLVESGGGLVQAGGSLRLSCVASgipfstrtMAWYRQPPGNERDWVATirsgapVYADSVKGRFTVSRDNAKNTLYLQMNSLEPEDTATYYCwdvngdiWGQGTLVTV |  | H1-21, 13-13, H35-48, 8-8, H57-93, 7-7, H101-109 |
| ctrl_vhh_9MTE_B | PDB 9MTE chain B | VHH IGHV3-like, FR1 GGGVVQPGG | VQLVESGGGVVQPGGSLRLSCAASgsifgfndMAWYRQAPGKQRELVALisrvgvtSSADSVKGRFTISRVNAKDTVYLQMNSLRPEDTALYYCymdqrldgstlayWGQGTLVTVS |  | H1-21, 13-13, H35-48, 9-9, H58-94, 13-13, H108-117 |
| ref_vhh_caplacizumab | apheris/rfd4/Data/03-vhh/pmvk_caplacizumab.cif.zst | VHH IGHV3-like, FR2 WFRQAPGKGRELVA | EVQLVESGGGLVQPGGSLRLSCAASgrtfsynpMGWFRQAPGKGRELVAAisrtygnTYPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCaaagvraedgrvrtWGQGTQVTVSS |  | H1-22, 13-13, H36-49, 9-9, H59-95, 14-14, H110-120 |
| vhh_2P45_B | PDB 2P45 chain B | VHH IGHV3-like, FR2 WFRQ/EREGVA | QVQMVESGGGLVQAGGSLRLSCAASgyaytyiyMGWFRQAPGKEREGVAAmdsggggtLYADSVKGRMTISRDKGKNTVYLQMDSLKPEDTATYYCaaggyelrdrtygqWGQGTQVTVSS |  | H1-22, 13-13, H36-49, 10-10, H60-96, 14-14, H111-121 |
| vhh_2X1O_A | PDB 2X1O chain A | VHH IGHV3-like, FR2 WFRQ/EREFVG | VQLQESGGGLVQAGGSLRLSCAAAgrnlrmyrMGWFRQAPGKEREFVGTmvwssdtiYYADSVKGRFIISRDNAKNTVYLQMNSLKPEDTAVYYCaagagwagtmtdynyWGQGTQVTVS |  | H1-21, 13-13, H35-48, 10-10, H59-95, 15-15, H111-120 |
| vhh_4B41_A | PDB 4B41 chain A | VHH IGHV3-like, FR2 WYRQ/QRELVA | VQLVESGGGLVQPGGSLRLSCAASrsiisnnaMGWYRQAPGKQRELVARissggrtTYADSVKGRFTISRDNAKTTVYLQMNSLKPEDTAVYYCnaaslvrgpldhWGQGTQVTVSS |  | H1-21, 13-13, H35-48, 9-9, H58-94, 12-12, H107-117 |
| vhh_4B50_A | PDB 4B50 chain A | VHH IGHV3-like, human-like FR1 EVQLVESGGGLVQPGGSLRLSCAAS | EVQLVESGGGLVQPGGSLRLSCAASgsissvdvMSWYRQAPGKQRELVAFitdrgrtNYKVSVKGRFTISRDNSKNMVYLQMNSLKPEDTADYLCraesrtswsspspldvWGRGTQVTVSS |  | H1-22, 13-13, H36-49, 9-9, H59-95, 16-16, H112-122 |
| vhh_5IVN_A | PDB 5IVN chain A | VHH IGHV3-like, FR2 WFRQ/EREGVSC | QVQLVESGGGLVQPGGSLTLSCTASgftldhydIGWFRQAPGKEREGVSCinnsdddtYYADSVKGRFTIFMNNAKDTVYLQMNSLKPEDTAIYYCaeargckrgryeydfWGQGTQVTVSS |  | H1-22, 13-13, H36-49, 10-10, H60-96, 15-15, H112-122 |
| vhh_5J1S_C | PDB 5J1S chain C | VHH IGHV3-like, FR1 QVQLVETGGG | QVQLVETGGGLVQAGGSLRLSCAASgnifsfnvMGWYRQAPGKQRELVAAitsgdttTYADSVQGRFTISRDNAKNAVYLQMNSLTPEDTAVYFCnarrnpingpyyttayWGQGTQVTVSS |  | H1-22, 13-13, H36-49, 9-9, H59-95, 16-16, H112-122 |
| vhh_5M2W_B | PDB 5M2W chain B | VHH IGHV3-like, FR1 LVQAGGTLKLSC | QVQLVESGGGLVQAGGTLKLSCAASgsisgivvMAWYRQAPGKQRELVASitsggttNYADSVKGRFTISKDNAENTLYLRMNSLKPEDTAVYYCkaffrrdyvgydyWGQGTQVTVSS |  | H1-22, 13-13, H36-49, 9-9, H59-95, 13-13, H109-119 |
| vhh_5VLV_A | PDB 5VLV chain A | VHH IGHV3-like, FR1 LVQTGDSLRLSC | EVQLVESGGGLVQTGDSLRLSCAASgrtytpyaMAWFRQAPGKEREFVAGiggidgtaAYADSVRGRATISRDSAKKTVYLQMNSLKPEDTAVYSCatrasmqvltsprvypiWGRGTQVTVSS |  | H1-22, 13-13, H36-49, 10-10, H60-96, 17-17, H114-124 |
| vhh_6APO_A | PDB 6APO chain A | VHH IGHV3-like, atypical FR2 WYRQVPLGKGRELLS | KLQESGGGLVQAGESLRLSCAVPpevfdirtvAWYRQVPLGKGRELLSSitpwnktTYEDSVKDRFTISRDNAKYTVYLQMNDLKPEDTAVYYCaqgwgiasmryWGQGTQVTVSS |  | H1-20, 14-14, H35-48, 9-9, H58-94, 11-11, H106-116 |
| vhh_6HER_B | PDB 6HER chain B | VHH IGHV3-like, FR2 WFRQ/GREFVA | QVQLQESGGGLVQPGGSLRLSCAASgrtfssynMGWFRQAPGKGREFVASitssgdksDYTDSVKGRFTISRDNAKNTMYLQMNNLKPEDTATYYCarglgiyiirarggydhWGQGTQVTVSS |  | H1-22, 13-13, H36-49, 10-10, H60-96, 17-17, H114-124 |
| vhh_6JB9_A | PDB 6JB9 chain A | VHH IGHV3-like, FR1 DVQLVESGGGSVQAGG | DVQLVESGGGSVQAGGSLRLSCAASgstdsieyMTWFRQAPGKAREGVAAlythtgntYYTDSVKGRFTISQDKAKNMAYLRMDSVKSEDTAIYTCgatrkyvpvrfaldqssydyWGQGTQVTVSS |  | H1-22, 13-13, H36-49, 10-10, H60-96, 20-20, H117-127 |
| vhh_6XYF_A | PDB 6XYF chain A | VHH IGHV3-like, FR2 WFRQ/AREFAA | VQLQESGGGLVQAGGSLRLSCAASggtfksggMAWFRQAPGKAREFAAGiswsggstDYEDSVKGRFTISRDNAKNTMYLQMNSLKPEDTAVYYCaaarrfragvvtraddvdyWGKGTQVTVSS |  | H1-21, 13-13, H35-48, 10-10, H59-95, 19-19, H115-125 |
| vhh_7AQZ_C | PDB 7AQZ chain C | VHH IGHV3-like, FR2 WFRQAPEKEREFVA | QVQLQESGGGLVQAGGSLRLSCEASgltfsnyaMAWFRQAKEREFVAGiswtgsrtYYADSVRGRFTTSRDGHKNTVYLQMNDLKPEDTAVYLCaadllgsgkdgtsvyeyWGQGTQVTVSS |  | H1-22, 13-13, H36-47, 10-10, H58-94, 17-17, H112-122 |
| vhh_7R20_B | PDB 7R20 chain B | VHH IGHV3-like, FR2 WYRQ/QREVVA | EVQLQESGGGLVQAGGSLRLSCAASgsifsgnaMGWYRQAPGKQREVVAVisagnssNYVDSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCnvvkrgpqwgmeWGKGTLVTVSS |  | H1-22, 13-13, H36-49, 9-9, H59-95, 12-12, H108-118 |
| vhh_7UTG_A | PDB 7UTG chain A | VHH IGHV3-like, FR2 WYRQ/ERELVA | VQLVESGGGLVQAGGSLRLSCAASgsifvnnaMGWYRQAPGKERELVAAisasggstYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCaadqdgypyeyWGQGTQVTVSS |  | H1-21, 13-13, H35-48, 10-10, H59-95, 11-11, H107-117 |
| vhh_8CDS_A | PDB 8CDS chain A | VHH, most divergent FR1 GGGVVQTGDSLMLSCTVS | QVQLVESGGGVVQTGDSLMLSCTVSghtidnwaMGWFRQTPRKQREFVAAidkrgtgaIYGNSVKGRFTVSRDNAKNMVYLRMNSLKPEDTAVYFCavdqlnaglgdvsydydyWGQGTQVTVSS |  | H1-22, 13-13, H36-49, 10-10, H60-96, 18-18, H115-125 |
| vhh_8CDT_A | PDB 8CDT chain A | VHH IGHV3-like but CLASSIC VH3 FR2 WVRQAPGKGLEWVS | DVQLVESGGGLVQPGGSLRLSCVASgftfsrvgMGWVRQAPGKGLEWVSDinasggtGYADSVKGRFAISRDNAKNTLYLQMNRLKPEDTAVYYCakmmdtamieagiikpAGQGTQVTVS |  | H1-22, 13-13, H36-49, 9-9, H59-95, 16-16, H112-121 |
| vhh_8E2N_B | PDB 8E2N chain B | VHH IGHV3-like, divergent FR1 LVPAGG | EVQLVESGGGLVPAGGSLRLSCTTSerafrsnaMGWFRQAPGKEREFVAAvsvlswsgdsaVVADSVAGRFTIFRDNAKNTVYLQMNSLKPEDTAVYYCngasdigalqsgasswsWGHGTQVTVSS |  | H1-22, 13-13, H36-49, 13-13, H63-99, 17-17, H117-127 |
