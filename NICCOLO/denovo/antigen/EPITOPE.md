# HER2 domain IV — antigen-only design input and measured trastuzumab epitope

P0-worker (b), de novo antibody campaign vs HER2 at the trastuzumab epitope.
Everything below marked **[measured]** was computed here from the crystal coordinates.
Everything marked **[assumed]** is inherited from the brief / upstream docs and was *not*
independently verified. Nothing here is a literature value unless it carries a citation.

Scripts that produced this file (both re-runnable, both in this directory):

| File | What it does |
| --- | --- |
| `epitope_analysis.py` | contacts, per-residue burial, writes `epitope_contacts.csv` |
| `verify_antigen.py` | re-parses `her2_dIV_CD.cif` and asserts every property claimed below |

Run with `/home/hackathon/.venvs/abmodels/bin/python -I` (the only interpreter on this
box with numpy + biotite; the system `python3` has neither).

---

## 1. Deliverable: antigen-only mmCIF

**`/home/hackathon/hackathon-team-3/NICCOLO/denovo/antigen/her2_dIV_CD.cif`**

Built from `/home/hackathon/hackathon-team-3/apheris/rfd4/Data/12-her2/her2_fv_1n8z.cif`
(PDB 1N8Z, cropped) by **surgical text filtering**, not a library round-trip:

* `_atom_site` rows with `label_asym_id` ∉ {C, D} deleted → trastuzumab VL (A) and VH (B) gone.
* `_entity`, `_entity_poly`, `_entity_poly_seq` rows for entities 1 and 2 (the two
  trastuzumab chains) deleted, so the header no longer advertises an antibody. Entity ids
  **3 and 4 were deliberately not renumbered** — renumbering would have forced an edit to
  `label_entity_id` in every atom row for no gain.
* `_atom_site.id` renumbered 1..837 and `_atom_site.atom_array_index` 0..836 so the file is
  internally self-consistent.
* Everything else — `_entry`, `_chem_comp`, `_chem_comp_bond`, `_cell`, `_chem_comp_atom`,
  and the custom `ig_region` column — copied byte-for-byte.

### Verified inventory [measured] — `verify_antigen.py`, all assertions pass

| Property | Value |
| --- | --- |
| atom_site rows | 837 (= all 837 C/D rows in the source; none dropped, none added) |
| chains, `label_asym_id` | `['C', 'D']` — no A/B atoms survive |
| chains, `auth_asym_id` | `['C', 'D']` — checked separately, both columns clean |
| `label_entity_id` | `['3', '4']` only |
| chain C | **92 residues, 489–580, contiguous, no gaps** |
| chain D | **17 residues, 591–607, contiguous, no gaps** |
| `label_seq_id == auth_seq_id` | true for all 837 rows (brief's premise confirmed) |
| NaN / non-finite coordinates | **0** |
| coordinate bbox | min `[1.69, 80.18, 88.11]`, max `[68.95, 116.51, 125.75]` |
| coordinate strings vs source | **0 differ** — no silent reformatting of any number |
| header residue | no `DIQMTQSPSS` / `EVQLVESGGG` sequence, no `1 polymer` / `2 polymer`, no strand id A/B |
| independent re-parse | biotite `pdbx.get_structure` → 837 atoms, chains C/D, 109 residues, 0 NaN |

Trastuzumab is absent from both the coordinates and the metadata, so the model cannot
rebuild it from the input.

---

## 2. Measured epitope [measured]

Heavy atoms only, hydrogens excluded, 4.5 Å any-heavy-atom-to-any-heavy-atom cutoff.
mmCIF reader and the chunked distance kernel are adapted from
`/home/hackathon/wt/jahn/JAHN/pipeline/abpipe/structure.py` (`read_chains`,
`contact_residues`), copied into `epitope_analysis.py` — that file was not modified.

### 2a. HER2 residues within 4.5 Å of trastuzumab Fv (chains A+B) — **n = 19**

```
C/557 P   C/558 E   C/560 D   C/561 Q   C/569 K   C/570 D   C/571 P   C/572 P
C/573 F   C/579 P   D/591 I   D/593 K   D/596 D   D/598 E   D/600 A   D/602 Q
D/603 P   D/604 C   D/605 P
```

### 2b. HER2 residues within 4.5 Å of CDRH3 only (chain B, 99–108) — **n = 7**

```
C/570 D   C/571 P   C/572 P   C/573 F   C/579 P   D/591 I   D/593 K
```

CDRH3 sequence confirmed in the file as `W99 G100 G101 D102 G103 F104 Y105 A106 M107 D108`
= `WGGDGFYAMD`, matching trastuzumab. So CDRH3 reaches **7 of the 19** epitope residues;
the other 12 come from the remaining five CDRs.

### 2c. Agreement with the published 19-residue 557–605 epitope

**Full agreement on count and span.** My independent measurement gives exactly **19**
residues, lowest 557, highest 605 — identical to the published figure quoted in the brief
(Jahn et al., 4.5 Å). **No disagreements found.**

Caveat on the strength of that agreement: the brief supplied only the *count* and *span*
of the published set, not its residue list, so this is an agreement on two summary
statistics, not a residue-by-residue reconciliation. Given the same structure, the same
cutoff and the same heavy-atom convention, that is expected rather than surprising — the
useful conclusion is that no methodological discrepancy (hydrogen handling, altloc,
chain-split artefact) has crept in.

Two residues sit right on the boundary and would move under a slightly different cutoff:
`C/579` (`d_min` = **4.50 Å**, exactly at the cutoff — it is in only because of float
rounding) and `D/605` (4.15 Å, only 10.7 % buried). A 4.0 Å cutoff would drop both; a 5.0 Å
cutoff would add `C/575` (5.35 Å) and `D/601` (5.12 Å). Any downstream recall/precision
metric against a reference epitope should state its cutoff.

One trap in the reference code worth recording: `contact_residues()` returns bare residue
*numbers* and discards the chain id. On this split antigen that happens to be safe because
C (489–580) and D (591–607) do not overlap in numbering — but it would silently merge the
two chains on any antigen where they did. My rollup keys on `(chain, resnum)`.

### 2d. Buried surface area [measured — real SASA, not a proxy]

Shrake–Rupley SASA via `biotite.structure.sasa`, `vdw_radii="Single"`, `point_number=500`.
Per-residue burial = SASA(HER2 alone, chains C+D) − SASA(HER2 inside the C+D+A+B complex).
**This is a measured ΔSASA, not a geometric proxy** — a SASA implementation was available,
so no proxy was needed.

Total HER2-side interface area buried by the Fv: **747 Å²**.

| res | nFv | nH3 | d_min Fv | d_min H3 | ΔSASA Å² | ΔSASA sc | free SASA | % buried | closest atom |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| C/558 E | 22 | 0 | 3.13 | 7.63 | **90.1** | 89.6 | 101.6 | **88.6** | OE1 |
| C/571 P | 25 | 5 | 2.89 | 3.32 | **76.4** | 47.3 | 82.0 | **93.2** | O |
| C/573 F | 27 | 16 | 3.42 | 3.44 | **74.4** | 74.4 | 85.3 | **87.3** | CZ |
| C/572 P | 32 | 8 | 3.11 | 4.00 | **70.3** | 60.4 | 78.1 | **90.1** | N |
| D/603 P | 20 | 0 | 3.56 | 6.65 | **68.5** | 61.8 | 84.3 | **81.2** | O |
| C/560 D | 15 | 0 | 2.60 | 6.02 | 45.1 | 36.5 | 67.4 | 66.8 | OD2 |
| C/557 P | 12 | 0 | 3.51 | 12.96 | 40.6 | 38.4 | 68.8 | 59.0 | CB |
| D/602 Q | 11 | 0 | 3.05 | 7.88 | 40.0 | 40.0 | 61.8 | 64.8 | OE1 |
| C/569 K | 9 | 0 | 2.96 | 7.28 | 37.7 | 37.7 | 87.9 | 42.9 | NZ |
| D/591 I | 4 | 4 | 3.74 | 3.74 | 34.3 | 34.3 | 165.4 | *20.7* | CG2 |
| D/598 E | 3 | 0 | 3.52 | 15.00 | 33.9 | 33.9 | 170.9 | *19.8* | OE1 |
| D/593 K | 7 | 7 | 2.42 | 2.42 | 27.2 | 27.2 | 47.0 | 57.9 | NZ |
| C/561 Q | 9 | 0 | 3.71 | 10.58 | 26.2 | 26.2 | 47.0 | 55.6 | NE2 |
| C/570 D | 12 | 12 | 3.50 | 3.50 | 21.6 | 18.6 | 27.9 | 77.3 | O |
| D/600 A | 3 | 0 | 3.82 | 10.26 | 13.3 | 13.3 | 19.8 | 67.1 | CB |
| D/605 P | 2 | 0 | 4.15 | 12.14 | 13.3 | 3.6 | 124.2 | *10.7* | CA |
| C/579 P | 1 | 1 | 4.50 | 4.50 | 8.7 | 8.7 | 43.0 | *20.2* | CG |
| D/596 D | 4 | 0 | 2.99 | 12.78 | 8.6 | 8.6 | 53.9 | *15.9* | OD2 |
| D/604 C | 3 | 0 | 3.81 | 9.04 | 6.9 | 0.0 | 56.7 | *12.1* | O |

Full table incl. the three sub-threshold residues (C/532, C/575, D/601, which bury area
but make no 4.5 Å contact) is in **`epitope_contacts.csv`**.

**The contact-count trap, confirmed on this structure.** Ranking by `nFv` and ranking by
ΔSASA give different answers. `C/572` has the most contacts (32) but is 4th by ΔSASA.
`C/558` buries the most area of any residue in the epitope (90 Å², 89 % of its free
surface) yet has *zero* CDRH3 contacts — contact-count ranking against CDRH3 would have
discarded the single most buried residue at the interface. Conversely `D/591` and `D/598`
look substantial by raw ΔSASA (34 Å² each) but are only ~20 % buried: they are large,
mostly-solvent-exposed side chains that happen to graze the antibody, i.e. periphery, not
core. The `% buried` column is what separates the two.

### 2e. Sub-patch decomposition by CDR [measured]

| Antibody region | HER2 residues contacted (contact count) |
| --- | --- |
| VL, chain A (all) | C/560:6 C/569:9 C/571:20 C/572:24 D/596:4 D/598:3 D/600:3 D/602:11 D/603:20 D/604:3 D/605:2 |
| VH, chain B (all) | C/557:12 C/558:22 C/560:9 C/561:9 C/570:12 C/571:5 C/572:8 C/573:27 C/579:1 D/591:4 D/593:7 |
| CDRH3 B99–108 | C/570:12 C/571:5 C/572:8 C/573:16 C/579:1 D/591:4 D/593:7 |
| CDRH2 B50–58 | C/557:12 C/558:19 C/560:5 C/573:2 |
| CDRH1 B31–35 | C/558:1 C/573:9 |
| CDRL3 A89–97 | C/560:6 C/569:9 C/571:18 C/572:24 |
| CDRL1 A24–34 | C/571:2 D/596:4 D/598:3 D/600:3 D/602:11 |

The epitope splits cleanly into **two sub-patches 14.2 Å apart (residue-centroid to
residue-centroid)**:

* **VH-side patch**, the chain-C loops 557–561 + 569–573. Buries **491 Å²**, two thirds of
  the interface. Hit by CDRH1, CDRH2 and CDRH3.
* **VL-side juxtamembrane patch**, chain D 596–605. Buries **185 Å²**. Hit by CDRL1 and
  CDRL2/3 only — **no CDRH3 contact at all**.
* `C/569`/`C/571`/`C/572` are the shared hinge, touched by both CDRL3 and CDRH3.

Full epitope diameter: 22.7 Å (C/558 to D/603 centroids).

### 2f. The unresolved 581–590 loop is a design hazard [measured]

`C/580` CA → `D/591` CA is only **6.7 Å**, yet **10 residues (581–590) must bridge that
gap**. Those residues are real protein that is simply not in the crystal, and they have to
bulge out into the volume beside the interface. Distance from each candidate patch centroid
to the midpoint of the chain break:

| patch | distance to gap midpoint |
| --- | --- |
| Set B core (VH-side) | **16.8 Å** — safest |
| Set C patch (VL-side) | 13.2 Å |
| Set A, Apheris's 4 atoms | **9.0 Å** — closest to the hole |

A de novo binder steered into the gap region will look clash-free against this CIF and
clash against the real receptor. Set B is the most defensible on this criterion; designs
driven by Set A should be re-checked against a model that rebuilds 581–590 before anything
is ordered.

---

## 3. Three hotspot hypotheses for `C_HOT`

All atoms below were checked to **exist in `her2_dIV_CD.cif`** — see §4. Selector syntax
follows `apheris/rfd4/foundry.md` and `02-miniprotein.md`: `sel('CHAIN/*/RESNUM/ATOM')`,
and selections use the original mmCIF `label_seq_id`, which here equals `auth_seq_id`.

### Set A — control: reproduce Apheris's four atoms verbatim

```yaml
    C_HOT:
      true:
        - select: "sel('D/*/593/NZ')"
        - select: "sel('C/*/573/CD1')"
        - select: "sel('D/*/591/CG2')"
        - select: "sel('C/*/572/CD')"
```
```json
"C_HOT": {"true": [{"select": "sel('D/*/593/NZ')"}, {"select": "sel('C/*/573/CD1')"},
                   {"select": "sel('D/*/591/CG2')"}, {"select": "sel('C/*/572/CD')"}]}
```

**Why:** the control. Identical to `Data/12-her2/her2_cdrh3.yaml`, so any difference in
outcome between A and B/C is attributable to the hotspot choice and not to the input file.
Measured check: all four atoms are genuinely CDRH3-proximal (`d_min` to CDRH3 = 2.42, 3.44,
3.74, 4.00 Å), so the upstream description of them is accurate.

**Measured weaknesses — the reason A is only a control.** (i) It is a *CDRH3*-centric set,
tuned for redesigning one loop on a fixed framework, not for a de novo binder that has to
build a whole interface. (ii) It omits `C/558`, the most buried residue in the epitope
(90 Å², 89 %), and `C/571` (93 % buried). (iii) `D/591/CG2` is a weak hotspot by every
measure taken here: 20.7 % buried, 4 contacts, and it is the N-terminal residue of chain D
sitting right at the unresolved-loop break, so its large apparent free SASA (165 Å²) is
partly a cropping artefact. (iv) Patch span 15.7 Å / mean 10.9 Å — the loosest of the
three, and the closest to the 581–590 hole (9.0 Å).

### Set B — burial-ranked core of the VH-side patch (**recommended primary**)

```yaml
    C_HOT:
      true:
        - select: "sel('C/*/558/OE1')"
        - select: "sel('C/*/560/OD2')"
        - select: "sel('C/*/571/CB')"
        - select: "sel('C/*/572/CB')"
        - select: "sel('C/*/573/CZ')"
```
```json
"C_HOT": {"true": [{"select": "sel('C/*/558/OE1')"}, {"select": "sel('C/*/560/OD2')"},
                   {"select": "sel('C/*/571/CB')"}, {"select": "sel('C/*/572/CB')"},
                   {"select": "sel('C/*/573/CZ')"}]}
```

**Why:** these are the five residues trastuzumab buries hardest, and the atom chosen in each
is that residue's most-buried atom. Per-residue: `C/558` 90.1 Å² / 88.6 %, `C/571`
76.4 Å² / 93.2 %, `C/573` 74.4 Å² / 87.3 %, `C/572` 70.3 Å² / 90.1 %, `C/560` 45.1 Å² /
66.8 %. **Every member is ≥ 67 % buried**, which is the property that distinguishes an
interface core from interface periphery — and the one the ΔSASA-vs-contact-count comparison
in §2d shows contact counting misses. Together they account for **356 Å², 48 % of the
whole 747 Å² interface**, from five atoms. Patch span 13.0 Å, mean 6.8 Å — the most compact
of the three, so a single binder loop or helix face can reach all five. It is also the
furthest from the unresolved 581–590 loop (16.8 Å).

Atom-level ΔSASA for the picks: `C/558/OE1` 43.5 Å² (the most buried single atom in the
epitope, 3.13 Å contact), `C/560/OD2` 32.8 Å² (**2.60 Å** — the shortest polar contact on
the chain-C side, almost certainly a salt bridge or strong H-bond), `C/573/CZ` 27.3 Å²
(ring tip; note Apheris picks `CD1` at only 13.8 Å²), `C/572/CB` 26.1 Å², `C/571/CB`
25.1 Å². Chemically this set is well mixed — two carboxylates to pair with binder
Arg/Lys/Ser, and three apolar carbons (two proline rings plus a Phe edge) for a hydrophobic
anchor — which is what the campaign wants if the designed binder is to carry its own
affinity rather than copy CDRH3.

Backbone alternative if a hydrogen-bond-driven design is wanted instead: `C/571/O` is the
single most buried backbone atom (27.6 Å², 2.89 Å contact) and can be swapped in for
`C/571/CB`.

### Set C — VL-side juxtamembrane sub-patch (orthogonal competitor)

```yaml
    C_HOT:
      true:
        - select: "sel('D/*/596/OD2')"
        - select: "sel('D/*/598/OE1')"
        - select: "sel('D/*/600/CB')"
        - select: "sel('D/*/602/OE1')"
        - select: "sel('D/*/603/CG')"
```
```json
"C_HOT": {"true": [{"select": "sel('D/*/596/OD2')"}, {"select": "sel('D/*/598/OE1')"},
                   {"select": "sel('D/*/600/CB')"}, {"select": "sel('D/*/602/OE1')"},
                   {"select": "sel('D/*/603/CG')"}]}
```

**Why this is a genuinely different hypothesis, not a reshuffle of B.** Every residue here
is contacted by trastuzumab's **light** chain and by **zero** CDRH3 atoms (`d_min` to CDRH3
= 12.78, 15.00, 10.26, 7.88, 6.65 Å). Its centroid is **14.2 Å** from Set B's, and it is
all chain D where B is all chain C. So a binder driven to this patch would **sterically
compete with trastuzumab** — it occupies the VL footprint, within the same 22.7 Å epitope,
so the two cannot bind at once — while making an essentially **non-overlapping set of
residue contacts**. That is the useful experiment: it tests whether the trastuzumab epitope
can be blocked through a different chemistry, and it gives the campaign a candidate whose
escape-mutation profile should differ from trastuzumab's.

Burial support: `D/603` is the 5th most buried residue in the entire epitope (68.5 Å²,
81.2 %) and `D/602` is 40.0 Å² / 64.8 %; `D/600/CB` is 67.1 % buried. The patch buries
185 Å² in total — a quarter of the interface, enough footprint for a miniprotein or nanobody
paratope. Span 14.0 Å, mean 7.3 Å.

**Honest weaknesses.** Two members, `D/596` (15.9 %) and `D/598` (19.8 %), are
peripheral by fractional burial — they are in because they sit in this patch and make close
polar contacts (2.99 and 3.52 Å), not because they are buried. Set C is therefore a
*weaker-anchored* hypothesis than B by design, and 13.2 Å from the 581–590 hole. Expect
lower success rate; run it for differentiation, not as the primary.

### Recommended order

**B first** (most buried, most compact, furthest from the missing loop), **A as the
control**, **C for an orthogonal competitor**. Keep hotspot resampling off when a specific
contact is essential (`apheris/rfd4/foundry.md`, antibody notes).

---

## 4. Atom-existence verification [measured]

A `C_HOT` selector naming a missing atom is a silent failure, so every atom in §3 was
looked up by `(label_asym_id, label_seq_id, label_atom_id)` in the **output** CIF, not the
input. `verify_antigen.py` re-runs this check.

| Set | Selector atom | Exists in `her2_dIV_CD.cif` |
| --- | --- | --- |
| A | `C/572/CD` | yes |
| A | `C/573/CD1` | yes |
| A | `D/591/CG2` | yes |
| A | `D/593/NZ` | yes |
| B | `C/558/OE1` | yes |
| B | `C/560/OD2` | yes |
| B | `C/571/CB` | yes |
| B | `C/572/CB` | yes |
| B | `C/573/CZ` | yes |
| C | `D/596/OD2` | yes |
| C | `D/598/OE1` | yes |
| C | `D/600/CB` | yes |
| C | `D/602/OE1` | yes |
| C | `D/603/CG` | yes |

Also verified as existing (named above as alternatives): `C/571/O`.
**14 / 14 present. No silent-failure selectors in any of the three sets.**

---

## 5. Literature check — AMASS cache

Searched the read-only cache at
`/home/hackathon/hackathon-team-3/LEO/tasks/experiment-ledger/cache/` (32 cached responses).
No new network calls were made.

**No residue-level HER2 evidence is available there.** Every cached response is a
BiomedCore/PatentCore record carrying only `title`, `abstract`, `pmid`, `journal`,
`meshTerms` etc. — no full text. A regex sweep of all 32 files for HER2 domain-IV residue
numbers (550–607) in any title or abstract returned **four hits, all false positives** (a
567-sequence Tm dataset, a 600-fold rate-constant ratio). So the cache **neither supports
nor contradicts** the specific atom picks in §3; the hotspot sets rest on the structure
measured here, not on literature.

What the cache *does* contain that is relevant, with provenance:

* `biomedcore__trastuzumab-fab-her2-crystal-structure-contact-residues-1n8z.json`, record 2
  — **Cho et al., "Structure of the extracellular region of HER2 alone and in complex with
  the Herceptin Fab"** (Nature; the 1N8Z source paper). Abstract: *"Herceptin binds to the
  juxtamembrane region of HER2, identifying this site as a target for anticancer
  therapies."* This **supports Set C**: the juxtamembrane segment (chain D here, HER2
  591–607, immediately N-terminal to the transmembrane helix) is called out by the authors
  of the structure as the therapeutically relevant site, and Set C is the sub-patch inside
  it. It is the only citation in the cache that speaks to *where* on HER2 to aim.
* `biomedcore__trastuzumab-alanine-scanning-mutagenesis-her2-binding-hot-sp.json` — contains
  **no trastuzumab/HER2 alanine scan**. The query returned generic hot-spot-prediction
  methodology papers. Two are still useful as *criteria*:
  * record 1, *"Relationship between hot spot residues and ligand binding hot spots in
    protein-protein interfaces"*: *"a residue that is found by alanine scanning to
    contribute little to binding rarely interacts with hot spot regions"*, and hot spots are
    *"a residue or cluster of residues that makes a major contribution to the binding free
    energy."* The **cluster** framing is why Sets B and C are compact patches (span 13.0 and
    14.0 Å) rather than the top-5 ΔSASA residues globally, which would have straddled two
    patches 14 Å apart.
  * record 3 (PMID in file), *"Hot spot occlusion from bulk water…"*: *"a hypothesis … that
    considers that inaccessibility to the solvent must be a necessary condition to define a
    residue as a binding hot spot"* (the **O-ring hypothesis**). This is the published
    justification for ranking on **fractional** burial rather than contact count or raw
    ΔSASA — i.e. for Set B's ≥ 67 %-buried rule, and for flagging `D/591` (20.7 %),
    `D/598` (19.8 %), `D/596` (15.9 %) and `D/605` (10.7 %) as peripheral.
* `biomedcore__her2-antibody-interface-buried-surface-area-shape-complement.json` — generic
  interface-metric papers, no HER2 numbers. Nothing here to compare the measured 747 Å²
  against.

**Standing caveat:** no experimental ΔΔG for any HER2 residue was obtained. The three sets
are structure-derived hypotheses. "Hot spot" in this document means "buried and closely
contacted", which is a structural proxy for the energetic definition, and the two are
correlated but not identical.

---

## 6. Premises in the brief — checked

| Premise | Status |
| --- | --- |
| Source = 1N8Z cropped; A = VL 1–107, B = VH 1–120, C = HER2 489–580, D = HER2 591–607 | **[measured] correct** — C has 92 residues 489–580, D has 17 residues 591–607, A/B present in the source with 1747 heavy atoms between them |
| 581–590 unresolved | **[measured] correct** — and the C/580→D/591 CA distance is only 6.7 Å for 10 missing residues, so the loop bulges into the interface neighbourhood (§2f) |
| `label_seq_id == auth_seq_id` | **[measured] correct** — true for all 837 antigen rows |
| Apheris's CDRH3 hotspots are `D/593/NZ`, `C/573/CD1`, `D/591/CG2`, `C/572/CD` | **[measured] correct** — read back from `her2_cdrh3.yaml`; all four exist and all four are within 4.5 Å of CDRH3 |
| Published epitope = 19 residues, 557–605 | **[measured] reproduced exactly** — 19 residues, span 557–605 |
| Contact count ≠ hot spot; ranking by contacts is a trap | **[measured] confirmed** — §2d: most-contacts (`C/572`, 32) ≠ most-buried (`C/558`, 90 Å²), and `C/558` has zero CDRH3 contacts |
| Input validation rejects a numbering gap inside one chain | **[assumed] not verified** — this is a comment in `her2_cdrh3.yaml`; no validator was run. The C/D split is preserved regardless, so nothing depends on it being true |
| RFD4-Proteina rejects PDB, mmCIF only | **[assumed]** from `~/.claude/skills/rfd4-proteina/SKILL.md` ("The input is one annotated CIF, not a plain sequence, PDB, server path") — consistent, but no upload was attempted |
| SASA tool availability | **[measured] resolved in our favour** — the system `python3` has no numpy/biotite/Bio/freesasa/mdtraj/gemmi at all, but `/home/hackathon/.venvs/abmodels/bin/python` has numpy 2.5.4, biotite 1.8.0, biopython 1.88 and scipy 1.18.1. Real Shrake–Rupley SASA was used; **no geometric proxy was needed anywhere in this document** |

Nothing in the brief was found to be wrong. The one correction to upstream material is a
judgement, not a fact: Apheris's four atoms are a correct *CDRH3* hotspot set but a poor
*de novo binder* hotspot set, for the reasons in §3 Set A.

### Not verified / out of scope here

* No RFD4 run was submitted, so the CIF is verified as *well-formed and faithful*, not as
  *accepted by the live validator*. `Data/03-fv/target_F.cif` — an existing accepted
  target-only RFD4 input — carries only `data_`, the `_atom_site` loop and `_cell`, i.e. a
  strict subset of what `her2_dIV_CD.cif` has, which is the reason for confidence.
* No `contig` string is proposed; that belongs with whoever sets binder length. Note chain
  C is 92 residues and chain D is 17, and a two-chain target needs both started in the
  contig (cf. `/0` in the miniprotein spec).
* `C_CRD` / `C_SEQ` / `C_CRDGRP` conditions for the antigen are not written here; the
  analogous miniprotein pattern is `C_CRD: {true: [{select: "sel('C/*/*/[N,CA,C,O]')"}]}`
  plus the same for D, and `C_SEQ: {true: [{select: "sel('C')"}, {select: "sel('D')"}]}`.
