# Refolding-oracle calibration against HER2 domain IV

**Status:** rf3 (RoseTTAFold3) MEASURED, calibrated, and **usable** — clean separation,
positive 0.585 vs best negative 0.268 on pair ipTM. OpenFold3 (second, independent family)
MEASURED and **FAILS validation** on this panel: it ranks a known non-domain-IV VHH 4.4x
above trastuzumab, so it gets no vote in ranking (§6). We are therefore still effectively
single-sourced on rf3. All numbers in this file are **measured**, from real jobs on
the hosted JapanFold API, unless a cell says `pending` or `n/a`. Nothing here is predicted
or carried over from another experiment.

Harness: `refold.py` · panel: `controls.csv` (built by `make_controls.py`) ·
durable per-job/per-design ledger: `runs.jsonl` · structures: `structures/`.

---

## 1. The antigen construct (read this before comparing any number to any other)

**HER2 domain IV, residues 489-607, submitted as ONE contiguous chain of 119 aa.**

1N8Z is deposited with residues 581-590 unresolved, and the rfd4 input file
(`apheris/rfd4/Data/12-her2/her2_fv_1n8z.cif`) therefore splits the antigen into chain C
(489-580) and chain D (591-607). The brief framed this as a binary choice — one chain with
the gap "closed", or two chains. **It is a false dilemma.** The true 581-590 sequence is
recoverable from PDB **1S78** entity 1 (the full HER2 ectodomain), which reads
`...FCVARCPS|GVKPDLSYMP|IWKFPDEEG...`. So we fill the gap with the real ten residues
`GVKPDLSYMP` and submit one biologically correct, sequence-contiguous chain. No artificial
junction and no artificial chain break.

This is not cosmetic. We ran the identical trastuzumab Fv against both representations:

| antigen representation | trastuzumab WT pair ipTM | verdict |
|---|---|---|
| **one chain, 489-607, gap filled from 1S78** | **0.585** | clean separation from decoys |
| two chains, 489-580 + 591-607 (the rfd4 split) | 0.247 | **indistinguishable from the decoy floor (0.231)** |

**The split representation destroys the oracle.** Splitting the antigen at the loop costs
trastuzumab 0.34 pair ipTM and drops the true binder onto the non-binder pile. Always use
the single-chain construct. (`HER2_DIV_SINGLE` in `refold.py`.)

Chain ids used: `H` = heavy/VHH, `L` = light, `A` = antigen. Input order is `H, L, A`,
deliberately **different** from sorted order `A, H, L` — that difference is what makes the
indexing check in §3 meaningful.

---

## 2. Control panel — MEASURED rf3 (RoseTTAFold3) scores

Job `9a889c9cca62558160402d2806e54ce1`, model `rf3`, seed 1, 1 sample, 10 recycles, MSA on
(depth 1024), 7 targets in one job, **$0.020954 total**.

`pair_max` = max antibody-chain-to-antigen pair ipTM, correctly indexed (§3). This is the
headline metric. `ab_plddt` = mean pLDDT over the antibody chains only, from the returned
mmCIF.

| sequence_id | role | epitope on HER2 (measured from the PDB) | overall ipTM | pTM | **pair ipTM (max)** | pair ipTM (mean) | ab pLDDT | complex pLDDT | clash |
|---|---|---|---|---|---|---|---|---|---|
| `tras_wt_fv` | **POSITIVE** | domain IV 557-605 — **present** in construct | 0.7630 | 0.8274 | **0.5850** | 0.5815 | 79.57 | 0.8250 | no |
| `tras_scram_cdr` | **DECOY FLOOR** | n/a (CDRs shuffled) | 0.5623 | 0.6846 | **0.2312** | 0.2250 | 76.78 | 0.8019 | no |
| `pert_wt_fv` | negative, epitope absent | domain II 245-315 — **absent** from construct | 0.5975 | 0.7122 | **0.2482** | 0.2460 | 79.02 | 0.8218 | no |
| `vhh_27A05` | negative, epitope mismatch | domain I + III (35-500) | 0.1979 | 0.5607 | **0.1979** | 0.1979 | 78.75 | 0.8059 | no |
| `vhh_2Rs15d` | negative, epitope mismatch | domain I 109-217 | 0.2681 | 0.5842 | **0.2681** | 0.2681 | 76.69 | 0.8043 | no |
| `vhh_cablys3` | negative, irrelevant antigen | hen lysozyme (not HER2) | 0.2185 | 0.5916 | **0.2185** | 0.2185 | 76.05 | 0.7984 | no |
| `tras_wt_fv_2chainAg` | representation check | domain IV, antigen split in two | 0.6135 | 0.7299 | **0.2472** | 0.2267 | 79.63 | 0.8215 | no |

Per-chain detail for the positive: `H-A` 0.5779, `L-A` 0.5850, and the internal `H-L`
obligate VH/VL pair 0.9163.

### What the panel establishes

* **Separation is clean on pair ipTM.** Positive 0.585; the single best negative is 0.268
  (`vhh_2Rs15d`). Gap = **0.317**, i.e. the positive is 2.2x the decoy floor. Every
  negative, including three different *kinds* of negative, lands in a tight band
  0.198-0.268.
* **Overall ipTM must NOT be used to rank designs.** It is inflated by intra-antibody
  pairs. For an Fv it folds in the VH/VL pair (~0.88-0.92 for anything with a real
  framework), so the decoy scores 0.562 and the epitope-absent pertuzumab 0.598 against the
  true binder's 0.763 — a 0.17 gap where pair ipTM gives 0.34. Worse, **overall ipTM is not
  comparable across formats at all**: a VHH is a single chain, so its overall ipTM *is* its
  interface ipTM (0.198-0.268), while any Fv's overall ipTM is dominated by VH/VL
  (0.56-0.76). Ranking a mixed VHH + Fv pool by overall ipTM would reject every VHH on
  format alone. Pair ipTM is format-independent and is the only metric to rank on.
* **pLDDT does not discriminate** and must not be used as an affinity proxy. All seven
  targets sit at ab pLDDT 76-80 and complex pLDDT 0.80-0.83, binder and decoy alike. It is
  a *fold-sanity* filter only.
* **pTM tracks format, not binding** (0.56-0.59 for VHH complexes, 0.68-0.83 for Fv
  complexes). Do not threshold on it across formats either.

---

## 3. How the `pair_chains_iptm` indexing was verified

`pair_chains_iptm` comes back keyed by chain **position** as a string (`{"0":{"1":...}}`),
not by chain id, so the position-to-chain mapping has to be established, not assumed. The
brief's premise — **indexed by SORTED chain id** — is **CONFIRMED**, and here is the proof
rather than an assertion.

The `tras_wt_fv_2chainAg` target was built as the probe. Chains submitted in **input order
`H`(120 aa), `L`(107 aa), `A`(92 aa), `B`(17 aa)**; sorted order is `A, B, H, L`. The two
orders disagree on every position, so one job discriminates.

**Proof 1 — the 17-mer's per-chain pTM.** Returned `chains_ptm`:

```
{"0": 0.795626, "1": 0.111576, "2": 0.923673, "3": 0.927762}
```

Position 1 has pTM **0.1116**. A 17-residue fragment is the only chain in that complex that
can possibly have a pTM of 0.11; the other three are well-folded domains of 92-120 aa.
Chain `B` is the 17-mer.
* Sorted order `A,B,H,L` puts `B` at position **1**. ✅ matches.
* Input order `H,L,A,B` puts `B` at position 3, which scored 0.9278 — impossible for a
  17-mer, and it puts the 107-aa light chain at the 0.11 slot. ❌ excluded.

**Proof 2 — the two halves of one domain.** Cell `[0][1]` = **0.731**, far and away the
largest off-diagonal entry in that target (all others 0.205-0.247). Chains `A` (489-580)
and `B` (591-607) are the two covalently contiguous halves of a single HER2 domain, so they
must pack tightly. Sorted order puts them at positions 0 and 1. ✅ Under input order,
positions 0 and 1 would be `H` and `L`.

**Proof 3 — the obligate VH/VL pair.** Cell `[2][3]` = **0.917**, consistent with an
obligate VH/VL pairing. Sorted order puts `H` at 2 and `L` at 3. ✅ For the 3-chain Fv
targets the same holds: `[1][2]` = 0.916 (`H-L`) while position 0 carries the antigen's
lower pTM (0.74-0.79 vs 0.88-0.93 for the Ig domains).

**Proof 4 — chain ids round-trip through the structure.** `verify_chain_indexing()`
re-reads every returned mmCIF and checks each chain id's real residue count against what
was submitted. All 7/7 targets report `chain_id_verified = True`; every chain in the panel
has a distinct length, so this pins chain identity independently of any ordering
convention. Per-chain mean pLDDT corroborates it (chain `B`, the 17-mer, is the
lowest-confidence chain at 68.64).

`refold.py` does not hardcode positions: `chain_order()` returns `sorted(chain_spec)` and
`pair_iptm_lookup()` resolves chain ids through it, so a design with different chain ids
still indexes correctly.

---

## 4. Decision rule for scoring real designs — RELATIVE, never absolute

Absolute ipTM thresholds are not transferable (the published 0.5 "pass line" would accept
the scrambled decoy's own overall ipTM of 0.562). **Re-run the controls in the same batch as
the designs and score relative to them.** Define, per batch:

* `P` = `pair_iptm_ab_ag_max` of `tras_wt_fv` (the positive) — measured today: **0.585**
* `N` = max `pair_iptm_ab_ag_max` over all negative/decoy controls — today: **0.268**
  (`vhh_2Rs15d`)
* `D` = `pair_iptm_ab_ag_max` of the design's **own** scrambled-CDR decoy

**ACCEPT a design when all four hold:**

| # | criterion | today's value |
|---|---|---|
| 1 | `pair_iptm_ab_ag_max >= 0.85 * P` | >= 0.497 |
| 2 | `pair_iptm_ab_ag_max >= N + 0.15` | >= 0.418 |
| 3 | `pair_iptm_ab_ag_max >= D + 0.20` (own scrambled-CDR decoy) | >= decoy + 0.20 |
| 4 | `ab_mean_plddt >= 75` **and** `has_clash == false` | fold sanity only |

A convenient single number for triage is the midpoint `(P + N) / 2` = **0.43**; criteria 1
and 2 are deliberately stricter than that and bracket it. Criterion 3 is **per design and
mandatory** — Smorodina et al. (`hackathon-data/zenodo-18390239`) report ipTM often failing
to separate real from shuffled pairs, and the panel bears that out for *overall* ipTM
(decoy 0.562), so every design ships with its own shuffled-CDR twin. The 0.15/0.20 margins
are roughly half the 0.317 positive-to-negative gap measured here; tighten them if a
later, larger panel shows a narrower gap.

**Report `pair_iptm_ab_ag_max` for ranking. Never rank on overall ipTM, pTM or pLDDT**
(§2). If `P` in a later batch comes back below ~0.5, treat the whole batch as suspect
(wrong antigen representation, MSA failure) rather than re-baselining onto it.

---

## 5. Second oracle (independent family)

**OpenFold3** (`openfold3`), the OpenFold Consortium's AlphaFold3 reproduction — a
different lineage from RoseTTAFold3, and a genuine second family rather than a second seed.
It takes the identical chain-list YAML, so it is a drop-in: `refold.py run --model
openfold3`. It is set up, measured and reproducible — but it **failed validation on the
control panel and must not be used for ranking**; see §6 for the numbers and the rank
inversion. The ledger holds both models side by side and `refold.py table` prints them
together.

Chosen over the alternatives for this task:

* `af2ig` (AF2 initial-guess) is the best-suited scorer for an antibody-antigen pair, but it
  takes a **different input** — a *structure* plus the binder sequence, not a chain list —
  so it cannot score a sequence-only design until that design has a backbone. It is the
  right third oracle once designs come with coordinates; it needs its own target builder.
* `protenix` — viable hosted third family, not run here.
* `opendde-abag` — the antibody-antigen checkpoint, attractive for this exact problem.
  **For validation it must be run UNGUIDED:** epitope guidance pulls the antibody onto the
  epitope, so a guided run cannot confirm that the antibody binds there. Record this on any
  future OpenDDE result.
* **Boltz-2 ipTM is excluded from ranking.** Jahn's calibration on 273 SPR-labelled
  trastuzumab CDRH3 variants gave AUC 0.53 — no signal. Our panel independently shows *why*
  an overall-ipTM metric misleads (§2). Not retested here.

### Local docker path — verified invocation, deliberately NOT run

The single H100 is **shared and busy** (`nvidia-smi`: 78-99% utilisation, a teammate's
`quay.io/apheris/oslo-hackathon:apheris-protenix-v1.3.0-v2` container up), so no local job
was started, per the brief. Images confirmed present locally:

```
quay.io/apheris/oslo-hackathon:apheris-openfold3-v1.3.0-of3ob   11.1GB
quay.io/apheris/oslo-hackathon:apheris-protenix-v1.3.0-v2       11.2GB
quay.io/apheris/oslo-hackathon:apheris-openfold3-v1.3.0         8.8GB
quay.io/apheris/oslo-hackathon:apheris-mol1-v1.3.0              8.56GB
quay.io/apheris/oslo-hackathon:apheris-esmfold2-v1.4.0          9GB
```

Shape of the invocation (run only when the GPU is free, and check `docker ps` first):

```bash
nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv   # must be idle
docker ps                                                          # must be empty
docker run --rm --gpus all \
  -v "$PWD/in:/in:ro" -v "$PWD/out:/out" \
  quay.io/apheris/oslo-hackathon:apheris-openfold3-v1.3.0-of3ob \
  --input /in/target.yaml --output /out
```

Confirm the entrypoint's real flags with
`docker run --rm <image> --help` before trusting that command line.

---

## 6. Control panel — MEASURED OpenFold3 scores

Job `72e7cbe1a6110279ba40be1d62613d24`, model `openfold3`, 7 targets in one job, **$0.013992**
(193.7 billable s). Same panel, same single-chain antigen, same harness.

| sequence_id | role | overall ipTM | pTM | **pair ipTM (max)** | ab pLDDT | complex pLDDT | rf3 pair ipTM (for comparison) |
|---|---|---|---|---|---|---|---|
| `tras_wt_fv` | **POSITIVE** | 0.5519 | 0.6809 | **0.1413** | 84.53 | 0.8530 | 0.5850 |
| `tras_scram_cdr` | **DECOY FLOOR** | 0.4680 | 0.6177 | **0.0919** | 77.53 | 0.8119 | 0.2312 |
| `pert_wt_fv` | negative, epitope absent | 0.5022 | 0.6463 | **0.0935** | 83.99 | 0.8490 | 0.2482 |
| `vhh_27A05` | negative, epitope mismatch | 0.6174 | 0.6947 | **0.6174** | 87.68 | 0.8947 | 0.1979 |
| `vhh_2Rs15d` | negative, epitope mismatch | 0.0846 | 0.4925 | **0.0846** | 74.19 | 0.7945 | 0.2681 |
| `vhh_cablys3` | negative, irrelevant antigen | 0.0806 | 0.4992 | **0.0806** | 73.69 | 0.7859 | 0.2185 |
| `tras_wt_fv_2chainAg` | representation check | 0.5374 | 0.6773 | **0.0847** | 86.84 | 0.8807 | 0.2472 |

### OpenFold3 FAILS validation on this panel — do not rank designs with it

**Rank inversion.** OpenFold3's top-scoring target is `vhh_27A05` at pair ipTM **0.617**,
which is **4.4x the true binder's 0.141** — and `vhh_27A05` is a negative, a VHH whose real
epitope is HER2 domain I + III. OpenFold3 is also *confident* about this wrong answer
(ab pLDDT 87.7, complex pLDDT 0.895, both the highest in its panel), so no confidence
filter catches it.

`vhh_27A05` is the panel's hardest negative: its measured epitope does clip the very start
of our construct (497-500, 4 residues at the domain III/IV boundary), giving it a genuine
but far-too-small toehold. rf3 handled it correctly and in fact scored it *lowest* of all
seven (0.198). OpenFold3 built a large spurious interface on it instead.

The positive does clear the decoy floor (0.141 vs 0.092), but by 0.049 — against rf3's
0.354 — so even ignoring the inversion there is no usable margin. **Applying the §4 rule
to OpenFold3 would accept a non-binder and reject trastuzumab.**

**Consequence: we are still effectively single-sourced on rf3 for ranking.** OpenFold3 is
wired up, measured and reproducible, and it is worth keeping as a *reporter* (it is cheap,
and a second number is useful for spotting weirdness), but it gets no vote. The next
second-family candidates to try are `af2ig` (once designs have backbones) and
`opendde-abag` **run unguided**; `protenix` is the remaining hosted option.

**One thing OpenFold3 did corroborate: the chain indexing.** Independently of rf3, its
`chains_ptm` for the 4-chain probe is
`{"0": 0.838, "1": 0.070, "2": 0.920, "3": 0.929}` — position 1 again carries the
impossible-for-anything-else value, the two HER2 halves again pair at `[0][1]` = 0.733, and
VH/VL again at `[2][3]` = 0.920. Sorted-chain-id ordering is a property of the platform, not
of one model, and `chain_id_verified = True` for all 7/7 OpenFold3 targets too (§3).

---

## 7. Cost actually spent

| job | model | targets | billable s | USD |
|---|---|---|---|---|
| `9a889c9cca62558160402d2806e54ce1` | rf3 | 7 | 290.1 | **0.020954** |
| `72e7cbe1a6110279ba40be1d62613d24` | openfold3 | 7 | 193.7 | **0.013992** |
| | | | **total** | **$0.034946** |

Spend cap for this task was $15. Actual spend is **$0.0349**, ~0.2% of the cap. Balance
after: $1098.56. Wall-clock: rf3 160 s for all 7 targets, OpenFold3 493 s (it queued behind
other accounts). A 7-target control panel costs about $0.014-0.021 per oracle, so re-running
the full panel alongside every design batch — which §4 requires — is effectively free and
there is no reason ever to skip it.

Re-collecting results from a finished job costs nothing (results are cached server-side), so
`refold.py collect --job <id>` is safe to re-run; both jobs above were collected twice.

---

## 8. Premises in the task brief that turned out to be WRONG

1. **"pertuzumab ... a good oracle should still call it a binder to HER2 overall."**
   Wrong for this antigen. Pertuzumab's epitope, measured from 1S78, is HER2 **domain II**,
   residues 245-315 (23 of 25 heavy-chain contacts in domain II). Our construct is domain IV
   only (489-607), so **pertuzumab's epitope is not in the antigen at all** and it cannot
   bind it. It is an epitope-**absent negative**, not an other-epitope positive, and rf3
   correctly scored it 0.248. This makes it a *better* control than intended — a "right
   protein, wrong epitope" negative is much harder than an irrelevant one — but anyone
   expecting it to score as a binder would have concluded the oracle was broken.

2. **"a known anti-HER2 VHH ... binds trastuzumab's epitope at sub-nM ... the single most
   informative positive control for the VHH route."** No such structure is public. All
   structural anti-HER2 nanobodies in the PDB were checked by computing their epitopes at a
   5 Å contact cutoff:

   | VHH | PDB | measured epitope | domain |
   |---|---|---|---|
   | 2Rs15d | 5MY6 | 109-217 | I |
   | 47D05 | 9MTE | 103-257 | I (14/16 contacts) |
   | 27A05 | 9MTE | 35-500 | I + III; only 497-500 in IV |
   | 29E09 | 9MTX | 35-498 | I + III; only 497-498 in IV |

   **None binds domain IV**, so **the VHH route has no structural positive control.** The
   VHH entries in this panel are necessarily negatives. The sub-nM VHH of PMID:37632647 has
   no PDB entry. If a VHH positive control is needed, it has to come from a non-structural
   source (a published sequence) and will be unvalidated structurally.

3. **"RF3 scores wild-type trastuzumab at only ipTM 0.59."** Only true for the *split*
   antigen. With the correct single-chain construct, rf3 gives trastuzumab WT overall ipTM
   **0.763** and pair ipTM 0.585. Our split-antigen control reproduces 0.6135 ≈ the quoted
   0.59, which strongly suggests that earlier figure was measured on the two-chain
   representation. The brief's conclusion — never use absolute ipTM thresholds, always score
   relative to controls — remains correct and is reinforced, but the 0.59 figure itself is
   an artifact of antigen representation, not a property of trastuzumab or of rf3.

4. **The caplacizumab negative control is not usable from the repo file.**
   `apheris/rfd4/Data/03-vhh/pmvk_caplacizumab.cif.zst` has its CDRs masked to single `X`
   residues (`...SCXWFRQAPGKGRELVAXYPDSVEGRF...`) — it is an RFD4 *design template*, not a
   sequence — and its antigen is phosphomevalonate kinase, not von Willebrand factor.
   Substituted **cAb-Lys3** (PDB 1MEL chain B, anti-hen-egg-lysozyme, HA tag trimmed): a
   real, structurally validated, genuinely irrelevant VHH.

5. **The 581-590 gap is a false dilemma** (one chain with the gap closed *vs* two chains).
   The real sequence is recoverable from 1S78, so the gap can simply be filled correctly.
   See §1 — and the choice matters more than anything else measured here.

6. **Confirmed, not wrong:** `pair_chains_iptm` really is ordered by sorted chain id (§3);
   it is keyed by chain *position*, so the mapping must be derived, and mis-indexing it does
   fail silently. Also confirmed: a per-design shuffled-CDR decoy is mandatory — overall
   ipTM does not separate the scrambled decoy (0.562) from the true binder (0.763) nearly
   well enough, exactly as Smorodina et al. report.

---

## 9. Reproducing

```bash
cd NICCOLO/denovo/oracles
python3 -I make_controls.py                              # -> controls.csv
python3 -I refold.py run --csv controls.csv --model rf3 --seed 1
python3 -I refold.py run --csv controls.csv --model openfold3
python3 -I refold.py table                               # both models, from runs.jsonl
```

To score real designs, write a CSV with the same columns
(`sequence_id, kind, heavy_sequence, light_sequence, antigen_mode, expect, note`), set
`antigen_mode=single`, include a scrambled-CDR twin per design, **and append the control
panel rows to the same CSV** so `P` and `N` of §4 are measured in the same batch.
`light_sequence` empty = VHH. Designs already present in `runs.jsonl` are skipped, so an
interrupted run resumes without paying twice; `--force` overrides.

To route submissions through the central broker when it lands:

```python
import refold
refold.SUBMITTER = my_broker_submit      # (targets, model, name, params) -> (job_id, idem_key)
```

`submit`, `poll` and `collect` are independent; `poll` waits on exactly one job per call.
