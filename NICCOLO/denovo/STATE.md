# STATE — de novo HER2 campaign

**Updated:** 2026-11-14 13:05 Oslo (12:05 UTC)
**Branch:** `niccolo/denovo-loop` · **Deadline:** 18:00 Oslo (17:00 UTC)

## Schedule (revised — user gave final deadline 18:00, not the assumed 19:30)

| Oslo | Phase |
|---|---|
| 13:45 | start |
| 14:30 | P0 artifacts in |
| 15:00 | P1 pilot, lane split chosen |
| 15:00–16:30 | P2 campaign, successive halving |
| 16:30–17:00 | P3 refinement |
| 17:00 | **shortlist frozen** |
| 17:00–17:20 | P4 verify + red team |
| 17:20–17:50 | REPORT.md |
| 18:00 | done |

## Right now

- 6 P0 workers running in parallel (gate, antigen/epitope, scaffolds, oracles, ledger/literature, developability).
- Local RFD4 dry run (1 VHH design, bundled caplacizumab example) launched to de-risk the local lane and measure s/design.
- Orchestrator writing the spec generator and compute broker.

## Done

- Read root AGENTS.md + the 5 required RFD4 docs + challenge doc.
- Branch `niccolo/denovo-loop` created; workspace `NICCOLO/denovo/` laid out.
- **Decoded the RFD4 framework-CIF convention** (this was the main unknown). In `Data/03-fv/framework_adalimumab.cif` the CDRs are *deleted* from the framework and a single skipped residue number marks each insertion point: chain H is `1-22, [23 absent], 24-37, [38 absent], 39-75, [76 absent], 77-87`, and the contig `H1-22, 13-15, H24-37, 9-10, H39-75, 11-22, H77-87` sets each CDR's designed length. Caplacizumab's VHH framework segments the same way (22+14+37+10 = 83 kept residues), so one recipe generalises to any VHH scaffold.
- Confirmed the bundled VHH example (`Data/03-vhh/pmvk_caplacizumab.cif.zst`) is an **annotated** CIF — conditioning rides in extra `_atom_site.condition_*` columns and CDRs appear as `X`/UNK placeholders. Hand-writing that format is error-prone, so the campaign uses the documented plain `input` + `contig` + `conditions` route instead (as in `Data/03-fv/design.json`).

## Blocked / watching

- **Local GPU is NOT idle** — a teammate's `apheris-boltz-v1.3.0-boltz2` container is running (37% util, 4.3 GB). One heavy job at a time still holds; the broker serialises our side only, so local throughput is a measured unknown until the pilot.
- `AMASS_API_KEY` unset, credits unclaimed → literature work is cache-only, zero new searches.

## Facts corrected so far (vs. the briefing premises)

- Deadline is **18:00 Oslo**, not the assumed 17:30/18:30/19:30. Whole schedule compressed by ~1.5 h.
- Novelty reference count is quoted as **525,063** in root AGENTS.md vs 610,109 in the briefing. Being measured directly from the installed bundle; will record the real number.

## Next

1. Spec generator + arm matrix (orchestrator).
2. broker.py with durable `jobs.jsonl` (orchestrator).
3. Pilot one VHH arm per lane; measure cost/time/gate-pass; pick the lane split.

## Log

- 13:45 start; read docs; branch + workspace created.
- 12:50 UTC fanned out 6 P0 workers.
- 12:52 UTC local VHH dry run launched.
- 13:05 Oslo framework-CIF convention decoded; spec route chosen (plain spec, not annotated CIF).
