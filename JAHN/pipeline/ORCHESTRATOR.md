# The main loop: generate → filter → score → orchestrator

`abpipe iterate` runs a fixed loop; `abpipe orchestrate` runs the same stages and **steers them
between rounds**. That fourth component is the difference: a policy that reads the run history and
the machine, adapts the next round's configuration inside a hard allowlist, and stops when the time
budget says one more round will not fit.

```bash
cd JAHN/pipeline
# optimisation track (trastuzumab CDRH3 window, Niccolo's pre-filter)
~/.venvs/abpipe/bin/abpipe orchestrate configs/orchestrated_loop.yaml \
    --rounds 6 --select 60 --out runs/orchestrated --deadline 17:30
# de novo track (ANDERS' VHH generator, PR #11)
~/.venvs/abpipe/bin/abpipe orchestrate configs/orchestrated_denovo.yaml \
    --rounds 4 --select 20 --out runs/orchestrated_denovo
# see what it would decide, run nothing
... orchestrate configs/orchestrated_loop.yaml --dry-run
```

| Component | Where | What it is |
|---|---|---|
| 1. Generate | `generators/scored_evolution.py` (optimisation), `generators/vhh_variants.py` (de novo, wraps `ANDERS/vhh_variants.py`) | CDRH3 windows from the pre-filter prior + CEM feedback from the previous round's elites; or VHH recombine/extend/mutate/graft from a ranked parent list |
| 2. Filter | `filters/` (7 modules) | Cheap first, expensive last; ends at `track_check` = BASE VALIDATION with the team's official scripts |
| 3. Score | `scorers/` | `niccolo_prefilter` for the optimisation track; `pass_through(info__parent_score)` for the de novo lane until a validated scorer exists there |
| 4. Orchestrator | `orchestrate.py` + `policy.py` + `probe.py` | Round driver with elitism + diversity selection, allowlisted config patches, deadline/budget stops, full decision journal |

## What a round does

1. **Probe** the machine: GPU util/mem, heavy Docker containers, disk, load, JapanFold balance
   (`GET /v1/credits`, key from the environment — never written anywhere), AMASS key presence.
2. **Decide**: the policy returns a *patch* — loop knobs and per-module parameters, every one
   allowlisted with hard bounds in `policy.py`. The bounds encode the eligibility floors, so no
   policy can produce a config that breaks them:
   - `min_subs >= 4` — the 70% CDRH3-identity rule binds before the mutation count (AGENTS.md)
   - `min_dist == 3` (novelty) — more than 2 substitutions from every provided window
   - `min_pair_hamming >= 3` — the same rule inside the selected list
   - `identity_max <= 0.68` / `cdrh3_dist_min >= 3` — margins inside the official de novo gate
   Loop knobs: `select`, `min_pair_hamming`, `max_sub_share`, `pool_frac` (how deep into the ranking
   the selection may walk — jahn's first-full-run fix: a short selection means the pool top is too
   mutually similar, and the walk must widen rather than pad the list).
3. **Stop** if the deadline cannot fit another measured round, the disk is < 2 GB, or (rules
   policy) the primary score has plateaued with every exploration lever at its cap.
4. **Run** the round: ordinary `runner.run` stages into `round_NN/`, cached per module; survivors
   are ranked **together with the previous round's selection** (elitism — a round can only improve
   the selected set under the configured scorers), then the top `select` are kept with the
   diversity constraint. The same round-directory layout as `iterate`, so the analysis tooling is
   shared.

## The two policies

- **rules** (default, deterministic): starved funnel → bigger batch (+relaxed `min_p_binder` if
  nothing survived); flooded funnel under time pressure → smaller batch; rich region with time →
  larger `select`; selection itself coming back short → widen `pool_frac`; crowding (low pairwise
  Hamming) → more fresh samples, weaker elite pull; plateau (no score improvement and no new blood
  for `patience` rounds) → explore harder (temperature up, eta down, wider band), and stop as
  *converged* when every lever is capped.
- **agent** (`--policy agent`, default command `claude -p`): the same state digest plus the
  allowlist is sent to an external LLM agent, which answers with a JSON patch **through the same
  validator**. Garbage, timeouts, or out-of-bounds answers fall back to the rules decision, and the
  fallback is recorded. Nothing an agent proposes can touch a knob outside the allowlist or a bound
  inside it.

## The journal (scored in the pitch: "clear reasoning")

Per round, in three places under `--out`:

- `decisions.jsonl` — one machine-readable line per round: probe, budget, patch, reasoning, funnel.
- `DECISIONS.md` — the same, human-readable, one section per round.
- `orchestration.json` — the full state: every round's stats, score medians, `mean_pairwise_hamming`,
  origin mix, wall seconds, and the stopped reason if any.

## Compute and external resources

- The optimisation lane is CPU-only (the pre-filter screens ~1.5M sequences/s; no GPU needed).
- The probe reports the GPU and heavy containers but the loop does not itself start Docker jobs;
  the de novo GPU lane (RFD4-Proteina) is run by NICCOLO/denovo's broker. One H100 for the whole
  team still holds.
- JapanFold: the balance is in every probe snapshot, so an agent policy can weigh "fold the top-k
  remotely" against the local GPU being busy. Refold the shortlist as a plausibility check only —
  measured evidence says fold metrics do not screen or rank CDRH3-window designs
  (NICCOLO/eda/AGENT_GUIDE.md §6).
- AMASS: the probe reports key presence. The key is installed and valid, but the account is out of
  API credits (403 as of 14 Nov 14:45) — the $100 hackathon credits must be claimed by a person at
  platform.amass.tech/settings/credits. Once claimed, use it for prior-art/literature checks on the
  shortlist (PatentCore/BiomedCore), not per round.

## Known gaps (do not read past these)

- The loop optimises the **surrogate's** score. Both pre-filter models are linear one-hot; their
  optimum is nearly reachable in one shot, so the rules mostly adapt the funnel, the exploration
  mix and the stopping time. The loop earns its keep when scorers for thermostability, epitope or
  expression are added — two of the three judged properties are not scored today.
- `pass_through(info__parent_score)` on the de novo lane ranks on an inherited prior, not evidence.
- Two `needs_data` tests in `tests/test_niccolo_loop.py` fail on 14 Nov **before** any of this work
  (`KeyError: 'p_facs_high'` inside `prefilter.py` with the live `~/nick-work/splits`); verified on
  a pristine tree. Everything else (56 tests incl. the 20 new orchestrator tests) passes.
- The de novo cold start needs a real VHH shortlist; `configs/orchestrated_denovo.yaml` points at
  `ANDERS/demo/ranked.csv` as a placeholder.
