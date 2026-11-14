"""The orchestrator's brain: turn round results and a resource snapshot into the next round's config.

A policy never edits the config directly. It returns a PATCH: loop-level knobs plus per-module
parameters, each one allowlisted with hard bounds in this file. The bounds encode the challenge's
eligibility floors, so no policy — rule-based or LLM — can produce a config that breaks them:

  min_subs >= 4        the 70% CDRH3-identity rule binds before the mutation count (AGENTS.md)
  min_dist == 3        more than 2 substitutions from every provided window
  min_pair_hamming >= 3  diversity floor, same rule applied inside the selected list

`AgentPolicy` runs an external agent (default `claude -p`) against a state digest and applies its
JSON answer through the same allowlist as the rules; anything invalid falls back to `rules_decision`.
"""
from __future__ import annotations

import json
import shlex
import subprocess
from dataclasses import dataclass, field

LOOP_BOUNDS = {
    "select": (int, 10, 200),
    "min_pair_hamming": (int, 3, 5),          # floor 3: the >2-mutations-from-provided rule
    "max_sub_share": (float, 0.05, 0.5),
    "pool_frac": (float, 0.02, 0.5),          # how deep into the ranking selection may walk
}

MODULE_BOUNDS = {
    "scored_evolution": {
        "n": (int, 100, 200_000),
        "min_subs": (int, 4, 10),              # floor 4: the CDRH3 identity rule at 13-aa parent length
        "max_subs": (int, 1, 10),
        "eta": (float, 0.1, 1.0),
        "temperature": (float, 0.25, 3.0),
        "n_elite": (int, 5, 500),
        "frac_sample": (float, 0.0, 1.0),
        "frac_mutate": (float, 0.0, 1.0),
        "frac_cross": (float, 0.0, 1.0),
        "w_clf": (float, 0.0, 3.0),
        "w_rank": (float, 0.0, 3.0),
        "w_facs": (float, 0.0, 3.0),
        "min_dist_provided": (int, 3, 3),      # locked: eligibility
    },
    "vhh_variants": {
        "n": (int, 4, 200),
        "parents_top": (int, 2, 24),
        "per_strategy": (int, 2, 100),
        "identity_max": (float, 0.40, 0.68),    # margin, always strictly under the 0.70 cap
        "cdrh3_dist_min": (int, 3, 8),          # margin, always over the official floor of 2
    },
    "prefilter_gate": {"min_p_binder": (float, 0.5, 1.0)},
    "window_novelty": {"min_dist": (int, 3, 3)},   # locked: eligibility
    "new_liabilities": {"max_new": (int, 0, 2)},
    "dedupe": {"min_cdrh3_distance": (int, 0, 3)},
}

OUTPUT_CONTRACT = """Reply with ONLY a JSON object, no prose before or after:
{"reasoning": ["one short sentence per change, in your own words"],
 "patch": {"loop": {<subset of the loop knobs>}, "modules": {"<module name>": {<param>: <value>}}},
 "stop": false, "stop_reason": ""}"""


class PolicyError(ValueError):
    """The patch is not allowlisted; the message says which key and what is allowed."""


@dataclass
class Decision:
    patch: dict = field(default_factory=lambda: {"loop": {}, "modules": {}})
    reasoning: list = field(default_factory=list)
    stop: bool = False
    stop_reason: str = ""
    source: str = "rules"


def _check(kind, name, bounds, where):
    if not isinstance(name, str) or name not in bounds:
        raise PolicyError(f"{where}: unknown {kind} {name!r}; allowed {sorted(bounds)}")
    for param, (typ, lo, hi) in bounds[name].items():
        yield param, typ, lo, hi


def validate_patch(patch: dict) -> None:
    """Shape + bounds check. Raises PolicyError; never mutates the input."""
    if not isinstance(patch, dict):
        raise PolicyError("the patch must be a JSON object")
    unknown = set(patch) - {"loop", "modules"}
    if unknown:
        raise PolicyError(f"patch keys {sorted(unknown)}; only 'loop' and 'modules' are allowed")
    loop = patch.get("loop") or {}
    if not isinstance(loop, dict):
        raise PolicyError("patch.loop must be an object")
    for k, v in loop.items():
        typ, lo, hi = LOOP_BOUNDS.get(k, (None, None, None))
        if typ is None:
            raise PolicyError(f"loop knob {k!r} is not tunable; allowed {sorted(LOOP_BOUNDS)}")
        if not isinstance(v, typ) or isinstance(v, bool) or not lo <= v <= hi:
            raise PolicyError(f"loop {k} = {v!r}; must be {typ.__name__} in [{lo}, {hi}]")
    mods = patch.get("modules") or {}
    if not isinstance(mods, dict):
        raise PolicyError("patch.modules must be an object")
    for name, params in mods.items():
        if not isinstance(params, dict) or not params:
            raise PolicyError(f"patch.modules[{name!r}] must be a non-empty object of param: value")
        for param, typ, lo, hi in _check("module", name, MODULE_BOUNDS, f"module {name!r}"):
            if param not in params:
                continue
            v = params[param]
            if not isinstance(v, typ) or isinstance(v, bool) or not lo <= v <= hi:
                raise PolicyError(f"{name}.{param} = {v!r}; must be {typ.__name__} in [{lo}, {hi}]")


def apply_loop_patch(loop_params: dict, patch: dict) -> dict:
    validate_patch(patch)
    out = dict(loop_params)
    out.update(patch.get("loop") or {})
    return out


def apply_module_patch(raw: dict, patch: dict) -> dict:
    """A deep-copied config `raw` with the patch's module params applied to every spec with that name.

    'n' is a generator-spec key (not a param); everything else lands in the module's params.
    Cross-field rules are re-checked against the MERGED values, so a patch cannot leave the config
    self-inconsistent even when it only touches one side of min_subs/max_subs.
    """
    import copy
    validate_patch(patch)
    mods = patch.get("modules") or {}
    if not mods:
        return copy.deepcopy(raw)
    out = copy.deepcopy(raw)
    present = {name for specs in (out.get("generators") or [], out.get("filters") or [], out.get("scorers") or [])
               for name in ({s["name"] if isinstance(s, dict) else s for s in specs} or set())}
    missing = set(mods) - present
    if missing:
        raise PolicyError(f"patch names modules the config does not have: {sorted(missing)}")
    for section, kind in (("generators", "generator"), ("filters", "filter"), ("scorers", "scorer")):
        for spec in out.get(section) or []:
            if isinstance(spec, str):
                spec = {"name": spec}
            name = spec["name"]
            if name not in mods:
                continue
            spec.setdefault("params", {})
            for param, value in mods[name].items():
                if kind == "generator" and param == "n":
                    spec["n"] = value
                else:
                    spec["params"][param] = value
    se = next((s.get("params", {}) for s in out.get("generators") or []
               if (s if isinstance(s, dict) else {"name": s}).get("name") == "scored_evolution"), None)
    if se is not None:
        if se.get("min_subs", 4) > se.get("max_subs", 8):
            raise PolicyError(f"min_subs {se.get('min_subs')} > max_subs {se.get('max_subs')}")
        fr = [se.get(k, 0.4) for k in ("frac_sample", "frac_mutate", "frac_cross")]
        if sum(fr) <= 0:
            raise PolicyError("frac_sample + frac_mutate + frac_cross must stay positive")
    return out


# --------------------------------------------------------------------------- rules

def _clamp(value, lo, hi):
    return max(lo, min(hi, value))


def rules_decision(state: dict) -> Decision:
    """Deterministic adaptation. Order matters: safety stops first, then the funnel, then learning."""
    reasoning: list = []
    loop: dict = {}
    mods: dict = {}
    stop, why = False, ""
    probe, eff = state["probe"], state["effective"]
    hist = state.get("rounds") or []
    rem, est = state.get("deadline_remaining_s"), state.get("est_round_seconds")
    lp = state["loop"]

    if rem is not None:
        if rem <= 60:
            return Decision(reasoning=[f"deadline is in {rem / 60:.0f} min"], stop=True, stop_reason="deadline")
        if est is not None and rem < est + 60:
            return Decision(reasoning=[f"only {rem / 60:.0f} min left but a round costs about {est / 60:.0f} min"],
                            stop=True, stop_reason="deadline")
    free = probe.get("disk_free_gb")
    if free is not None and free < 2.0:
        return Decision(reasoning=[f"only {free:.1f} GB free"], stop=True, stop_reason="disk")
    if probe.get("heavy_job_running"):
        reasoning.append("a heavy GPU job is running; the loop is CPU-only and continues, GPU advice is unaffected")

    if not hist:
        reasoning.append("cold start: run round 0 on the priors, adapt from its funnel")
        return Decision({"loop": {}, "modules": {}}, reasoning, source="rules")

    last = hist[-1]
    select = lp["select"]
    se = eff.get("scored_evolution") or {}
    pg = (eff.get("prefilter_gate") or {}).get("min_p_binder")

    _GEN_N = {g: MODULE_BOUNDS[g]["n"] for g in ("scored_evolution", "vhh_variants")}
    if last["survivors"] < select:
        note = f"funnel is starved: {last['survivors']} survivors for {select} slots"
        for gname, ncap in _GEN_N.items():
            g = eff.get(gname) or {}
            if "n" not in g:
                continue
            n2 = _clamp(int(g["n"]) * 3 // 2, *ncap[1:])
            if n2 > g["n"]:
                mods.setdefault(gname, {})["n"] = n2
                note += f"; {gname} batch {g['n']} -> {n2}"
        if last["survivors"] == 0 and pg is not None and pg > 0.5:
            pg2 = round(pg - 0.05, 2)
            mods.setdefault("prefilter_gate", {})["min_p_binder"] = pg2
            note += f"; nothing survived, relax the gate {pg:g} -> {pg2:g}"
        reasoning.append(note)
    elif last["survivors"] > 6 * select and rem is not None and est is not None and rem < 6 * est:
        for gname, ncap in _GEN_N.items():
            g = eff.get(gname) or {}
            if "n" in g:
                n2 = max(ncap[1], int(g["n"]) * 3 // 4)
                if n2 < g["n"]:
                    mods.setdefault(gname, {})["n"] = n2
                    reasoning.append(f"funnel is flooded ({last['survivors']} survivors) and time is short; "
                                     f"{gname} batch {g['n']} -> {n2}")
    elif last["survivors"] > 4 * select and (rem is None or est is None or rem > 3 * est):
        sel2 = _clamp(select * 2, *LOOP_BOUNDS["select"][1:])
        if sel2 > select:
            loop["select"] = sel2
            reasoning.append(f"rich region ({last['survivors']} survivors) with time to spare; select {select} -> {sel2}")

    sel_n = last.get("selected")
    if sel_n is not None and sel_n < select and last["survivors"] >= select:
        pf = round(_clamp(float(lp.get("pool_frac", 0.10)) * 1.5, *LOOP_BOUNDS["pool_frac"][1:]), 4)
        if pf > lp.get("pool_frac", 0.10):
            loop["pool_frac"] = pf
            reasoning.append(f"selection itself came back short ({last['selected']} of {select}): the top of the "
                             f"pool is too mutually similar; widen the walk to the best {pf:.0%} of the ranking")

    mph = last.get("mean_pairwise_hamming")
    if mph is not None and mph <= lp["min_pair_hamming"] + 0.5:
        bumped = False
        if "frac_sample" in se:
            fs = _clamp(float(se.get("frac_sample", 0.4)) + 0.15, 0.05, 1.0)
            if fs > se.get("frac_sample", 0.4):
                mods.setdefault("scored_evolution", {})["frac_sample"] = fs
                bumped = True
        if "eta" in se:
            eta = _clamp(float(se.get("eta", 0.5)) - 0.1, *MODULE_BOUNDS["scored_evolution"]["eta"][1:])
            if eta < se.get("eta", 0.5):
                mods.setdefault("scored_evolution", {})["eta"] = eta
                bumped = True
        if bumped:
            reasoning.append(f"selected windows are crowding each other (mean pairwise Hamming {mph:g}); "
                             f"more fresh samples, weaker pull towards the elites")

    patience = int(state.get("patience", 2))
    col, direction = state.get("primary"), state.get("primary_direction")
    if col and len(hist) >= patience + 1:
        window = [h for h in hist[-(patience + 1):]
                  if (h.get("selected_medians") or {}).get(col) is not None]
        if len(window) == patience + 1 and (window[-1].get("selected_medians") or {}).get(col) is not None:
            med = [h["selected_medians"][col] for h in window]
        better = (lambda a, b: a < b) if direction == "lower" else (lambda a, b: a > b)
        improved = better(med[-1], med[0])
        new_blood = any(h.get("new_in_selected", 0) > 0 for h in window[1:])
        if not improved and not new_blood:
            caps, bumps = [], {}
            if "temperature" in se:
                t2 = _clamp(float(se.get("temperature", 1.0)) + 0.25, *MODULE_BOUNDS["scored_evolution"]["temperature"][1:])
                if t2 > se.get("temperature", 1.0):
                    bumps["temperature"] = t2
                else:
                    caps.append("temperature")
            if "eta" in se:
                e2 = _clamp(float(se.get("eta", 0.5)) - 0.15, *MODULE_BOUNDS["scored_evolution"]["eta"][1:])
                if e2 < se.get("eta", 0.5):
                    bumps["eta"] = e2
                else:
                    caps.append("eta")
            if "max_subs" in se:
                m2 = _clamp(int(se.get("max_subs", 8)) + 1, int(se.get("min_subs", 4)),
                            MODULE_BOUNDS["scored_evolution"]["max_subs"][2])
                if m2 > se.get("max_subs", 8):
                    bumps["max_subs"] = m2
                else:
                    caps.append("max_subs")
            if bumps:
                mods.setdefault("scored_evolution", {}).update(bumps)
                reasoning.append(f"plateau: {col} median has not improved and no new window entered the "
                                 f"selection for {patience} rounds; explore harder ({bumps})")
            elif not mods:
                stop, why = True, "converged"
                reasoning.append(f"plateau and every exploration lever is at its cap {caps}; the surrogate's "
                                 f"optimum is reached. Add a scorer that sees a different property, or stop.")
    return Decision({"loop": loop, "modules": mods}, reasoning, stop, why, source="rules")


# --------------------------------------------------------------------------- agent

def _allowlist_doc() -> str:
    def fmt(bounds):
        return "\n".join(f"    {k}: {typ.__name__} in [{lo}, {hi}]" for k, (typ, lo, hi) in bounds.items())
    return (f"loop knobs:\n{fmt(LOOP_BOUNDS)}\nmodule params "
            + "(generators may also take 'n' at spec level):\n"
            + "\n".join(f"  {name}:\n{fmt(p)}" for name, p in MODULE_BOUNDS.items()))


def build_prompt(state: dict) -> str:
    digest = {k: v for k, v in state.items() if k != "prompt_extra"}
    return f"""You are the orchestrator of an antibody-design loop (generate -> filter -> score -> select).
Decide the NEXT round's configuration adjustments from the run history and the machine state.

Rules you cannot change and must not propose values for: min_subs stays >= 4, novelty distances stay
at 3 (the >2-mutations-from-provided rule), min_pair_hamming stays >= 3. The allowlist below is the
complete set of knobs you may touch; anything outside it is rejected and the loop falls back to the
built-in rules.

{_allowlist_doc()}

State (JSON):
{json.dumps(digest, indent=2, default=str)}

Guidance: adapt the funnel (batch size, gates) when survivors are far from the select target; increase
exploration (temperature up, eta down, wider substitution band) on score plateaus; stop when the
surrogate has converged and no lever is left, or when the time budget cannot fit another round.
Every value must be absolute, not a delta. Keep reasoning short and factual.

{OUTPUT_CONTRACT}"""


def _extract_json(text: str) -> dict:
    i = text.find("{")
    if i < 0:
        raise PolicyError("the agent answer contains no JSON object")
    obj, _ = json.JSONDecoder().raw_decode(text[i:])
    return obj


class AgentPolicy:
    """Ask an external agent; anything it proposes still goes through the same allowlist."""

    def __init__(self, agent_cmd: str = "claude -p", timeout_s: int = 600):
        self.agent_cmd, self.timeout_s = agent_cmd, timeout_s

    def __call__(self, state: dict) -> Decision:
        prompt = build_prompt(state)
        failure = None
        try:
            proc = subprocess.run(shlex.split(self.agent_cmd), input=prompt, capture_output=True,
                                  text=True, timeout=self.timeout_s)
            if proc.returncode != 0:
                failure = f"exit {proc.returncode}: {(proc.stderr or '')[-300:]}"
            else:
                answer = _extract_json(proc.stdout)
                patch = {"loop": answer.get("patch", {}).get("loop") or {},
                         "modules": answer.get("patch", {}).get("modules") or {}}
                validate_patch(patch)
                reasoning = list(answer.get("reasoning") or [])
                if not reasoning:
                    raise PolicyError("the agent gave no reasoning")
                stop = bool(answer.get("stop"))
                stop_reason = str(answer.get("stop_reason") or ("stopped by the agent policy" if stop else ""))
                return Decision(patch, reasoning, stop, stop_reason, source="agent")
        except PolicyError as e:
            failure = f"rejected patch: {e}"
        except Exception as e:                                            # timeout, missing binary, bad JSON
            failure = f"{type(e).__name__}: {e}"
        fallback = rules_decision(state)
        note = f"agent policy failed ({failure}); fell back to the rules policy"
        return Decision(fallback.patch, [note, *fallback.reasoning], fallback.stop, fallback.stop_reason,
                        source="rules(agent-failed)")


def make_policy(kind: str, agent_cmd: str | None = None):
    if kind == "rules":
        return rules_decision
    if kind == "agent":
        return AgentPolicy(agent_cmd or "claude -p")
    raise ValueError(f"unknown policy {kind!r}; use 'rules' or 'agent'")
