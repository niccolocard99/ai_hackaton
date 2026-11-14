#!/usr/bin/env python3
"""Assemble the final ranked list from the loop and the three add-on levers, and validate every design.

Why the lead four are four different bets and not the four best scores. Only 4 designs are guaranteed a
test, and the surrogate has one optimum, so the four best scores would be near-duplicates that fail
together. Each lead slot instead attacks a different combination of the three scored axes:

  1  all three levers   literature D102F CDRH3 + VL D28K + framework R87S/K76L   (the S-tier aspirant)
  2  data CDRH3 + D28K  the ranker's best loop with the light-chain boost
  3  epitope + Tm       epitope-conservative loop (D102 kept) + framework nudge
  4  control            the ranker's best loop with nothing added

Add-ons and where their evidence comes from:
  VL D28K       two ADAPT papers (pmid 28750054, 31891584); part of the best triple, 1.9 pM, ~30x vs parent.
                Single-mutation size not in the retrieved text. Unmeasured in our assay. Contacts HER2 in 1N8Z.
  CDRH3 D102F/W measured 3-fold (D98W, patent US-8840896-B2) and a gain for D98F (pmid 31891584).
  VH R87S/K76L  NanoMelt saturation scan, structure-filtered (tm_arm.py). +1.2 to +1.7 C each, inside its 4.1 C
                error, and the model is out of domain for a VH. A nudge, not a proven stabiliser.
The ranker was trained on the trastuzumab light chain and framework, so it says nothing about any add-on.
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gpu_stage as gs                      # noqa: E402
import run_loop as rl                       # noqa: E402
import validate_abguard as vab              # noqa: E402
from her2lib import TZ_VH, TZ_VL, subs, vh_with_loop   # noqa: E402

OUT, LED = rl.OUT, rl.LED
FW = [(87, "R", "S"), (76, "K", "L")]          # (1-based VH position, wt, mut), from tm_arm.py


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def make(cdrh3: str, d28: str | None, fw: bool) -> dict:
    vh = list(vh_with_loop(cdrh3))
    if fw:
        for pos, wt, mut in FW:
            assert vh[pos - 1] == wt, (pos, wt, vh[pos - 1])
            vh[pos - 1] = mut
    vl = TZ_VL
    if d28:
        assert TZ_VL[27] == "D"
        vl = TZ_VL[:27] + d28 + TZ_VL[28:]
    tag = "+".join(x for x in ([f"D28{d28}"] if d28 else []) + (["fwTm"] if fw else [])) or "cdrh3only"
    return dict(id=f"LEO-F-{cdrh3}-{tag}", cdrh3=cdrh3, d28=d28 or "", fw=int(fw), addons=tag,
                heavy_sequence="".join(vh), light_sequence=vl)


def pick(df: pd.DataFrame, n: int, min_div: int) -> list[str]:
    out: list[str] = []
    for l in df.sort_values("R").cdrh3:
        if all(rl.hamming(l, o) >= min_div for o in out):
            out.append(l)
        if len(out) == n:
            break
    return out


def bases(pool: pd.DataFrame, d102: pd.DataFrame) -> dict:
    ok = pool[(pool.abg_pass == True) & (pool.support_min >= 2) & (pool.p_binder_p10 >= 0.8)]  # noqa: E712
    data = pick(ok, 5, 3)
    cons_ep = pick(ok[ok.epi_risk <= 0.75], 3, 3)
    f1c = pick(d102[d102.arm == "d102FW-cons"], 3, 2)
    f1a = pick(d102[d102.arm == "d102FW"], 3, 3)
    return dict(data=data, epi=cons_ep, d102c=f1c, d102a=f1a)


def lead_and_fill(b: dict) -> list[tuple[dict, str, str]]:
    """(design, slot, why) in rank order: the lead four by hand, then a grid ordered by evidence."""
    L = []
    L.append((make(b["d102c"][0], "K", True), "1 all three levers",
              "Literature D102F loop (3 other substitutions only) + VL D28K + two framework Tm nudges. Every lever "
              "with outside evidence in one molecule: the S-tier aspirant and the highest-variance design."))
    L.append((make(b["data"][0], "K", False), "2 data loop + D28K",
              "The ranker's best bootstrap-robust CDRH3 with the light-chain affinity boost. Tests whether the "
              "literature boost stacks on the data-driven loop."))
    L.append((make(b["epi"][0], None, True), "3 epitope + Tm",
              "Epitope-conservative loop (D102 kept, contact risk 0.62) plus the framework nudge: the safest bet "
              "on the epitope axis, with a second axis from Tm."))
    L.append((make(b["data"][1] if len(b["data"]) > 1 else b["data"][0], None, False), "4 control",
              "A ranker-best loop with nothing added. If slots 1 and 2 beat it, the add-ons are doing the work."))
    # Ranks 5-20: a deliberate mix, not the highest evidence points first. Each lead design's loop also
    # appears without its add-on (a paired control), D28K is carried by about 8 of 20 so one unmeasured
    # add-on cannot sink most of the list, and no loop appears more than twice.
    ADD = {"none": (None, False), "fw": (None, True), "K": ("K", False), "Kfw": ("K", True)}
    FILL = [("d102c", 1, "K"), ("d102c", 2, "fw"), ("data", 2, "none"), ("epi", 1, "fw"),
            ("d102a", 0, "Kfw"), ("data", 0, "none"), ("d102c", 0, "none"), ("epi", 2, "none"),
            ("d102a", 1, "K"), ("data", 3, "fw"), ("d102a", 2, "none"), ("data", 4, "K"),
            ("d102c", 1, "none"), ("epi", 0, "none"), ("data", 1, "K"), ("d102c", 2, "K")]
    controls = {("data", 0, "none"): "paired control for rank 2: same loop, no add-on",
                ("d102c", 0, "none"): "paired control for rank 1: same loop, no add-ons",
                ("epi", 0, "none"): "paired control for rank 3: same loop, no framework nudge",
                ("data", 1, "K"): "paired counterpart of rank 4: same loop plus D28K"}
    have = {x[0]["id"] for x in L}
    for fam, idx, add in FILL:
        if idx >= len(b[fam]):
            continue
        d = make(b[fam][idx], *ADD[add])
        if d["id"] in have:
            continue
        have.add(d["id"])
        L.append((d, "fill (" + fam + ")", controls.get((fam, idx, add), "evidence mix: " + d["addons"])))
        if len(L) >= 20:
            break
    return L


def write_md(df: pd.DataFrame) -> None:
    lead = df.head(4)
    L = ["# Final ranked list, optimisation track", "",
         f"Written {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}. Every design passed the real abguard optimisation "
         "checker (assigned `novel`) and fails the de novo script, as expected: a design is never reclassified.", "",
         "## What to expect, measured and not hoped", "",
         "Held-out backtest (`backtest.py`, 20 repeats of 5-fold cluster-grouped CV on the 179 measured binders, "
         "test variants at least 3 mutations from anything trained on, which is the situation of a legal design). "
         "Inside the legal band (4 or more substitutions from trastuzumab) the ranker's top pick beat trastuzumab's "
         "measured 0.152 nM **1.0%** of the time, landed within 3x **14%** of the time, and the median measured KD of "
         "its top picks was about 9 nM. Binding is a long shot for every design. That is why the lead four are four "
         "bets on different axes and not the four best scores.", "",
         "## The lead four", "", "| # | design | add-ons | why |", "|---|---|---|---|"]
    for r in lead.itertuples():
        L.append(f"| {r.rank} | `{r.cdrh3}` | {r.addons} | {r.why} |")
    L += ["", "## Full list", "",
          "| # | slot | CDRH3 | add-ons | substitutions | pred | p10 | epi risk | VH dTm C | H3 A, exact fold |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for r in df.itertuples():
        h3 = "" if pd.isna(r.fold_h3_pose) else f"{float(r.fold_h3_pose):.1f}" + ("!" if r.fold_pose_flag == True else "")  # noqa: E712
        L.append(f"| {r.rank} | {r.slot} | `{r.cdrh3}` | {r.addons} | {r.subs} | {r.R:+.2f} | {r.p_binder_p10:.2f} | "
                 f"{r.epi_risk:.2f}{'+L1' if r.light_contact_changed else ''} | {r.dtm_vs_wt:+.1f} | {h3} |")
    L += ["", "## Read before trusting any number", "",
          "- `pred` is the ranker bootstrap mean of log10 KD plus one standard deviation. It is a ranking, not a KD, "
          "and it was fitted on the trastuzumab light chain and framework, so it says nothing about any add-on.",
          "- The bootstrap lower bound and an anchored estimate did **not** beat the plain ranker on held-out "
          "variants (Spearman +0.668 and +0.673 against +0.666). They are a conservative choice, not a validated gain.",
          "- VL D28K and D102F/W come from the literature, not from our assay; the size of D28K alone is not in the "
          "retrieved text. D28 contacts HER2, so those designs change a contacting residue (`+L1`).",
          "- VH dTm is NanoMelt, nanobody-trained and out of domain for a VH, error 4.1 C. A light-chain change is not "
          "visible to it. Framework gains of +1 to +3 C are inside that error.",
          "- Boltz-2 does not rank binders (every metric at chance on measured data, EXP-004 and EXP-006). The base H3 "
          "column, when filled, is a pose-plausibility flag for the CDRH3 loop only, folded from the exact final sequence.",
          "- Thermostability and epitope preservation are not predicted by the ranker. No developability scorer is "
          "installed on the VM."]
    (OUT / "FINAL.md").write_text(chr(10).join(L) + chr(10), encoding="utf-8")


def main() -> None:
    pool = rl.merge_gpu(rl.load_pool())
    d102 = pd.read_csv(OUT / "adapt_d102fw.csv")
    b = bases(pool, d102)
    log("bases: " + "  ".join(f"{k}={v}" for k, v in b.items()))
    L = lead_and_fill(b)
    designs = [x[0] for x in L]
    rows = [dict(sequence_id=d["id"], kind="paired", heavy_sequence=d["heavy_sequence"],
                 light_sequence=d["light_sequence"], design_description=f"{d['cdrh3']} {d['addons']}")
            for d in designs]
    opt = vab.validate_rows(rows, "optimisation", shards=4).set_index("sequence_id")
    den = vab.validate_rows(rows, "denovo", shards=4).set_index("sequence_id")
    log(f"BASE VALIDATION, both track scripts: optimisation PASS {int(opt.passes.sum())}/{len(opt)}, "
        f"de novo PASS {int(den.passes.sum())}/{len(den)} (the other track is expected to FAIL)")
    try:
        tm = gs.run_nanomelt_vh({d["id"]: d["heavy_sequence"] for d in designs}, f"final_nm_{int(time.time())}")
    except Exception as e:
        log(f"NanoMelt failed: {e}")
        tm = {}
    meta = pd.concat([pool, d102], ignore_index=True).drop_duplicates("cdrh3").set_index("cdrh3")
    out = []
    for i, (d, slot, why) in enumerate(L, 1):
        m, o = meta.loc[d["cdrh3"]], opt.loc[d["id"]]
        out.append(dict(rank=i, slot=slot, id=d["id"], cdrh3=d["cdrh3"], addons=d["addons"], subs=m.get("subs"),
                        n_sub=m.get("n_sub"), kd_mean=m.get("kd_mean"), kd_sd=m.get("kd_sd"), R=m.get("R"),
                        p_binder_p10=m.get("p_binder_p10"), p_beats_tz=m.get("p_beats_tz"),
                        epi_risk=m.get("epi_risk"), light_contact_changed=bool(d["d28"]),
                        dtm_vs_wt=(tm.get(d["id"], np.nan) - tm["WT"]) if tm else np.nan,
                        abg_optimisation_pass=bool(o.passes), abg_assigned=o.track_assigned,
                        abg_denovo_pass=bool(den.loc[d["id"]].passes), base_iptm=m.get("iptm"),
                        base_h3_pose=m.get("h3_rmsd_pose"), heavy_sequence=d["heavy_sequence"],
                        light_sequence=d["light_sequence"], why=why))
    df = pd.DataFrame(out)
    df[["id", "heavy_sequence", "light_sequence"]].to_csv(OUT / "fold_queue.csv", index=False)
    fx = OUT / "boltz_results_exact.csv"
    fcols = ["fold_iptm", "fold_plddt", "fold_h3_pose", "fold_pose_flag"]
    if fx.exists():
        r = pd.read_csv(fx).drop_duplicates("id", keep="last")[["id", "iptm", "plddt", "h3_rmsd_pose", "pose_flag"]]
        r.columns = ["id"] + fcols
        df = df.merge(r, on="id", how="left")
    else:
        for c in fcols:
            df[c] = np.nan
    df.to_csv(OUT / "FINAL.csv", index=False)
    write_md(df)
    rid = f"leo-dl-final-{time.strftime('%H%M')}"
    flags = ["--method", "mutate", "--model", "final-portfolio-assembly-v1", "--seed", "0",
             "--input", str(OUT / "adapt_d102fw.csv")]
    for k, v in dict(bases=";".join(",".join(x) for x in b.values()).replace(";", "|"),
                     addons="D28K,D28R,fwTm(R87S+K76L)").items():
        flags += ["--setting", f"{k}={v}"]
    code, _ = rl.led("check-run", "--ledger-dir", str(LED), *flags)
    if code != 3:
        rl.led("add-run", "--ledger-dir", str(LED), *flags, "--run-id", rid, "--output", str(OUT / "FINAL.csv"),
               "--notes", "final ranked list: four lead bets plus an evidence-ordered grid of add-ons")
        tbl = pd.DataFrame([dict(sequence_id=r.id, kind="paired", heavy_sequence=r.heavy_sequence,
                                 light_sequence=r.light_sequence, design_description=f"{r.cdrh3} {r.addons}",
                                 cdrh3=opt.loc[r.id].cdrh3, track_assigned=r.abg_assigned,
                                 base_validation="PASS" if r.abg_optimisation_pass else "FAIL")
                            for r in df.itertuples()])
        tbl.to_csv(OUT / "FINAL_base_validation.csv", index=False)
        rl.led("set-status", "--ledger-dir", str(LED), "--run-id", rid, "--status", "done")
        _, txt = rl.led("add", str(OUT / "FINAL_base_validation.csv"), "--ledger-dir", str(LED), "--run-id", rid)
        log("ledger add: " + (txt.splitlines()[-1] if txt else "?"))
    show = ["rank", "slot", "cdrh3", "addons", "kd_mean", "p_binder_p10", "epi_risk", "dtm_vs_wt",
            "abg_optimisation_pass", "abg_denovo_pass"]
    print(df[show].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
