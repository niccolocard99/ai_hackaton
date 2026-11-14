#!/usr/bin/env python3
"""Literature-informed affinity arm: the measured ADAPT moves on trastuzumab, on top of our CDRH3 designs.

Evidence (LEO/tasks/litreview-optimisation/findings/08-hard-numbers-verified.md, sections 1 and 11):
  * CDRH3 D102 (Kabat H-D98): D98W is a measured 3-fold SPR affinity gain on hu4D5-8 Fab (US-8840896-B2,
    most frequent substitution Trp, 23%); D98F also gains affinity (pmid 31891584) by removing electrostatic
    repulsion with HER2 D586/E598. Our numbering: D102.
  * Light chain D28 (CDR-L1): D28K appears in both ADAPT papers (28750054: D28R/K; 31891584: D28K only), and
    the best triple mutant, D102F + a second single + L-D28R or L-D28K, reached KD 1.9 pM, about 30x better
    than parent. The single-mutation effect of D28K is NOT in the retrieved text, so no size is claimed here.
  * Developability warning from the same papers: D102F + H-D31R was heterogeneous on purification. D31 is
    therefore not used.

What this arm does: (1) a CDRH3 search that forces D102 to F or W and lets the team ranker pick the rest,
under the same gates as the main loop; (2) adds VL D28K (and D28R for the lead designs) to the best CDRH3
designs from both the main portfolio and (1). The ranker was trained on the trastuzumab light chain, so it
says nothing about the D28 change: that effect rests on the literature alone and is unmeasured in our assay.
D28 contacts HER2 in 1N8Z (CDR-L1: ASP28, ASN30, THR31, ALA32), so this changes a contacting residue.
"""
from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, "/home/hackathon/wt/leo-design-loop/NICCOLO/eda")
import designloop as dl                       # noqa: E402
import gpu_stage as gs                        # noqa: E402
import validate_abguard as vab                # noqa: E402
from her2lib import TZ_VL, vh_with_loop       # noqa: E402

OUT = Path(os.environ.get("DL_OUT", HERE / "out"))
LED = Path(os.environ.get("DL_LEDGER", HERE / "ledger"))
dl.OUT = OUT
AA = dl.AA


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def led(*args: str) -> tuple[int, str]:
    p = subprocess.run([sys.executable, str(LED / "ledger.py"), *args], capture_output=True, text=True)
    return p.returncode, (p.stdout + p.stderr).strip()


def d102_arm(ctx) -> pd.DataFrame:
    """CDRH3 designs with D102 forced to F or W; everything else chosen by the ranker."""
    A = [ctx.top_alphabet(p, 7) for p in range(10)]
    A[3] = np.array(sorted(AA.index(c) for c in "FW"), np.int16)
    codes, _, n = dl.product_topk(list(range(10)), A, ctx.G, dl.TZ_CODES, dl.MIN_SUB, 400_000)
    g, d, f = dl.gate(codes, ctx)
    log(f"d102 F/W arm: space {n:,} -> top 400,000 -> novel+scope {f.get('in_scope', 0):,}")
    df = dl.score(g, d, ctx)
    df = df[(df.p_binder_p10 >= dl.P_THR) & (df.support_min >= dl.MIN_SUPPORT)]
    log(f"  after classifier gate and support >= {dl.MIN_SUPPORT}: {len(df):,}")
    sel = dl.select_diverse(df, 30, 2)
    sel["arm"] = "d102FW"
    # the literature-faithful, epitope-safe version: D102 -> F/W plus only the three never-contacting
    # positions (106-108) and the two lowest-contact ones (100, 101, 104); nothing else moves
    P = [1, 2, 3, 5, 7, 8, 9]
    A2 = [np.array([dl.TZ_CODES[p]], np.int16) for p in range(10)]
    for p in P:
        A2[p] = ctx.observed_alphabet(p, 3)
    A2[3] = A[3]
    c2, _, n2 = dl.product_topk(P, A2, ctx.G, dl.TZ_CODES, dl.MIN_SUB, 400_000)
    g2, d2, f2 = dl.gate(c2, ctx)
    df2 = dl.score(g2, d2, ctx)
    df2 = df2[(df2.p_binder_p10 >= dl.P_THR) & (df2.support_min >= dl.MIN_SUPPORT)]
    log(f"  D102F/W conservative arm: space {n2:,} -> in scope {f2.get('in_scope', 0):,} -> classifier {len(df2):,}")
    sel2 = dl.select_diverse(df2, 15, 2)
    sel2["arm"] = "d102FW-cons"
    return pd.concat([sel, sel2], ignore_index=True).drop_duplicates("cdrh3")


def stage1():
    OUT.mkdir(exist_ok=True)
    port = pd.read_csv(OUT / "portfolio.csv")
    ctx = dl.Context(B=100, seed=0)
    flags = ["--method", "prefilter", "--model", "prefilter-niccolo-eda-with-leo-design-loop-v1-d102FW",
             "--seed", "0", "--input", str(HERE / "data/trastuzumab_her2_spr_all_768.csv"),
             "--input", str(HERE / "ref/tz_her2_facs_3classes_524346.csv"),
             "--input", str(HERE / "designloop.py"), "--input", str(HERE / "adapt_arm.py")]
    for k, v in dict(force_D102="F,W", conservative_arm="106-108+100,101,104", top_alphabet=7, keep=400000, min_sub=dl.MIN_SUB, d_max=dl.D_MAX,
                     p_thr=dl.P_THR, k_lcb=dl.K_LCB, min_support=dl.MIN_SUPPORT, n=30).items():
        flags += ["--setting", f"{k}={v}"]
    code, _ = led("check-run", "--ledger-dir", str(LED), *flags)
    rid = f"leo-dl-d102fw-{time.strftime('%H%M')}"
    f102 = OUT / "adapt_d102fw.csv"
    if code == 3 and f102.exists():
        log("ledger: this exact D102 F/W generation is already logged; reusing its output")
        sel = pd.read_csv(f102)
    else:
        led("add-run", "--ledger-dir", str(LED), *flags, "--run-id", rid, "--output", str(f102),
            "--notes", "D102 forced to F or W (measured ADAPT / phage affinity moves), rest chosen by the ranker")
        sel = d102_arm(ctx)
        sel.to_csv(f102, index=False)
        led("set-status", "--ledger-dir", str(LED), "--run-id", rid, "--status", "done" if len(sel) else "empty")
    if len(sel):
        b = sel.iloc[0]
        log(f"  best D102F/W design {b.cdrh3} {b.subs}  pred {b.kd_mean:+.2f}+/-{b.kd_sd:.2f}  p10 {b.p_binder_p10:.2f}")
    sel.to_csv(f102, index=False)
    return sel, port, ctx


def main() -> None:
    sel, port, ctx = stage1()
    cons = sel[sel.arm == "d102FW-cons"].head(4)
    bases = pd.concat([port.head(6)[["cdrh3"]].assign(source="portfolio"),
                       cons[["cdrh3"]].assign(source="d102FW-cons"),
                       sel[sel.arm == "d102FW"].head(4)[["cdrh3"]].assign(source="d102FW")]).drop_duplicates("cdrh3").reset_index(drop=True)
    meta = pd.concat([port, sel], ignore_index=True).drop_duplicates("cdrh3").set_index("cdrh3")
    assert TZ_VL[27] == "D", "VL position 28 is not the CDR-L1 aspartate"
    rows = []
    for r in bases.itertuples():
        for aa in ("K", "R"):
            if aa == "R" and r.Index >= 4:      # D28R only on the first four bases; D28K on all
                continue
            rows.append(dict(id=f"LEO-AD-{r.cdrh3}-LD28{aa}", base=r.cdrh3, source=r.source, light=f"D28{aa}",
                             vh=vh_with_loop(r.cdrh3), vl=TZ_VL[:27] + aa + TZ_VL[28:]))
    flags = ["--method", "mutate", "--model", "adapt-light-chain-D28-literature-v1", "--seed", "0",
             "--input", str(OUT / "portfolio.csv"), "--input", str(OUT / "adapt_d102fw.csv")]
    for k, v in dict(light_mutations="D28K,D28R", bases=",".join(bases.cdrh3),
                     evidence="pmid-28750054+pmid-31891584").items():
        flags += ["--setting", f"{k}={v}"]
    code, _ = led("check-run", "--ledger-dir", str(LED), *flags)
    rid = f"leo-dl-adapt-{time.strftime('%H%M')}"
    if code == 3 and (OUT / "adapt_arm.csv").exists():
        log("ledger: this exact light-chain run is already logged; reusing out/adapt_arm.csv")
        print(pd.read_csv(OUT / "adapt_arm.csv").round(2).to_string(index=False))
        return
    led("add-run", "--ledger-dir", str(LED), *flags, "--run-id", rid, "--output", str(OUT / "adapt_arm.csv"),
        "--notes", "VL D28K/R on the best CDRH3 designs: literature-informed affinity arm, unmeasured in our assay")
    try:
        nm = gs.run_nanomelt(list(bases.cdrh3), f"nm_adapt_{int(time.time())}").set_index("cdrh3")
    except Exception as e:
        log(f"NanoMelt failed: {e}")
        nm = pd.DataFrame(index=bases.cdrh3, columns=["tm_pred", "dtm_vs_wt"])
    rep = vab.validate_rows([dict(sequence_id=r["id"], kind="paired", heavy_sequence=r["vh"], light_sequence=r["vl"],
                                  design_description=f"{r['base']} + VL {r['light']}") for r in rows],
                            "optimisation", shards=4).set_index("sequence_id")
    keep = ["subs", "n_sub", "kd_mean", "kd_sd", "R", "p_binder_p10", "p_beats_tz", "epi_risk", "support_min",
            "nearest_binder", "nearest_binder_kd_nM", "d_meas"]
    out = []
    for r in rows:
        x, m = rep.loc[r["id"]], meta.loc[r["base"]]
        out.append(dict(id=r["id"], base_cdrh3=r["base"], source=r["source"], light_change=r["light"],
                        **{c: m.get(c, np.nan) for c in keep},
                        dtm_vs_wt_vh=nm.loc[r["base"], "dtm_vs_wt"], abg_pass=bool(x.passes),
                        abg_track=x.track_assigned, abg_min_dist=int(x.cdrh3_min_distance),
                        light_contact_changed=True, heavy_sequence=r["vh"], light_sequence=r["vl"]))
    df = pd.DataFrame(out)
    df.to_csv(OUT / "adapt_arm.csv", index=False)
    tbl = pd.DataFrame([dict(sequence_id=r.id, kind="paired", heavy_sequence=r.heavy_sequence,
                             light_sequence=r.light_sequence, design_description=f"{r.base_cdrh3} + VL {r.light_change}",
                             cdrh3=rep.loc[r.id].cdrh3, track_assigned=r.abg_track,
                             base_validation="PASS" if r.abg_pass else "FAIL") for r in df.itertuples()])
    tbl.to_csv(OUT / "adapt_arm_base_validation.csv", index=False)
    led("set-status", "--ledger-dir", str(LED), "--run-id", rid, "--status", "done")
    _, txt = led("add", str(OUT / "adapt_arm_base_validation.csv"), "--ledger-dir", str(LED), "--run-id", rid)
    log("ledger add: " + (txt.splitlines()[-1] if txt else "?"))
    show = ["id", "source", "subs", "kd_mean", "p_binder_p10", "epi_risk", "abg_pass", "abg_min_dist"]
    print(df[show].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
