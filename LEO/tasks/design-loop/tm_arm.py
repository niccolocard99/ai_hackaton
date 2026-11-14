#!/usr/bin/env python3
"""Tm bet: framework-stabilised variants of the best CDRH3 designs.

Why. The held-out backtest (backtest.py) puts the chance that a model-picked LEGAL design beats
trastuzumab's measured KD at about 1-2%, and the chance it lands within 3x at about 13%. Binding is a
long shot for every design. Thermostability is the one scored property we can engineer on purpose:
CDRH3 is Tm-inert (EXP-005), framework residues are not.

Position choice is structural, not only NanoMelt's ranking. From the Fv + HER2 cropped 1N8Z file the
RFD4 spec uses, a framework residue is eligible only if
  - the file's own ig_region annotation says framework,
  - no heavy atom is within 6 A of HER2 (the epitope contact surface is untouched),
  - no heavy atom is within 5 A of the VL (the VH/VL interface, which a VH-only model cannot see),
  - no heavy atom is within 5 A of CDRH3 (the designed loop is not perturbed),
  - it is not one of NICCOLO's fragile positions, and not 98, the CDRH3 junction.
Substitutions to C and P are never used. D and E are excluded as well: NanoMelt has a known pull
toward Glu (EXP-005: 11 of the top 12 singles), so those gains are charge artefacts, not stability.

Caveat that stays attached to every number: NanoMelt is trained on nanobodies, so for a VH it is out of
domain. Its predictions pick candidates; they do not prove stability.
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
import gpu_stage as gs                      # noqa: E402
import validate_abguard as vab              # noqa: E402
from her2lib import TZ_VH, TZ_VL, vh_with_loop   # noqa: E402

OUT = Path(os.environ.get("DL_OUT", HERE / "out"))
LED = Path(os.environ.get("DL_LEDGER", HERE / "ledger"))
THERMO = next(p for p in (Path("/home/hackathon/wt/leo-design-loop/NICCOLO/thermo"),
                          Path("/home/hackathon/hackathon-team-3/NICCOLO/thermo"))
              if (p / "nanomelt_sat.csv").exists())
CIF = Path("/home/hackathon/hackathon-team-3/apheris/rfd4/Data/12-her2/her2_fv_1n8z.cif")
FRAGILE = {22, 36, 96, 20, 83, 48, 73, 89}      # NICCOLO/thermo README: do not touch
HER2_CUT, VL_CUT, H3_CUT, MIN_DTM, THREE = 6.0, 5.0, 5.0, 1.0, {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def read_atoms(path: Path) -> pd.DataFrame:
    cols, rows = [], []
    for line in path.read_text().splitlines():
        if line.startswith("_atom_site."):
            cols.append(line.strip().split(".", 1)[1])
        elif line.startswith(("ATOM", "HETATM")):
            f = line.split()
            if len(f) == len(cols):
                rows.append(f)
    d = pd.DataFrame(rows, columns=cols)
    for c in ("Cartn_x", "Cartn_y", "Cartn_z"):
        d[c] = d[c].astype(float)
    d["seq"] = d.label_seq_id.astype(int)
    return d


def position_table() -> pd.DataFrame:
    """For every VH residue: region, and the distance to HER2, to the VL, and to CDRH3."""
    a = read_atoms(CIF)
    xyz = lambda d: d[["Cartn_x", "Cartn_y", "Cartn_z"]].to_numpy()          # noqa: E731
    vh = a[a.label_asym_id == "B"]
    seq = "".join(THREE[r] for r in vh.drop_duplicates("seq").sort_values("seq").label_comp_id)
    assert seq == TZ_VH, "structure chain B is not the trastuzumab VH numbered 1-120"
    her2, vl = xyz(a[a.label_asym_id.isin(["C", "D"])]), xyz(a[a.label_asym_id == "A"])
    h3 = xyz(vh[vh.seq.between(99, 109)])
    rows = []
    for s, g in vh.groupby("seq"):
        p = xyz(g)
        dm = lambda q: float(np.sqrt(((p[:, None, :] - q[None, :, :]) ** 2).sum(2)).min())  # noqa: E731
        rows.append(dict(pos=int(s), wt=THREE[g.label_comp_id.iloc[0]], region=g.ig_region.iloc[0],
                         d_her2=dm(her2), d_vl=dm(vl), d_h3=dm(h3)))
    t = pd.DataFrame(rows)
    t["eligible"] = ((t.region == "framework") & (t.d_her2 > HER2_CUT) & (t.d_vl > VL_CUT)
                     & (t.d_h3 > H3_CUT) & ~t.pos.isin(FRAGILE) & (t.pos != 98))
    return t


def candidates(t: pd.DataFrame) -> pd.DataFrame:
    sat = pd.read_csv(THERMO / "nanomelt_sat.csv").rename(columns={"NanoMelt Tm (C)": "tm"})
    dsat = pd.read_csv(THERMO / "data_saturation.csv")
    wt = float(sat.loc[sat.ID == "WT", "tm"].iloc[0])
    m = dsat[dsat.name != "WT"].merge(sat[["ID", "tm"]], left_on="name", right_on="ID")
    m["dtm"] = m.tm - wt
    m = m.merge(t[["pos", "region", "d_her2", "d_vl", "d_h3", "eligible"]], on="pos")
    ok = m[m.eligible & ~m.mut.isin(list("CPED")) & (m.dtm >= MIN_DTM)]
    return ok.sort_values("dtm", ascending=False).drop_duplicates("pos").reset_index(drop=True)


def led(*args: str) -> tuple[int, str]:
    p = subprocess.run([sys.executable, str(LED / "ledger.py"), *args], capture_output=True, text=True)
    return p.returncode, (p.stdout + p.stderr).strip()


def build(bases: list[str], muts: pd.DataFrame) -> list[dict]:
    sets = {r.name: [r] for r in muts.itertuples()}
    if len(muts) >= 2:
        sets["+".join(r.name for r in muts.head(2).itertuples())] = list(muts.head(2).itertuples())
    rows = []
    for loop in bases:
        for label, ms in sets.items():
            vh = list(vh_with_loop(loop))
            for m in ms:
                assert vh[m.pos - 1] == m.wt, (loop, m.name)
                vh[m.pos - 1] = m.mut
            rows.append(dict(id=f"LEO-FW-{loop}-{label}", base=loop, fw=label, vh="".join(vh)))
    return rows


def main() -> None:
    port = pd.read_csv(OUT / "portfolio.csv")
    bases = port.cdrh3.head(4).tolist()
    t = position_table()
    cand = candidates(t)
    muts = cand.head(2)
    log(f"eligible framework positions {int(t.eligible.sum())}/84; candidates {list(cand.name)}; "
        f"using {list(muts.name)} on bases {bases}")
    inputs = [OUT / "portfolio.csv", THERMO / "nanomelt_sat.csv", THERMO / "data_saturation.csv", CIF]
    flags = ["--method", "mutate", "--model", "framework-picker-nanomelt-sat-1n8z-v1", "--seed", "0"]
    for f in inputs:
        flags += ["--input", str(f)]
    for k, v in dict(min_dtm=MIN_DTM, her2_cut=HER2_CUT, vl_cut=VL_CUT, h3_cut=H3_CUT,
                     fragile=",".join(map(str, sorted(FRAGILE))), excluded_subs="CPED",
                     bases=",".join(bases), muts=",".join(muts.name)).items():
        flags += ["--setting", f"{k}={v}"]
    code, _ = led("check-run", "--ledger-dir", str(LED), *flags)
    rid = f"leo-dl-tmarm-{time.strftime('%H%M')}"
    if code == 3:
        log("ledger: this exact framework run is already logged; reusing out/tm_arm.csv, not regenerating")
        print(pd.read_csv(OUT / "tm_arm.csv").round(2).to_string(index=False))
        return
    led("add-run", "--ledger-dir", str(LED), *flags, "--run-id", rid, "--output", str(OUT / "tm_arm.csv"),
        "--notes", "framework-stabilised variants of the top portfolio CDRH3 designs: Tm bet")
    rows = build(bases, muts)
    seqs = {r["id"]: r["vh"] for r in rows} | {f"BASE-{l}": vh_with_loop(l) for l in bases}
    tm = gs.run_nanomelt_vh(seqs, f"tmarm_nm_{int(time.time())}")
    rep = vab.validate_rows([dict(sequence_id=r["id"], kind="paired", heavy_sequence=r["vh"],
                                  light_sequence=TZ_VL, design_description=f"{r['base']} + framework {r['fw']}")
                             for r in rows], "optimisation", shards=4)
    rep = rep.set_index("sequence_id")
    out = []
    for r in rows:
        x = rep.loc[r["id"]]
        out.append(dict(id=r["id"], base_cdrh3=r["base"], framework=r["fw"], tm_pred=tm[r["id"]],
                        dtm_vs_wt=tm[r["id"]] - tm["WT"], dtm_vs_base=tm[r["id"]] - tm[f"BASE-{r['base']}"],
                        abg_pass=bool(x.passes), abg_track=x.track_assigned, abg_min_dist=int(x.cdrh3_min_distance),
                        heavy_sequence=r["vh"], light_sequence=TZ_VL))
    df = pd.DataFrame(out)
    df.to_csv(OUT / "tm_arm.csv", index=False)
    tbl = pd.DataFrame([dict(sequence_id=r.id, kind="paired", heavy_sequence=r.heavy_sequence,
                             light_sequence=r.light_sequence, design_description=f"{r.base_cdrh3} + {r.framework}",
                             cdrh3=rep.loc[r.id].cdrh3, track_assigned=r.abg_track,
                             base_validation="PASS" if r.abg_pass else "FAIL") for r in df.itertuples()])
    tbl.to_csv(OUT / "tm_arm_base_validation.csv", index=False)
    led("set-status", "--ledger-dir", str(LED), "--run-id", rid, "--status", "done")
    _, txt = led("add", str(OUT / "tm_arm_base_validation.csv"), "--ledger-dir", str(LED), "--run-id", rid)
    log("ledger add: " + (txt.splitlines()[-1] if txt else "?"))
    print(f"\ncontrol (trastuzumab VH) Tm {tm['WT']:.2f} C; NanoMelt error bar 4.1 C, so every dTm below is "
          "inside the error")
    print(df.drop(columns=["heavy_sequence", "light_sequence"]).round(2).to_string(index=False))


if __name__ == "__main__":
    main()
