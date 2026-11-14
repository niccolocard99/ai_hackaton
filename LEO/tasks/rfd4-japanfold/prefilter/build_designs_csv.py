"""Build the BASE VALIDATION table for every candidate we hold.

AGENTS.md: no design is a candidate until it passes base validation, and the
table needs columns sequence_id, kind, heavy_sequence, light_sequence,
design_description.

Candidates come from two independent sources that disagree about which CDRH3
positions to touch, and both are included on purpose so abguard judges them on
the same footing:
  LEO-*   geometry-driven, 4 substitutions at the contact-free positions
  RFD4-*  RFD4-Proteina via JapanFold, 4-9 substitutions, mutating W99/Y105

The parent VH/VL are taken from the SPR table's own trastuzumab row rather than
retyped, so the framework is exactly what the organisers' data says it is.
"""
import sys
from pathlib import Path

import pandas as pd

TZ = "WGGDGFYAMD"
SPR = Path("data/trastuzumab_her2_spr_all_768.csv")
CDRH3_SLICE = slice(98, 108)  # 0-based: VH positions 99-108

# geometry-driven candidates, as 10-mers (Y109 is outside the window)
LEO = {
    "WGGEGFYSLE": "A1/A2 geometry: D102E A106S M107L D108E",
    "WGAEGFYSLD": "C1/C2 geometry: G101A D102E A106S M107L",
    "WGAEGFYSME": "D1 geometry: G101A D102E A106S D108E",
    "WGGEGFYSIE": "E1 geometry: D102E A106S M107I D108E",
}


def main():
    spr = pd.read_csv(SPR)
    par = spr[spr.cdrh3 == TZ]
    if par.empty:
        sys.exit("trastuzumab row not found in the SPR table")
    vh, vl = par.iloc[0].vh_sequence, par.iloc[0].vl_sequence
    print(f"parent VH {len(vh)} aa, VL {len(vl)} aa")
    assert vh[CDRH3_SLICE] == TZ, f"CDRH3 not at 99-108: found {vh[CDRH3_SLICE]}"
    print(f"CDRH3 confirmed at VH[99:108] = {vh[CDRH3_SLICE]}")

    rows = [{
        "sequence_id": "TRASTUZUMAB-PARENT",
        "kind": "paired",
        "heavy_sequence": vh,
        "light_sequence": vl,
        "design_description": "unmodified parent, control (expected FAIL: 0 mutations)",
    }]

    def add(prefix, cdrh3, desc, i):
        rows.append({
            "sequence_id": f"{prefix}-{i:02d}",
            "kind": "paired",
            "heavy_sequence": vh[:98] + cdrh3 + vh[108:],
            "light_sequence": vl,
            "design_description": desc,
        })

    for i, (s, d) in enumerate(LEO.items(), 1):
        add("LEO", s, d, i)

    rfd4 = pd.read_csv("rfd4_designs_gated.csv")
    rfd4 = rfd4[rfd4.legal].sort_values("pred_kd_nM")
    for i, r in enumerate(rfd4.itertuples(), 1):
        subs = " ".join(f"{a}{99+k}{b}" for k, (a, b) in enumerate(zip(TZ, r.cdrh3)) if a != b)
        add("RFD4", r.cdrh3, f"RFD4-Proteina/JapanFold n={r.n_designs} seeds: {subs}", i)

    df = pd.DataFrame(rows)
    df["cdrh3"] = [h[CDRH3_SLICE] for h in df.heavy_sequence]
    bad = df[df.heavy_sequence.str.len() != len(vh)]
    assert bad.empty, f"length changed for {list(bad.sequence_id)}"
    assert df.sequence_id.is_unique, "duplicate sequence_id"
    print(f"\n{len(df)} rows; all heavy chains {len(vh)} aa; VH+VL = {len(vh)+len(vl)} aa")
    df.drop(columns="cdrh3").to_csv("designs.csv", index=False)
    print("wrote designs.csv")
    print(df[["sequence_id", "cdrh3", "design_description"]].to_string(index=False))


if __name__ == "__main__":
    main()
