"""Harvest RFD4-Proteina CDRH3 designs from JapanFold, then gate them.

Run from ~/leo-work/prefilter (needs data/ for the novelty set and the ranker).
Usage: python -I harvest_rfd4.py <designs-root> [<designs-root> ...]

Order matters. The novelty and identity gates are HARD RULES and are applied
BEFORE any score is reported, because AGENTS.md is explicit that the
conservative weights' best-scoring designs are exactly the ones that get
rejected -- filter for novelty first, never after.

The framework check is not a formality. The contig says B1-98 and B109-120 stay
fixed; this verifies RFD4 actually did that, rather than trusting it.
"""
import io
import json
import sys
import zipfile
from collections import OrderedDict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path.cwd()))
from prefilter import TZ, PreFilter, load  # noqa: E402

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
PARENT_CIF = Path.home() / "leo-work/japanfold/her2_fv_1n8z.cif"
FACS_ZIP = Path(
    "/home/hackathon/hackathon-team-3/apheris/hackathon-data/her2-trastuzumab"
    "/tz_her2_facs_3classes_524346.csv.zip"
)
CDRH3_RANGE = range(99, 109)  # B99-108, the designed window


def chain_seq(cif_path, chain):
    """{residue number: one-letter} for one chain, by mmCIF column name."""
    cols, out = [], OrderedDict()
    with open(cif_path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("_atom_site."):
                cols.append(s.split(".", 1)[1])
                continue
            if not s.startswith(("ATOM", "HETATM")) or not cols:
                continue
            f = s.split()
            ia = cols.index("auth_asym_id") if "auth_asym_id" in cols else cols.index("label_asym_id")
            ir = cols.index("auth_seq_id") if "auth_seq_id" in cols else cols.index("label_seq_id")
            ic = cols.index("label_comp_id")
            if f[ia] != chain:
                continue
            out[int(f[ir])] = THREE_TO_ONE.get(f[ic].upper(), "X")
    return out


def known_cdrh3():
    seen = set(pd.read_csv("data/trastuzumab_her2_spr_all_768.csv").cdrh3.dropna())
    with zipfile.ZipFile(FACS_ZIP) as z:
        name = next(n for n in z.namelist() if n.endswith(".csv"))
        with z.open(name) as fh:
            for chunk in pd.read_csv(io.TextIOWrapper(fh, "utf-8"), chunksize=200_000):
                col = "seq" if "seq" in chunk.columns else chunk.columns[0]
                seen |= set(chunk[col].dropna())
    return seen


def main():
    roots = [Path(p) for p in sys.argv[1:]] or [Path.home() / "leo-work/japanfold"]
    parent = chain_seq(PARENT_CIF, "B")
    parent_fw = "".join(parent[i] for i in sorted(parent) if i not in CDRH3_RANGE)
    print(f"parent VH {len(parent)} res; CDRH3 = {''.join(parent[i] for i in CDRH3_RANGE)}")
    assert "".join(parent[i] for i in CDRH3_RANGE) == TZ, "parent CDRH3 != TZ"

    rows, broken = [], 0
    for root in roots:
        for cif in sorted(root.rglob("design_*_generated.cif")):
            b = chain_seq(cif, "B")
            if not all(i in b for i in CDRH3_RANGE):
                broken += 1
                continue
            fw = "".join(b[i] for i in sorted(b) if i not in CDRH3_RANGE)
            cdrh3 = "".join(b[i] for i in CDRH3_RANGE)
            rows.append({
                "cdrh3": cdrh3,
                "job": cif.parent.name[:8],
                "design": cif.stem.replace("_generated", ""),
                "framework_intact": fw == parent_fw,
            })
    df = pd.DataFrame(rows)
    print(f"\nharvested {len(df)} designs from {len(roots)} run dir(s); {broken} unparseable")
    if df.empty:
        return
    bad_fw = int((~df.framework_intact).sum())
    print(f"framework held fixed in {len(df)-bad_fw}/{len(df)} designs"
          + ("  *** RFD4 CHANGED THE FRAMEWORK, contig not honoured ***" if bad_fw else "  (contig honoured)"))
    df = df[df.framework_intact]

    u = df.groupby("cdrh3").agg(n_designs=("design", "size"),
                                jobs=("job", lambda s: ",".join(sorted(set(s))))).reset_index()
    u["n_mut"] = [sum(a != b for a, b in zip(s, TZ)) for s in u.cdrh3]
    u["identity"] = (10 - u.n_mut) / 10 * 100
    u["has_cys"] = u.cdrh3.str.contains("C")
    u["has_asn"] = u.cdrh3.str.contains("N")
    print(f"{len(u)} distinct CDRH3s among {len(df)} designs "
          f"({len(df)-len(u)} duplicates -> diversity {len(u)/len(df):.0%})")

    print("\nchecking novelty against the provided HER2 tables ...")
    known = known_cdrh3()
    print(f"  {len(known):,} distinct CDRH3s loaded")
    u["in_tables"] = [s in known for s in u.cdrh3]

    # HARD GATES, applied before any score is reported
    u["legal"] = (u.identity < 70.0) & ~u.in_tables & ~u.has_cys
    print("\n=== HARD GATE ===")
    print(f"  identity < 70%      : {int((u.identity < 70).sum())}/{len(u)}")
    print(f"  novel vs tables     : {int((~u.in_tables).sum())}/{len(u)}")
    print(f"  no free cysteine    : {int((~u.has_cys).sum())}/{len(u)}")
    print(f"  PASS ALL            : {int(u.legal.sum())}/{len(u)}")

    spr, facs = load()
    pf = PreFilter().fit(spr, facs)
    sc = pf.score(u.cdrh3.tolist())
    u = u.merge(sc[["cdrh3", "p_binder", "pred_log10_kd"]], on="cdrh3")
    u["pred_kd_nM"] = 10 ** u.pred_log10_kd

    pd.set_option("display.width", 220)
    cols = ["cdrh3", "n_mut", "identity", "in_tables", "has_asn", "n_designs",
            "p_binder", "pred_kd_nM"]
    ok = u[u.legal].sort_values("pred_kd_nM")
    print(f"\n=== {len(ok)} LEGAL, NOVEL RFD4 DESIGNS (ranked; see caveat) ===")
    print(ok[cols].round(3).to_string(index=False) if len(ok) else "  none")
    rej = u[~u.legal]
    if len(rej):
        print(f"\n=== {len(rej)} REJECTED ===")
        print(rej[cols + ["has_cys"]].round(3).to_string(index=False))

    u.to_csv("rfd4_designs_gated.csv", index=False)
    print("\nwrote rfd4_designs_gated.csv")
    print("\nCAVEAT: pred_kd_nM is from a ranker whose floor is 1.135 nM and whose")
    print("Spearman below 1 nM is +0.17 (negative in the >=4 band). Use it to drop")
    print("disasters, NOT to claim anything beats trastuzumab's measured 0.152 nM.")


if __name__ == "__main__":
    main()
