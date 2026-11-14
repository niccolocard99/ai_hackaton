"""Score LEO's drafted CDRH3 candidates with NICCOLO's validated pre-filter.

WINDOW. The pre-filter's own constant is `TZ = "WGGDGFYAMD"`, L = 10 -- positions
99-108. Y109 is OUTSIDE the scored window, which collapses the 8 drafted 11-mers
to 6 distinct 10-mers: A1/A2 share one, C1/C2 share another, and F1/B1 drop to
3 substitutions = exactly 70.0% identity = FAIL under the "strictly under 70%"
rule. So the 8 drafts are 4 legal designs.

CONTROLS are included on purpose. A predicted KD means nothing without known
answers on the same scale, so trastuzumab itself and the two best measured legal
variants from EXP-008 are scored alongside. If the model cannot rank the
controls correctly, its verdict on the candidates is not usable.

NOVELTY. `score()` implements stages 1-2 only; stage 0 (novelty) is a hard rule
and is the single likeliest cause of outright rejection, so it is checked here
directly against the full HER2 tables rather than assumed.
"""
import io
import sys
import zipfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from prefilter import TZ, PreFilter, load, nmut  # noqa: E402

# 10-mer window (99-108) -> the drafted ids that collapse onto it
CANDIDATES = {
    "WGGEGFYSLE": "A1/A2",
    "WGGEGFYALE": "F1",
    "WGAEGFYSLD": "C1/C2",
    "WGAEGFYSME": "D1",
    "WGGEGFYSIE": "E1",
    "WGGEGFYSLD": "B1",
}
# measured KD in nM, from EXP-008's calibration table
CONTROLS = {
    "WGGDGFYAMD": ("trastuzumab (the bar)", 0.152),
    "WGGSGPYLFD": ("EXP-008 best legal", 0.095),
    "WGAEGFYYWD": ("EXP-008 2nd legal", 0.368),
}

FACS_ZIP = Path(
    "/home/hackathon/hackathon-team-3/apheris/hackathon-data/her2-trastuzumab"
    "/tz_her2_facs_3classes_524346.csv.zip"
)


def subs(seq):
    return [f"{a}{99+i}{b}" for i, (a, b) in enumerate(zip(TZ, seq)) if a != b]


def known_cdrh3():
    """Every CDRH3 in the provided HER2 tables, for the stage-0 novelty gate."""
    seen = set()
    spr = pd.read_csv("data/trastuzumab_her2_spr_all_768.csv")
    seen |= set(spr.cdrh3.dropna())
    with zipfile.ZipFile(FACS_ZIP) as z:
        name = next(n for n in z.namelist() if n.endswith(".csv"))
        with z.open(name) as fh:
            for chunk in pd.read_csv(io.TextIOWrapper(fh, "utf-8"), chunksize=200_000):
                col = "seq" if "seq" in chunk.columns else chunk.columns[0]
                seen |= set(chunk[col].dropna())
    return seen


def main():
    spr, facs = load()
    print(f"fitting on {len(spr)} SPR variants and {len(facs)} FACS rows ...")
    pf = PreFilter().fit(spr, facs)

    print("checking novelty against the provided HER2 tables ...")
    known = known_cdrh3()
    print(f"  {len(known):,} distinct CDRH3s loaded")

    seqs = list(CANDIDATES) + [c for c in CONTROLS if c not in CANDIDATES]
    df = pf.score(seqs)
    df["pred_kd_nM"] = 10 ** df.pred_log10_kd
    df["identity"] = (10 - df.n_mut) / 10 * 100
    df["legal_lt70"] = df.identity < 70.0
    df["in_tables"] = [s in known for s in df.cdrh3]
    df["who"] = [CANDIDATES.get(s, CONTROLS.get(s, ("", 0))[0]) for s in df.cdrh3]
    df["measured_kd_nM"] = [CONTROLS.get(s, ("", float("nan")))[1] for s in df.cdrh3]
    df["substitutions"] = [" ".join(subs(s)) or "-" for s in df.cdrh3]

    cols = [
        "cdrh3", "who", "n_mut", "identity", "legal_lt70", "in_tables",
        "p_binder", "pred_kd_nM", "measured_kd_nM", "substitutions",
    ]
    out = df[cols].sort_values("pred_kd_nM").reset_index(drop=True)
    pd.set_option("display.width", 200, "display.max_colwidth", 44)
    print("\n=== CONTROLS FIRST: can the model rank what we already measured? ===")
    print(out[out.measured_kd_nM.notna()].round(3).to_string(index=False))
    print("\n=== ALL, ranked by predicted KD (lower is better) ===")
    print(out.round(3).to_string(index=False))

    tz_pred = float(df.loc[df.cdrh3 == TZ, "pred_kd_nM"].iloc[0])
    print(f"\ntrastuzumab predicted {tz_pred:.3f} nM vs measured 0.152 nM")
    legal = out[out.legal_lt70 & ~out.in_tables & out.measured_kd_nM.isna()]
    print(f"\nlegal, novel candidates: {len(legal)} of {len(CANDIDATES)} distinct 10-mers")
    if len(legal):
        best = legal.iloc[0]
        print(f"best candidate {best.cdrh3} ({best.who}) predicted {best.pred_kd_nM:.3f} nM")
        print(f"  => {'BEATS' if best.pred_kd_nM < tz_pred else 'LOSES TO'} trastuzumab's predicted {tz_pred:.3f} nM")
    out.to_csv("leo_candidates_scored.csv", index=False)
    print("\nwrote leo_candidates_scored.csv")


if __name__ == "__main__":
    main()
