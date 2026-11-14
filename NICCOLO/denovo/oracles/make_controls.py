#!/usr/bin/env python3
"""Build the CONTROL PANEL csv for the HER2 refolding oracle.

All sequences are taken from real structures (provenance in PROV below), never from memory.
Run:  python3 -I make_controls.py  ->  controls.csv
"""
import csv, os, random

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- antigen
# HER2 domain IV. 1N8Z is deposited with an unresolved loop 581-590, which the
# rfd4 input file splits into chain C (489-580) and chain D (591-607).
# The TRUE 581-590 sequence is recoverable from PDB 1S78 entity 1 (full HER2
# ectodomain, residues ...CVARCPS|GVKPDLSYMP|IWKFPDEEG...), so we do NOT have to
# choose between "close the gap artificially" and "two chains": we fill the gap
# with the real sequence and submit ONE contiguous chain 489-607.
HER2_C = "CHQLCARGHCWGPGPTQCVNCSQFLRGQECVEECRVLQGLPREYVNARHCLPCHPECQPQNGSVTCFGPEADQCVACAHYKDPPFCVARCPS"  # 489-580
HER2_GAP = "GVKPDLSYMP"                                                                                    # 581-590, from 1S78
HER2_D = "IWKFPDEEGACQPCPIN"                                                                               # 591-607
HER2_DIV = HER2_C + HER2_GAP + HER2_D                                                                      # 489-607, 119 aa

# ---------------------------------------------------------------- antibodies
# trastuzumab Fv, from apheris/rfd4/Data/12-her2/her2_fv_1n8z.cif (1N8Z)
TRAS_H = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
TRAS_L = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"

# pertuzumab Fv, from PDB 1S78 entities 3 (heavy) / 2 (light), Fab trimmed to Fv
PERT_H = "EVQLVESGGGLVQPGGSLRLSCAASGFTFTDYTMDWVRQAPGKGLEWVADVNPNSGGSIYNQRFKGRFTLSVDRSKNTLYLQMNSLRAEDTAVYYCARNLGPSFYFDYWGQGTLVTVSS"
PERT_L = "DIQMTQSPSSLSASVGDRVTITCKASQDVSIGVAWYQQKPGKAPKLLIYSASYRYTGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYYIYPYTFGQGTKVEIK"

# anti-HER2 VHHs (His/FLAG tags trimmed at the WGQGT.VTVSS framework-4 end)
VHH_27A05 = "EVQLVESGGGLVQAGGSLRLSCVASGIPFSTRTMAWYRQPPGNERDWVATIRSGAPVYADSVKGRFTVSRDNAKNTLYLQMNSLEPEDTATYYCWDVNGDIWGQGTLVTVSS"   # 9MTE_1
VHH_2RS15D = "QVQLQESGGGSVQAGGSLKLTCAASGYIFNSCGMGWYRQSPGRERELVSRISGDGDTWHKESVKGRFTISQDNVKKTLYLQMNSLKPEDTAVYFCAVCYNLETYWGQGTQVTVSS"  # 5MY6_2

# irrelevant VHH negative control: cAb-Lys3, anti-hen-egg-lysozyme, PDB 1MEL chain B
# (HA tag GRYPYDVPDYGSGRA trimmed). Used INSTEAD of caplacizumab: the repo's
# apheris/rfd4/Data/03-vhh/pmvk_caplacizumab.cif.zst has its CDRs masked to single
# 'X' residues (it is an RFD4 design template, not a usable sequence) and its antigen
# is phosphomevalonate kinase, not von Willebrand factor.
VHH_CABLYS3 = "DVQLQASGGGSVQAGGSLRLSCAASGYTIGPYCMGWFRQAPGKEREGVAAINMGGGITYYADSVKGRFTISQDNAKNTVYLLMNSLEPEDTAIYYCAADSTIYASYYECGHGLSTGGYGYDSWGQGTQVTVSS"

PROV = {
    "tras_wt_fv": "1N8Z via apheris/rfd4/Data/12-her2/her2_fv_1n8z.cif chain B(VH)+A(VL)",
    "tras_scram_cdr": "tras_wt_fv framework, 6 CDRs shuffled in place, random.Random(0)",
    "pert_wt_fv": "PDB 1S78 entity 3(VH)+2(VL), Fab trimmed to Fv",
    "vhh_27A05": "PDB 9MTE entity 1, tag trimmed",
    "vhh_2Rs15d": "PDB 5MY6 entity 2",
    "vhh_cablys3": "PDB 1MEL chain B, HA tag trimmed",
}

# ------------------------------------------------- scrambled-CDR decoy
# Kabat CDRs located by exact substring so the framework is bit-identical.
CDR_H = ["DTYIH", "RIYPTNGYTRYADSVK", "SRWGGDGFYAMDY"]
CDR_L = ["RASQDVNTAVA", "SASFLYS", "QQHYTTPPT"]


def scramble(seq, cdrs, rng):
    """Shuffle each CDR in place; framework and composition unchanged."""
    out = seq
    for c in cdrs:
        i = out.index(c)                      # raises if the CDR moved -> fail loud
        chars = list(c)
        while True:
            rng.shuffle(chars)
            if "".join(chars) != c:           # guarantee it actually changed
                break
        out = out[:i] + "".join(chars) + out[i + len(c):]
    return out


def main():
    rng = random.Random(0)
    scr_h = scramble(TRAS_H, CDR_H, rng)
    scr_l = scramble(TRAS_L, CDR_L, rng)
    assert len(scr_h) == len(TRAS_H) and len(scr_l) == len(TRAS_L)
    assert sorted(scr_h) == sorted(TRAS_H), "scramble must preserve composition"
    assert scr_h != TRAS_H and scr_l != TRAS_L

    rows = [
        # sequence_id, kind, heavy_sequence, light_sequence, antigen_mode, expect, note
        dict(sequence_id="tras_wt_fv", kind="control_pos", heavy_sequence=TRAS_H,
             light_sequence=TRAS_L, antigen_mode="single",
             expect="binder", note="POSITIVE. Epitope = HER2 domain IV 557-605, fully inside our 489-607 construct (measured from 1N8Z)."),
        dict(sequence_id="tras_scram_cdr", kind="control_decoy", heavy_sequence=scr_h,
             light_sequence=scr_l, antigen_mode="single",
             expect="non-binder", note="DECOY FLOOR. Same framework + same CDR composition, CDR order shuffled."),
        dict(sequence_id="pert_wt_fv", kind="control_neg_epitope_absent", heavy_sequence=PERT_H,
             light_sequence=PERT_L, antigen_mode="single",
             expect="non-binder", note="Pertuzumab epitope = HER2 domain II 245-315 (measured from 1S78). ABSENT from a domain IV construct."),
        dict(sequence_id="vhh_27A05", kind="control_neg_epitope_mismatch", heavy_sequence=VHH_27A05,
             light_sequence="", antigen_mode="single",
             expect="non-binder", note="anti-HER2 VHH but epitope is domain I+III (measured from 9MTE); only 497-500 touches domain IV."),
        dict(sequence_id="vhh_2Rs15d", kind="control_neg_epitope_mismatch", heavy_sequence=VHH_2RS15D,
             light_sequence="", antigen_mode="single",
             expect="non-binder", note="anti-HER2 VHH, epitope domain I 109-217 (measured from 5MY6)."),
        dict(sequence_id="vhh_cablys3", kind="control_neg_irrelevant", heavy_sequence=VHH_CABLYS3,
             light_sequence="", antigen_mode="single",
             expect="non-binder", note="IRRELEVANT VHH: anti-lysozyme cAb-Lys3. True unrelated-antigen negative."),
        dict(sequence_id="tras_wt_fv_2chainAg", kind="control_repr_check", heavy_sequence=TRAS_H,
             light_sequence=TRAS_L, antigen_mode="split",
             expect="binder", note="Same as tras_wt_fv but antigen submitted as TWO chains (489-580 / 591-607) to test the representation choice and to anchor pair-matrix indexing."),
    ]
    out = os.path.join(HERE, "controls.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("wrote", out, len(rows), "controls")
    print("HER2 domain IV construct 489-607:", len(HER2_DIV), "aa")
    for r in rows:
        n = len(r["heavy_sequence"]) + len(r["light_sequence"])
        print(f"  {r['sequence_id']:<22} {r['kind']:<32} ab={n:>3}aa expect={r['expect']}")
    print("\nscrambled heavy:", scr_h)
    print("scrambled light:", scr_l)


if __name__ == "__main__":
    main()
