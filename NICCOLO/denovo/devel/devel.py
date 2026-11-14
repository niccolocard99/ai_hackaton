#!/usr/bin/env python3
"""
devel.py -- developability + thermostability proxy scoring for de novo VHH / Fv designs.

CAMPAIGN OUTPUT RULE
====================
Every column produced by a model, heuristic or calculation carries a `pred_`
prefix. NOTHING in this table is a measured value. If a lab number ever lands in
the same file it must use a different prefix (e.g. `meas_`). Predicted and
measured values must never share a column.

WHAT PRODUCED EACH COLUMN
=========================
All `pred_*` columns below are computed by THIS script from sequence (and, when a
structure is supplied, from that structure's coordinates). Provenance per block:

  pred_net_charge_ph74, pred_pI
      Henderson-Hasselbalch over side-chain + terminal pKa (EMBOSS pKa set).
      Deterministic calculation, not a learned model.

  pred_gravy, pred_mean_hydrophobicity_kd
      Kyte-Doolittle hydropathy, mean over the chain(s). Deterministic.

  pred_hydrophobic_patch_* (STRUCTURE)
      Shrake-Rupley SASA (biotite) -> exposed apolar residues (rel. SASA > 0.25)
      -> single-linkage clustering on CB-CB <= 8 A -> largest cluster.
      This is a PROXY FOR Aggrescan3D, not Aggrescan3D. A3D scores an
      aggregation free-energy per residue from a structure-corrected scale; this
      only measures apolar surface contiguity. Direction of effect agrees
      (bigger exposed apolar patch = worse), the units do not.

  pred_hydrophobic_patch_seq_proxy (SEQUENCE, LABELLED PROXY)
      Max mean Kyte-Doolittle over a 5-residue sliding window. Used when no
      structure is available. Strictly weaker than the structural proxy: it
      cannot know whether the window is buried. Treat as a screen, not a score.

  pred_cdr_*
      CDRs located by regex anchors on the conserved framework motifs
      (FR1 Cys, W[VFYGL]RQ, R/K-x-TIS, YxC ... WGxG for H; ITC, W[YLF]Q,
      FGxG for L). These are ESTIMATES, NOT ANARCI/IMGT numbering. CDR1/CDR2
      are Kabat-like; CDR3 is IMGT-like (conserved Cys+1 .. W of WGxG-1), so
      trastuzumab CDR-H3 comes out as the 13-mer SRWGGDGFYAMDY rather than the
      Kabat 11-mer WGGDGFYAMDY. pred_cdr_parse_ok flags rows where the anchors
      did not all hit -- do not trust CDR columns when it is 0.

  pred_liab_* and pred_liab_detail
      Regex motif scan with 1-based positions, per chain.

  pred_vhh_hallmark_n, pred_vhh_hallmark_detail, pred_disulfide_*
      FR2 tetrad read positionally from the conserved W (indexed as Kabat 36)
      at Kabat 37/44/45/47; conserved Cys pair located by framework anchors and,
      when a structure is supplied, confirmed by SG-SG distance. That pair is
      Kabat C22-C96 == IMGT C23-C104 (the brief's "C22-C96 (IMGT)" mislabels
      Kabat numbering as IMGT; the positions 22/96 are the Kabat ones and are
      what this script reports).

  pred_tm_c, pred_tm_source
      See TM STATUS below.

TM STATUS (as verified on this VM, 14 Nov 2026)
===============================================
There is NO installed thermostability predictor for a VHH on this machine.
  * TEMPRO (Alvarez & Dean, Sci Rep 2024, nanobody Tm, reported MAE 4.03 C) is
    the only published VHH Tm model. Its repo ships Keras models but requires
    ESM-2 embeddings we could not generate inside the time box -- see DEVEL.md.
  * Nebius `stab-ddg` predicts BINDING ddG for point mutants (an affinity
    proxy), NOT folding stability -- it is not a Tm predictor.
  * JapanFold `opendde`/`opendde-abag` are docking models, not Tm predictors.
Therefore `pred_tm_c` is EMPTY by default and `pred_tm_source` says
"no_tm_predictor_available". If TEMPRO is later made to run, feed its output in
via --tempro-csv (sequence_id,tm_c) and pred_tm_source becomes "TEMPRO_keras".
A crude labelled surrogate is always reported separately as
`pred_tm_surrogate_rank` -- a RANK, not a temperature, and not a Tm.

USAGE
=====
  devel.py --csv designs.csv --out devel_scores.csv
           [--structure-dir DIR]   # <sequence_id>.pdb / .cif looked up here
           [--structure PATH]      # single structure, applied to a single-row CSV
           [--tempro-csv PATH]     # sequence_id,tm_c from a working TEMPRO
Input CSV columns: sequence_id, kind, heavy_sequence, light_sequence
(`light_sequence` may be empty for a VHH; `kind` is free text, e.g. vhh / scfv.)
"""

import argparse
import csv
import math
import os
import re
import sys

# ----------------------------------------------------------------------------
# scales
# ----------------------------------------------------------------------------

# Kyte & Doolittle 1982 hydropathy
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# EMBOSS pKa set (as used by Biopython's IsoelectricPoint)
PKA_SIDE = {"C": 8.5, "D": 3.9, "E": 4.1, "H": 6.5, "K": 10.0, "R": 12.5, "Y": 10.1}
PKA_NTERM = 8.6
PKA_CTERM = 3.6
POSITIVE = {"K", "R", "H"}
NEGATIVE = {"D", "E", "C", "Y"}

APOLAR = set("AVILMFWYC")          # for patch / CDR hydrophobicity
AA20 = set(KD)

# Maximum residue SASA (Tien et al. 2013, theoretical, A^2) for relative SASA
MAX_SASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0, "Q": 225.0,
    "E": 223.0, "G": 104.0, "H": 224.0, "I": 197.0, "L": 201.0, "K": 236.0,
    "M": 224.0, "F": 240.0, "P": 159.0, "S": 155.0, "T": 172.0, "W": 285.0,
    "Y": 263.0, "V": 174.0,
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}

# ----------------------------------------------------------------------------
# basic physchem
# ----------------------------------------------------------------------------


def clean(seq):
    return "".join(c for c in (seq or "").strip().upper() if c.isalpha())


def charge_at_ph(seq, ph):
    """Net charge by Henderson-Hasselbalch. Termini included once per chain set."""
    q = 0.0
    # N-terminus (positive when protonated)
    q += 1.0 / (1.0 + 10 ** (ph - PKA_NTERM))
    # C-terminus (negative when deprotonated)
    q -= 1.0 / (1.0 + 10 ** (PKA_CTERM - ph))
    for aa in seq:
        pka = PKA_SIDE.get(aa)
        if pka is None:
            continue
        if aa in POSITIVE:
            q += 1.0 / (1.0 + 10 ** (ph - pka))
        else:
            q -= 1.0 / (1.0 + 10 ** (pka - ph))
    return q


def isoelectric_point(seq):
    lo, hi = 2.0, 13.0
    for _ in range(100):
        mid = (lo + hi) / 2.0
        if charge_at_ph(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def gravy(seq):
    vals = [KD[a] for a in seq if a in KD]
    return sum(vals) / len(vals) if vals else float("nan")


def max_hydrophobic_window(seq, w=5):
    """Sequence-only proxy for an exposed apolar patch. LABELLED PROXY."""
    vals = [KD.get(a, 0.0) for a in seq]
    if len(vals) < w:
        return float("nan"), -1
    best, best_i = -99.0, -1
    for i in range(len(vals) - w + 1):
        m = sum(vals[i:i + w]) / w
        if m > best:
            best, best_i = m, i + 1  # 1-based start
    return best, best_i


# ----------------------------------------------------------------------------
# CDR location (Kabat-like regex anchors -- NOT ANARCI)
# ----------------------------------------------------------------------------

RE_H_FR2 = re.compile(r"W[VFYGLIA]RQ")
RE_H_FR3 = re.compile(r"[RK][FLVSAM]T[IVLM][STF]")
RE_H_CYS2 = re.compile(r"[YFHVAC][YFHVCLM]C")
RE_H_FR4 = re.compile(r"WG[QKRPAE]G")
RE_H_CYS1 = re.compile(r"[LIVMF][STRQKE]C")

RE_L_CYS1 = re.compile(r"[LIVMF][TSRQ]C")
RE_L_FR2 = re.compile(r"W[YLFH][QLH]")
RE_L_CYS2 = re.compile(r"[YFH][YFHCV]C")
RE_L_FR4 = re.compile(r"FG[QGPSTA]G")


def _last(regex, seq, start=0):
    hits = [m for m in regex.finditer(seq, start)]
    return hits[-1] if hits else None


def find_cdrs_heavy(seq):
    """Return (dict of cdr -> (1-based start, end inclusive, subseq), ok_flag)."""
    out, ok = {}, True
    m_c1 = RE_H_CYS1.search(seq[:40])
    m_fr2 = RE_H_FR2.search(seq)
    if m_c1 and m_fr2 and m_fr2.start() > m_c1.end():
        s, e = m_c1.end() + 3, m_fr2.start()          # Kabat H1 ~ 31..35
        out["H1"] = (s + 1, e, seq[s:e])
    else:
        ok = False
    if m_fr2:
        s = m_fr2.start() + 14                         # end of FR2 (W36..A49)
        m_fr3 = RE_H_FR3.search(seq, s)
        if m_fr3:
            out["H2"] = (s + 1, m_fr3.start(), seq[s:m_fr3.start()])
        else:
            ok = False
    else:
        ok = False
    m_c2 = _last(RE_H_CYS2, seq[: int(len(seq) * 0.95)])
    m_fr4 = _last(RE_H_FR4, seq)
    if m_c2 and m_fr4 and m_fr4.start() > m_c2.end():
        out["H3"] = (m_c2.end() + 1, m_fr4.start(), seq[m_c2.end():m_fr4.start()])
    else:
        ok = False
    return out, ok


def find_cdrs_light(seq):
    out, ok = {}, True
    m_c1 = RE_L_CYS1.search(seq[:40])
    m_fr2 = RE_L_FR2.search(seq)
    if m_c1 and m_fr2 and m_fr2.start() > m_c1.end():
        out["L1"] = (m_c1.end() + 1, m_fr2.start(), seq[m_c1.end():m_fr2.start()])
    else:
        ok = False
    if m_fr2:
        s = m_fr2.start() + 17                         # FR2 L is ~15-17 long
        out["L2"] = (s + 1, s + 7, seq[s:s + 7])
    else:
        ok = False
    m_c2 = _last(RE_L_CYS2, seq[: int(len(seq) * 0.95)])
    m_fr4 = _last(RE_L_FR4, seq)
    if m_c2 and m_fr4 and m_fr4.start() > m_c2.end():
        out["L3"] = (m_c2.end() + 1, m_fr4.start(), seq[m_c2.end():m_fr4.start()])
    else:
        ok = False
    return out, ok


# ----------------------------------------------------------------------------
# liability motifs
# ----------------------------------------------------------------------------

LIAB_PATTERNS = [
    ("nglyc",      re.compile(r"N[^P][ST]"),  "N-glycosylation sequon NxS/NxT"),
    ("deamid_NG",  re.compile(r"N(?=G)"),     "deamidation NG"),
    ("deamid_NS",  re.compile(r"N(?=S)"),     "deamidation NS"),
    ("isomer_DG",  re.compile(r"D(?=G)"),     "isomerisation DG"),
    ("hydrol_DP",  re.compile(r"D(?=P)"),     "acid-labile DP"),
]


def scan_liabilities(seq, chain_tag, cdr_spans):
    """Return counts dict and a detail string with 1-based positions."""
    counts, details = {}, []
    for name, rx, _desc in LIAB_PATTERNS:
        pos = [m.start() + 1 for m in rx.finditer(seq)]
        counts[name] = len(pos)
        for p in pos:
            details.append("%s:%s@%d%s" % (chain_tag, name, p,
                                           "(CDR)" if in_cdr(p, cdr_spans) else ""))
    # cysteines
    cys = [i + 1 for i, a in enumerate(seq) if a == "C"]
    counts["cys_total"] = len(cys)
    counts["cys_odd"] = 1 if len(cys) % 2 else 0
    # cysteines beyond the two conserved framework ones are extra/free risk
    counts["cys_extra"] = max(0, len(cys) - 2)
    for p in cys:
        details.append("%s:CYS@%d%s" % (chain_tag, p,
                                        "(CDR)" if in_cdr(p, cdr_spans) else ""))
    # Trp oxidation -- CDR/surface Trp is the risk, framework Trp is buried
    trp = [i + 1 for i, a in enumerate(seq) if a == "W"]
    trp_cdr = [p for p in trp if in_cdr(p, cdr_spans)]
    counts["trp_total"] = len(trp)
    counts["trp_cdr"] = len(trp_cdr)
    for p in trp_cdr:
        details.append("%s:TRPox@%d(CDR)" % (chain_tag, p))
    met = [i + 1 for i, a in enumerate(seq) if a == "M"]
    met_cdr = [p for p in met if in_cdr(p, cdr_spans)]
    counts["met_total"] = len(met)
    counts["met_cdr"] = len(met_cdr)
    for p in met:
        details.append("%s:METox@%d%s" % (chain_tag, p,
                                          "(CDR)" if in_cdr(p, cdr_spans) else ""))
    return counts, details


def in_cdr(pos1, cdr_spans):
    for (s, e) in cdr_spans:
        if s <= pos1 <= e:
            return True
    return False


# ----------------------------------------------------------------------------
# VHH hallmarks + conserved disulfide
# ----------------------------------------------------------------------------

# VHH-like identity at Kabat 37/44/45/47 (camelid hallmark tetrad);
# classical VH is V37 / G44 / L45 / W47.
HALLMARK_VHH = {37: set("FYLH"), 44: set("EQD"), 45: set("RCLP"), 47: set("FGLS")}
HALLMARK_VH = {37: set("V"), 44: set("G"), 45: set("L"), 47: set("W")}


def vhh_hallmarks(seq):
    """Read Kabat 37/44/45/47 positionally from the conserved W (= Kabat 36)."""
    m = RE_H_FR2.search(seq)
    if not m:
        return None, "FR2_anchor_not_found", ""
    w = m.start()                      # index of W, Kabat 36
    res, n = {}, 0
    for kab in (37, 44, 45, 47):
        idx = w + (kab - 36)
        res[kab] = seq[idx] if idx < len(seq) else "?"
        if res[kab] in HALLMARK_VHH[kab]:
            n += 1
    detail = ";".join("K%d=%s%s" % (k, res[k],
                                    "*" if res[k] in HALLMARK_VHH[k] else
                                    ("=VH" if res[k] in HALLMARK_VH[k] else ""))
                      for k in (37, 44, 45, 47))
    return n, detail, seq[w:w + 14]


def conserved_cys(seq):
    """Locate the two conserved framework Cys (IMGT 23 / 104; Kabat 22 / 92-ish).

    Returns (pos1_fr1, pos1_fr3) 1-based, or None where not found.
    """
    c1 = RE_H_CYS1.search(seq[:40])
    p1 = c1.end() if c1 else None       # 1-based position of the Cys
    c2 = _last(RE_H_CYS2, seq[: int(len(seq) * 0.95)])
    p2 = c2.end() if c2 else None
    return p1, p2


# ----------------------------------------------------------------------------
# structure-based patch (Aggrescan3D PROXY)
# ----------------------------------------------------------------------------


def structure_patch(path, chains=None):
    """Largest exposed-apolar-surface cluster. PROXY for Aggrescan3D.

    Returns dict or {'error': ...}. Requires numpy + biotite.
    """
    try:
        import numpy as np
        import biotite.structure as struc
        from biotite.structure.io.pdb import PDBFile
        from biotite.structure.io.pdbx import CIFFile, get_structure
    except Exception as exc:                                    # noqa: BLE001
        return {"error": "biotite/numpy unavailable: %s" % exc}
    try:
        if path.lower().endswith((".cif", ".mmcif")):
            arr = get_structure(CIFFile.read(path), model=1)
        else:
            arr = PDBFile.read(path).get_structure(model=1)
        arr = arr[struc.filter_amino_acids(arr)]
        arr = arr[~np.isnan(arr.coord).any(axis=-1)]
        if chains:
            arr = arr[np.isin(arr.chain_id, list(chains))]
        if arr.array_length() == 0:
            return {"error": "no atoms after filtering"}
        sasa_atom = struc.sasa(arr, vdw_radii="Single", point_number=200)
        sasa_atom = np.nan_to_num(sasa_atom)
        res_sasa = struc.apply_residue_wise(arr, sasa_atom, np.sum)
        res_ids, res_names = struc.get_residues(arr)
        codes = [THREE_TO_ONE.get(n, "X") for n in res_names]
        # representative coordinate per residue: CB, else CA
        reps = []
        for rid, ch in zip(res_ids, [arr.chain_id[arr.res_id == r][0] for r in res_ids]):
            sel = arr[(arr.res_id == rid) & (arr.chain_id == ch)]
            cb = sel[sel.atom_name == "CB"]
            ca = sel[sel.atom_name == "CA"]
            reps.append((cb.coord[0] if cb.array_length() else
                         (ca.coord[0] if ca.array_length() else sel.coord[0])))
        reps = np.asarray(reps)
        exposed = []
        for i, (c, s) in enumerate(zip(codes, res_sasa)):
            if c in APOLAR and MAX_SASA.get(c) and s / MAX_SASA[c] > 0.25:
                exposed.append(i)
        if not exposed:
            return {"n_exposed_apolar": 0, "patch_n": 0, "patch_sasa": 0.0,
                    "total_apolar_sasa": 0.0}
        # single-linkage cluster on CB-CB <= 8 A
        parent = {i: i for i in exposed}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for a_i in range(len(exposed)):
            for b_i in range(a_i + 1, len(exposed)):
                ia, ib = exposed[a_i], exposed[b_i]
                if np.linalg.norm(reps[ia] - reps[ib]) <= 8.0:
                    ra, rb = find(ia), find(ib)
                    if ra != rb:
                        parent[ra] = rb
        groups = {}
        for i in exposed:
            groups.setdefault(find(i), []).append(i)
        best = max(groups.values(), key=lambda g: sum(res_sasa[i] for i in g))
        return {
            "n_exposed_apolar": len(exposed),
            "patch_n": len(best),
            "patch_sasa": float(sum(res_sasa[i] for i in best)),
            "total_apolar_sasa": float(sum(res_sasa[i] for i in exposed)),
            "patch_residues": ";".join("%s%d" % (codes[i], res_ids[i]) for i in sorted(best)),
        }
    except Exception as exc:                                    # noqa: BLE001
        return {"error": "%s: %s" % (type(exc).__name__, exc)}


def structure_disulfides(path, chains=None):
    """Count SG-SG pairs within 2.5 A; report whether the conserved pair is closed."""
    try:
        import numpy as np
        import biotite.structure as struc
        from biotite.structure.io.pdb import PDBFile
        from biotite.structure.io.pdbx import CIFFile, get_structure
    except Exception as exc:                                    # noqa: BLE001
        return {"error": str(exc)}
    try:
        if path.lower().endswith((".cif", ".mmcif")):
            arr = get_structure(CIFFile.read(path), model=1)
        else:
            arr = PDBFile.read(path).get_structure(model=1)
        arr = arr[struc.filter_amino_acids(arr)]
        arr = arr[~np.isnan(arr.coord).any(axis=-1)]
        if chains:
            arr = arr[np.isin(arr.chain_id, list(chains))]
        sg = arr[(arr.atom_name == "SG") & (arr.res_name == "CYS")]
        pairs = []
        for i in range(sg.array_length()):
            for j in range(i + 1, sg.array_length()):
                d = float(np.linalg.norm(sg.coord[i] - sg.coord[j]))
                if d <= 2.5:
                    pairs.append("%s%d-%s%d@%.2fA" % (sg.chain_id[i], sg.res_id[i],
                                                      sg.chain_id[j], sg.res_id[j], d))
        return {"n_ss": len(pairs), "ss_detail": ";".join(pairs),
                "n_free_cys": sg.array_length() - 2 * len(pairs)}
    except Exception as exc:                                    # noqa: BLE001
        return {"error": "%s: %s" % (type(exc).__name__, exc)}


# ----------------------------------------------------------------------------
# Tm surrogate (a RANK, never a temperature)
# ----------------------------------------------------------------------------


def tm_surrogate_rank(heavy, hallmark_n, ss_ok, cdr3_len, patch_metric, gravy_v):
    """Ordinal 0-100 'likely to be stable' rank. NOT a temperature. NOT a Tm.

    Built from the directions reported in the VHH literature: an intact
    C23-C104 disulfide and the hydrophilic FR2 tetrad raise Tm; long CDR3s and
    exposed apolar surface lower it. Weights are judgement, not a fit to data --
    this has never been regressed against a Tm dataset, so use it only to rank
    designs against each other inside this campaign.
    """
    score = 50.0
    score += 15.0 if ss_ok else -25.0
    score += 3.0 * (hallmark_n or 0)
    if cdr3_len:
        score -= max(0.0, (cdr3_len - 14)) * 1.2
    if patch_metric is not None and not math.isnan(patch_metric):
        score -= max(0.0, (patch_metric - 1.0)) * 6.0
    if not math.isnan(gravy_v):
        score -= max(0.0, (gravy_v + 0.35)) * 20.0
    return round(max(0.0, min(100.0, score)), 1)


# ----------------------------------------------------------------------------
# main per-design scoring
# ----------------------------------------------------------------------------


def score_design(row, struct_path=None, tempro_tm=None):
    sid = (row.get("sequence_id") or "").strip()
    kind = (row.get("kind") or "").strip()
    H = clean(row.get("heavy_sequence"))
    L = clean(row.get("light_sequence"))
    full = H + L
    rec = {"sequence_id": sid, "kind": kind,
           "len_heavy": len(H), "len_light": len(L)}

    bad = sorted(set(full) - AA20)
    rec["pred_nonstandard_residues"] = "".join(bad)

    rec["pred_net_charge_ph74"] = round(charge_at_ph(full, 7.4), 2)
    rec["pred_pI"] = round(isoelectric_point(full), 2)
    rec["pred_gravy"] = round(gravy(full), 3)
    rec["pred_mean_hydrophobicity_kd"] = rec["pred_gravy"]

    cdrH, okH = find_cdrs_heavy(H)
    cdrL, okL = (find_cdrs_light(L) if L else ({}, True))
    rec["pred_cdr_parse_ok"] = int(bool(okH and okL))
    spansH = [(s, e) for (s, e, _q) in cdrH.values()]
    spansL = [(s, e) for (s, e, _q) in cdrL.values()]

    for k in ("H1", "H2", "H3"):
        rec["pred_cdr%s" % k] = cdrH.get(k, (0, 0, ""))[2]
    for k in ("L1", "L2", "L3"):
        rec["pred_cdr%s" % k] = cdrL.get(k, (0, 0, ""))[2]
    rec["pred_cdrH3_len"] = len(rec["pred_cdrH3"])

    cdr_res = "".join(rec["pred_cdr%s" % k] for k in ("H1", "H2", "H3", "L1", "L2", "L3"))
    n_cdr = len(cdr_res)
    rec["pred_cdr_total_len"] = n_cdr
    rec["pred_cdr_frac_arg"] = round(cdr_res.count("R") / n_cdr, 3) if n_cdr else ""
    rec["pred_cdr_frac_lys"] = round(cdr_res.count("K") / n_cdr, 3) if n_cdr else ""
    rec["pred_cdr_frac_hydrophobic"] = round(
        sum(1 for a in cdr_res if a in APOLAR) / n_cdr, 3) if n_cdr else ""
    rec["pred_cdr_net_charge_ph74"] = round(
        sum(1.0 for a in cdr_res if a in "KR") + 0.1 * cdr_res.count("H")
        - sum(1.0 for a in cdr_res if a in "DE"), 2) if n_cdr else ""

    # liabilities
    allc, det = {}, []
    for tag, s, spans in (("H", H, spansH), ("L", L, spansL)):
        if not s:
            continue
        c, d = scan_liabilities(s, tag, spans)
        for k, v in c.items():
            allc[k] = allc.get(k, 0) + v
        det += d
    for k in ("nglyc", "deamid_NG", "deamid_NS", "isomer_DG", "hydrol_DP",
              "cys_total", "cys_odd", "cys_extra", "trp_total", "trp_cdr",
              "met_total", "met_cdr"):
        rec["pred_liab_%s" % k] = allc.get(k, 0)
    rec["pred_liab_total"] = sum(allc.get(k, 0) for k in
                                 ("nglyc", "deamid_NG", "deamid_NS", "isomer_DG",
                                  "hydrol_DP", "cys_extra", "trp_cdr", "met_cdr"))
    rec["pred_liab_detail"] = "|".join(det)

    # VHH hallmarks + conserved disulfide (sequence level)
    n_hm, hm_detail, fr2 = vhh_hallmarks(H)
    rec["pred_vhh_hallmark_n"] = "" if n_hm is None else n_hm
    rec["pred_vhh_hallmark_detail"] = hm_detail
    rec["pred_fr2_window"] = fr2
    rec["pred_vhh_like"] = "" if n_hm is None else int(n_hm >= 2)
    p1, p2 = conserved_cys(H)
    rec["pred_cys_fr1_pos"] = p1 or ""
    rec["pred_cys_fr3_pos"] = p2 or ""
    ss_seq_ok = bool(p1 and p2)
    rec["pred_canonical_ss_seq_ok"] = int(ss_seq_ok)

    # structure block
    rec["pred_struct_used"] = ""
    rec["pred_hydrophobic_patch_n_res"] = ""
    rec["pred_hydrophobic_patch_sasa_A2"] = ""
    rec["pred_hydrophobic_patch_residues"] = ""
    rec["pred_exposed_apolar_sasa_A2"] = ""
    rec["pred_canonical_ss_struct_closed"] = ""
    rec["pred_struct_n_disulfide"] = ""
    rec["pred_struct_free_cys"] = ""
    rec["pred_struct_note"] = ""
    if struct_path and os.path.exists(struct_path):
        rec["pred_struct_used"] = os.path.basename(struct_path)
        p = structure_patch(struct_path)
        if "error" in p:
            rec["pred_struct_note"] = "patch_failed:" + p["error"]
        else:
            rec["pred_hydrophobic_patch_n_res"] = p["patch_n"]
            rec["pred_hydrophobic_patch_sasa_A2"] = round(p["patch_sasa"], 1)
            rec["pred_hydrophobic_patch_residues"] = p.get("patch_residues", "")
            rec["pred_exposed_apolar_sasa_A2"] = round(p["total_apolar_sasa"], 1)
        d = structure_disulfides(struct_path)
        if "error" not in d:
            rec["pred_struct_n_disulfide"] = d["n_ss"]
            rec["pred_struct_free_cys"] = d["n_free_cys"]
            rec["pred_canonical_ss_struct_closed"] = int(d["n_ss"] >= 1)
            rec["pred_struct_note"] = (rec["pred_struct_note"] + " ss:" + d["ss_detail"]).strip()

    # sequence-only patch proxy (always reported, always labelled)
    mw, mi = max_hydrophobic_window(full, 5)
    rec["pred_hydrophobic_patch_seq_proxy"] = round(mw, 3) if not math.isnan(mw) else ""
    rec["pred_hydrophobic_patch_seq_proxy_start"] = mi
    rec["pred_hydrophobic_patch_seq_proxy_note"] = "PROXY:max mean KD over 5-res window; not A3D"

    # Tm
    if tempro_tm is not None:
        rec["pred_tm_c"] = round(float(tempro_tm), 1)
        rec["pred_tm_source"] = "TEMPRO_keras"
    else:
        rec["pred_tm_c"] = ""
        rec["pred_tm_source"] = "no_tm_predictor_available"
    patch_metric = mw if not math.isnan(mw) else None
    rec["pred_tm_surrogate_rank"] = tm_surrogate_rank(
        H, n_hm, ss_seq_ok, rec["pred_cdrH3_len"], patch_metric, rec["pred_gravy"])
    rec["pred_tm_surrogate_note"] = "RANK 0-100, NOT a temperature, never fitted to Tm data"

    return rec


COLUMNS = [
    "sequence_id", "kind", "len_heavy", "len_light",
    "pred_nonstandard_residues",
    "pred_net_charge_ph74", "pred_pI", "pred_gravy", "pred_mean_hydrophobicity_kd",
    "pred_cdr_parse_ok", "pred_cdrH1", "pred_cdrH2", "pred_cdrH3", "pred_cdrH3_len",
    "pred_cdrL1", "pred_cdrL2", "pred_cdrL3", "pred_cdr_total_len",
    "pred_cdr_frac_arg", "pred_cdr_frac_lys", "pred_cdr_frac_hydrophobic",
    "pred_cdr_net_charge_ph74",
    "pred_liab_nglyc", "pred_liab_deamid_NG", "pred_liab_deamid_NS",
    "pred_liab_isomer_DG", "pred_liab_hydrol_DP", "pred_liab_cys_total",
    "pred_liab_cys_odd", "pred_liab_cys_extra", "pred_liab_trp_total",
    "pred_liab_trp_cdr", "pred_liab_met_total", "pred_liab_met_cdr",
    "pred_liab_total", "pred_liab_detail",
    "pred_vhh_hallmark_n", "pred_vhh_hallmark_detail", "pred_fr2_window",
    "pred_vhh_like", "pred_cys_fr1_pos", "pred_cys_fr3_pos",
    "pred_canonical_ss_seq_ok",
    "pred_struct_used", "pred_hydrophobic_patch_n_res",
    "pred_hydrophobic_patch_sasa_A2", "pred_hydrophobic_patch_residues",
    "pred_exposed_apolar_sasa_A2", "pred_canonical_ss_struct_closed",
    "pred_struct_n_disulfide", "pred_struct_free_cys", "pred_struct_note",
    "pred_hydrophobic_patch_seq_proxy", "pred_hydrophobic_patch_seq_proxy_start",
    "pred_hydrophobic_patch_seq_proxy_note",
    "pred_tm_c", "pred_tm_source", "pred_tm_surrogate_rank", "pred_tm_surrogate_note",
]

HEADER_COMMENT = [
    "# devel.py developability + thermostability PROXY table.",
    "# ALL `pred_` columns are PREDICTED/COMPUTED. No column here is a measured value.",
    "# Producers: physchem = Henderson-Hasselbalch + Kyte-Doolittle (this script);",
    "#   patch = Shrake-Rupley SASA (biotite) + CB-CB 8A single-linkage -- PROXY for",
    "#   Aggrescan3D, not Aggrescan3D; CDRs = Kabat-like regex anchors, NOT ANARCI;",
    "#   pred_tm_c = TEMPRO if --tempro-csv given, else empty (no Tm predictor on this VM);",
    "#   pred_tm_surrogate_rank = hand-weighted ordinal rank, NEVER a temperature.",
]


def find_structure(sid, struct_dir):
    if not struct_dir:
        return None
    for ext in (".pdb", ".cif", ".mmcif", ".pdb.gz"):
        p = os.path.join(struct_dir, sid + ext)
        if os.path.exists(p):
            return p
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--structure", default=None)
    ap.add_argument("--structure-dir", default=None)
    ap.add_argument("--tempro-csv", default=None)
    ap.add_argument("--no-header-comment", action="store_true")
    args = ap.parse_args()

    tempro = {}
    if args.tempro_csv:
        with open(args.tempro_csv, newline="") as fh:
            for r in csv.DictReader(fh):
                tempro[r["sequence_id"].strip()] = r.get("tm_c")

    with open(args.csv, newline="") as fh:
        rows = [r for r in csv.DictReader(
            (ln for ln in fh if not ln.lstrip().startswith("#")))]
    if not rows:
        sys.exit("no rows in %s" % args.csv)

    recs = []
    for r in rows:
        sid = (r.get("sequence_id") or "").strip()
        sp = args.structure if (args.structure and len(rows) == 1) else \
            find_structure(sid, args.structure_dir)
        recs.append(score_design(r, sp, tempro.get(sid)))

    with open(args.out, "w", newline="") as fh:
        if not args.no_header_comment:
            for ln in HEADER_COMMENT:
                fh.write(ln + "\n")
        w = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction="ignore")
        w.writeheader()
        for rec in recs:
            w.writerow(rec)
    print("wrote %s (%d designs, %d columns)" % (args.out, len(recs), len(COLUMNS)))
    for rec in recs:
        print("  %-24s pI=%-6s q7.4=%-7s gravy=%-7s cdrH3=%-2s argCDR=%-6s "
              "liab=%-3s hm=%-2s ss=%s patch=%s tm=%s rank=%s" % (
                  rec["sequence_id"], rec["pred_pI"], rec["pred_net_charge_ph74"],
                  rec["pred_gravy"], rec["pred_cdrH3_len"], rec["pred_cdr_frac_arg"],
                  rec["pred_liab_total"], rec["pred_vhh_hallmark_n"],
                  rec["pred_canonical_ss_seq_ok"],
                  rec["pred_hydrophobic_patch_sasa_A2"] or "n/a",
                  rec["pred_tm_source"], rec["pred_tm_surrogate_rank"]))


if __name__ == "__main__":
    main()
