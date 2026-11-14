#!/usr/bin/env python3
"""Build the RFD4-Proteina inputs for de novo VHH design on the trastuzumab epitope.

Goal: a new VHH (nanobody) that binds HER2 domain IV on the trastuzumab epitope, so that it
competes with trastuzumab in the epitope-preservation assay.

Reads (never writes) two organiser files:
  apheris/rfd4/Data/12-her2/her2_fv_1n8z.cif        PDB 1N8Z cropped to Fv + HER2 domain IV
  apheris/rfd4/Data/03-vhh/pmvk_caplacizumab.cif.zst  caplacizumab VHH framework, CDRs collapsed

Writes into ./spec/ (the spec and the structures it references must sit side by side, and the
spec FILE, not its folder, is what `module-run design --input` takes):
  target_her2_d4.cif     HER2 chains C (489-580) and D (591-607), antibody chains removed
  framework_vhh.cif      VHH chain H with the three CDR slots (23, 38, 76) left empty
  design_fixedfw.json    framework sequence fixed (organiser-documented pattern)
  design_freefw.json     same, but exposed framework positions are redesigned too
  inputs_manifest.json   sources, hashes, thresholds and the freed positions

Why two variants: with a fixed framework the whole-domain identity to caplacizumab has a floor
of n_fixed / total_length, which sits right at the track's 70% line (see NOTES.md). Redesigning
exposed framework positions buys margin below it.

Pure standard library. Reads the .zst with the `zstandard` module if present, else `zstd -dc`.
Run from anywhere:  python3 build_inputs.py
"""
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]  # LEO/tasks/<task>/ -> repo root
DATA = REPO / "apheris" / "rfd4" / "Data"
HER2_CIF = DATA / "12-her2" / "her2_fv_1n8z.cif"
VHH_ZST = DATA / "03-vhh" / "pmvk_caplacizumab.cif.zst"
OUT = HERE / "spec"

TARGET_CHAINS = ("C", "D")  # HER2 split at its unresolved loop 581-590 (organiser's convention)
VHH_CHAIN = "H"

# CDR length ranges are the ones the organisers encode in the placeholder residues of their own
# VHH example (comp_id "13-13", "9-10", "12-21" at slots 23, 38, 76).
CDR_LENGTHS = ("13-13", "9-10", "12-21")
<<<<<<< HEAD
LONG_CDR_LENGTHS = ("13-13", "9-10", "17-21")  # the fixedfw_long variant: CDR3 17-21 only
=======
>>>>>>> origin/main

# Hotspots: HER2 atoms the VHH should contact. These are the four atoms the organisers use in
# Data/12-her2/her2_cdrh3.yaml ("the HER2 atoms nearest CDRH3 in PDB 1N8Z"), same atoms, same order.
# They follow the RFD4-Proteina definition (a target atom within 4.5 A of the binder, here trastuzumab's
# CDRH3), and check_spec.py re-verifies that against the structure. See NOTES.md.
# (chain, residue, expected comp_id, atom)
HOTSPOTS = (
    ("D", 593, "LYS", "NZ"),
    ("C", 573, "PHE", "CD1"),
    ("D", 591, "ILE", "CG2"),
    ("C", 572, "PRO", "CD"),
)

# Rule for which framework positions the free variant redesigns.
BURIAL_RADIUS = 10.0  # A, CB-CB (CA for Gly) neighbour count
FREE_MAX_NEIGHBOURS = 12  # exposed: at most this many neighbours
CDR_MARGIN = 3  # keep residues within this many positions of a CDR slot (loop anchors)
NTERM_KEEP = 3  # keep the first residues (chain-end, count is unreliable)
FR4_START = 77  # keep the conserved J-region fixed

AA3 = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
       "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}

# Column order of the organisers' framework/target files.
COLS = ("group_PDB", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
        "label_asym_id", "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code", "auth_seq_id",
        "auth_comp_id", "auth_asym_id", "auth_atom_id", "Cartn_x", "Cartn_y", "Cartn_z",
        "pdbx_PDB_model_num", "id")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_zst(path):
    try:
        import zstandard
        with open(path, "rb") as f:
            return zstandard.ZstdDecompressor().stream_reader(f).read().decode("utf8")
    except ImportError:
        return subprocess.run(["zstd", "-dc", str(path)], check=True, capture_output=True).stdout.decode("utf8")


def parse_atom_site(text):
    """Return a list of dicts, one per atom_site row, keyed by the column names in the file."""
    lines = text.split("\n")
    i = 0
    while i < len(lines):
        if lines[i].strip() == "loop_" and i + 1 < len(lines) and lines[i + 1].startswith("_atom_site."):
            j = i + 1
            cols = []
            while lines[j].lstrip().startswith("_"):
                cols.append(lines[j].strip().split(".", 1)[1])
                j += 1
            atoms = []
            while j < len(lines):
                s = lines[j].strip()
                if s in ("", "#") or s.startswith("_") or s == "loop_":
                    break
                tok = s.split()
                if len(tok) != len(cols):
                    raise ValueError(f"row has {len(tok)} fields, expected {len(cols)}: {s[:80]}")
                atoms.append(dict(zip(cols, tok)))
                j += 1
            return atoms
        i += 1
    raise ValueError("no _atom_site loop found")


def cell_block(text):
    return [l for l in text.split("\n") if l.startswith("_cell.")]


def write_cif(path, atoms, data_name, cell_lines=None):
    out = [f"data_{data_name}", "#", "loop_"] + [f"_atom_site.{c}" for c in COLS]
    for n, a in enumerate(atoms, 1):
        alt = "." if a["label_alt_id"] in ("?", ".") else a["label_alt_id"]
        if alt != ".":
            raise ValueError(f"alternate location {alt} on {a['label_asym_id']}{a['label_seq_id']}; resolve upstream")
        ins = "."  # none of these residues carries an insertion code
        out.append(" ".join([
            "ATOM", a["type_symbol"], a["label_atom_id"], ".", a["label_comp_id"], a["label_asym_id"],
            a["label_entity_id"], a["label_seq_id"], ins, a["auth_seq_id"], a["auth_comp_id"],
            a["auth_asym_id"], a["auth_atom_id"], a["Cartn_x"], a["Cartn_y"], a["Cartn_z"], "1", str(n)]))
    out.append("#")
    if cell_lines:
        out += cell_lines + ["#"]
    Path(path).write_text("\n".join(out) + "\n", encoding="utf8", newline="\n")


def residues(atoms, chain):
    """{label_seq_id: {"comp": ..., "atoms": {name: (x, y, z)}}} for one chain, in file order."""
    res = {}
    for a in atoms:
        if a["label_asym_id"] != chain:
            continue
        r = res.setdefault(int(a["label_seq_id"]), {"comp": a["label_comp_id"], "atoms": {}})
        r["atoms"][a["label_atom_id"]] = (float(a["Cartn_x"]), float(a["Cartn_y"]), float(a["Cartn_z"]))
    return res


def ranges(positions):
    """[1,2,3,5,7,8] -> '1-3, 5, 7-8' (compact selection syntax)."""
    pos = sorted(positions)
    out, start, prev = [], None, None
    for p in pos:
        if start is None:
            start = prev = p
        elif p == prev + 1:
            prev = p
        else:
            out.append((start, prev))
            start = prev = p
    if start is not None:
        out.append((start, prev))
    return ", ".join(f"{a}-{b}" if a != b else f"{a}" for a, b in out)


def free_positions(fw):
    """Exposed framework positions the free variant redesigns. Deterministic, see module docstring."""
    pts = {}
    for s, r in fw.items():
        pts[s] = r["atoms"].get("CB") or r["atoms"]["CA"]
    nbrs = {s: sum(1 for t, q in pts.items() if t != s and math.dist(p, q) <= BURIAL_RADIUS)
            for s, p in pts.items()}
    slots = [s for s in range(min(fw), max(fw) + 1) if s not in fw]  # the three CDR slots
    free = []
    for s, r in sorted(fw.items()):
        if nbrs[s] > FREE_MAX_NEIGHBOURS:
            continue
        if r["comp"] in ("CYS", "GLY", "PRO"):
            continue
        if any(abs(s - c) <= CDR_MARGIN for c in slots):
            continue
        if s <= NTERM_KEEP or s >= FR4_START:
            continue
        free.append(s)
    return free, nbrs


def sel(s):
    return {"select": f"sel('{s}')"}


<<<<<<< HEAD
def make_spec(target_ranges, fw_positions, free, fw_segments, cdr_lengths=CDR_LENGTHS):
=======
def make_spec(target_ranges, fw_positions, free, fw_segments):
>>>>>>> origin/main
    fixed = [p for p in fw_positions if p not in set(free)]
    contig = [f"C{target_ranges['C'][0]}-{target_ranges['C'][1]}", "/0",
              f"D{target_ranges['D'][0]}-{target_ranges['D'][1]}", "/0"]
    parts = []
    for k, (a, b) in enumerate(fw_segments):
        parts.append(f"{VHH_CHAIN}{a}-{b}")
<<<<<<< HEAD
        if k < len(cdr_lengths):
            parts.append(cdr_lengths[k])
=======
        if k < len(CDR_LENGTHS):
            parts.append(CDR_LENGTHS[k])
>>>>>>> origin/main
    contig_str = ", ".join(contig + parts)
    bb = "[N,CA,C,O]"
    return {
        "contig": contig_str,
        "conditions": {
            "C_CRD": {"true": [sel(f"C/*/*/{bb}"), sel(f"D/*/*/{bb}")]},
            "C_CRDGRP": {"1": [sel(f"{VHH_CHAIN}/*/[{ranges(fw_positions)}]/{bb}")]},
            "C_SEQ": {"true": [sel("C"), sel("D"), sel(f"{VHH_CHAIN}/*/[{ranges(fixed)}]")]},
            "C_HOT": {"true": [sel(f"{c}/*/{n}/{atom}") for c, n, _, atom in HOTSPOTS]},
            "C_CATH": {"2.60.40": [sel(VHH_CHAIN)]},
        },
    }


def main():
    for p in (HER2_CIF, VHH_ZST):
        if not p.is_file():
            sys.exit(f"missing organiser input: {p}")
    her2_text = HER2_CIF.read_text(encoding="utf8")
    her2 = parse_atom_site(her2_text)
    vhh = parse_atom_site(read_zst(VHH_ZST))

    # target: HER2 domain IV only
    target = [a for a in her2 if a["label_asym_id"] in TARGET_CHAINS]
    tres = {c: residues(target, c) for c in TARGET_CHAINS}
    for c, n, comp, atom in HOTSPOTS:
        got = tres[c].get(n)
        if got is None or got["comp"] != comp or atom not in got["atoms"]:
            sys.exit(f"hotspot {c}/{n}/{atom} ({comp}) not found in the target; refusing to continue")
    target_ranges = {c: (min(r), max(r)) for c, r in tres.items()}

    # framework: chain H without CDR placeholder residues (comp_id like '13-13')
    fw_atoms = [a for a in vhh if a["label_asym_id"] == VHH_CHAIN and a["label_comp_id"] in AA3]
    fw = residues(fw_atoms, VHH_CHAIN)
    positions = sorted(fw)
    slots = [s for s in range(positions[0], positions[-1] + 1) if s not in fw]
    if len(slots) != len(CDR_LENGTHS):
        sys.exit(f"expected {len(CDR_LENGTHS)} CDR slots in the framework, found {slots}")
    segments, start = [], positions[0]
    for s in slots:
        segments.append((start, s - 1))
        start = s + 1
    segments.append((start, positions[-1]))

    free, nbrs = free_positions(fw)

    OUT.mkdir(parents=True, exist_ok=True)
    write_cif(OUT / "target_her2_d4.cif", target, "structure", cell_block(her2_text))
    write_cif(OUT / "framework_vhh.cif", fw_atoms, "framework")

    specs = {
        "design_fixedfw.json": {"her2d4_vhh_fixedfw": {"input": ["target_her2_d4.cif", "framework_vhh.cif"],
                                                        **make_spec(target_ranges, positions, [], segments)}},
        "design_freefw.json": {"her2d4_vhh_freefw": {"input": ["target_her2_d4.cif", "framework_vhh.cif"],
                                                      **make_spec(target_ranges, positions, free, segments)}},
    }
<<<<<<< HEAD
    # Third variant, added after scoring batch 1: fixed framework, CDR3 constrained to 17-21. In batch 1 every
    # design with a CDR3 of 15 or fewer residues carried a HIGH estimated track risk (framework share of the
    # whole VHH >= 70%), while 17-21 had the same epitope quality (see NOTES.md section 8).
    specs["design_fixedfw_long.json"] = {"her2d4_vhh_fixedfw_long": {
        "input": ["target_her2_d4.cif", "framework_vhh.cif"],
        **make_spec(target_ranges, positions, [], segments, cdr_lengths=LONG_CDR_LENGTHS)}}
=======
>>>>>>> origin/main
    for name, spec in specs.items():
        (OUT / name).write_text(json.dumps(spec, indent=2) + "\n", encoding="utf8", newline="\n")

    manifest = {
        "sources": {str(HER2_CIF.relative_to(REPO)): sha256(HER2_CIF), str(VHH_ZST.relative_to(REPO)): sha256(VHH_ZST)},
        "outputs": {f.name: sha256(f) for f in sorted(OUT.glob("*")) if f.name != "inputs_manifest.json"},
        "target_ranges": target_ranges,
        "framework_residues": len(fw),
        "framework_segments": segments,
        "cdr_slots": slots,
        "cdr_lengths": list(CDR_LENGTHS),
        "hotspots": [f"{c}/{n}/{atom}" for c, n, _, atom in HOTSPOTS],
        "free_rule": {"burial_radius_A": BURIAL_RADIUS, "max_neighbours": FREE_MAX_NEIGHBOURS,
                      "cdr_margin": CDR_MARGIN, "keep_first": NTERM_KEEP, "keep_from": FR4_START,
                      "never_free": ["CYS", "GLY", "PRO"]},
        "free_positions": free,
    }
    (OUT / "inputs_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf8", newline="\n")

    print(f"target: C {target_ranges['C']}, D {target_ranges['D']}  ({len(target)} atoms)")
    print(f"framework: {len(fw)} residues, segments {segments}, CDR slots {slots}")
    print(f"free positions ({len(free)}): " + ", ".join(f"{fw[s]['comp'].title()}{s}(n={nbrs[s]})" for s in free))
    for name in specs:
        print("wrote", OUT / name)


if __name__ == "__main__":
    main()
