#!/usr/bin/env python3
"""specgen.py — build RFD4-Proteina framework CIFs and design specs for the
de novo HER2 campaign.

What was reverse-engineered to make this work
---------------------------------------------
The Apheris docs show a spec but not how the framework file must be shaped.
Decoded from the two bundled examples:

`Data/03-fv/framework_adalimumab.cif` has chain H numbered 1-22, 24-37, 39-75,
77-87 -- the CDRs are DELETED and a single skipped residue number marks each
insertion point. Its contig is
    "F1-76, /0, H1-22, 13-15, H24-37, 9-10, H39-75, 11-22, H77-87, /0, L1-23, ..."
so a bare range in the contig (13-15) is the designed CDR length, and it has no
relation to the size of the numbering gap (which is 1).

`Data/03-vhh/pmvk_caplacizumab.cif.zst` is the same thing in *annotated* form:
chain H is numbered 1-86 with no gaps, and positions 23, 38 and 76 are sentinel
residues whose `label_comp_id` is literally the length range -- `13-13`, `9-10`,
`12-21` -- carried by the `annotation_expsegmin_1_residue` /
`annotation_expsegmax_1_residue` columns. The RFD4 log calls these "expandable
segment sentinels". So caplacizumab's framework is 83 real residues at
1-22 / 24-37 / 39-75 / 77-86, and the bundle's own CDR lengths are
CDR1 = 13, CDR2 = 9-10, CDR3 = 12-21.

Both files therefore use ONE convention, and this module implements it for the
plain (non-annotated) route: drop the sentinels, keep the real framework
residues at their original numbering, and let the contig state each CDR length.
We use the plain route rather than hand-writing annotated CIFs because the
annotated format carries ~11 interdependent condition columns per atom, and a
silently wrong one is a failure we would not detect until after the GPU spend.

The identity lever
------------------
A design that keeps its scaffold's framework sequence is ~61-69% identical to
that scaffold's parent antibody by construction. Dropping the framework from
`C_SEQ` (keeping it only in `C_CRDGRP`, i.e. backbone-only) lets the model
redesign framework sequence while holding the fold, which is our main lever for
pushing identity down. `framework_mode` selects this:

  seqfix  framework sequence fixed (C_SEQ covers antigen + framework)
  bbonly  framework sequence free  (C_SEQ covers antigen only)

Usage
-----
  specgen.py framework --src <cif> --chain H --out <cif>      # strip sentinels
  specgen.py spec --arm-file arms.yaml --outdir specs/
  specgen.py arms --out arms.json                             # emit arm matrix
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

# Plain mmCIF atom_site column order we emit. RFD4 rejects PDB; it accepts
# ordinary mmCIF, and this is the exact column set the bundled plain examples
# (framework_adalimumab.cif, target_F.cif) use.
OUT_COLS = [
    "group_PDB", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
    "label_asym_id", "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code",
    "auth_seq_id", "auth_comp_id", "auth_asym_id", "auth_atom_id",
    "Cartn_x", "Cartn_y", "Cartn_z", "pdbx_PDB_model_num", "id",
]


def read_atom_site(path: Path) -> tuple[list[str], list[list[str]]]:
    """Read an mmCIF's atom_site loop as (column names, rows).

    Deliberately a small hand parser: the bundled files are single-block,
    single-model, whitespace-delimited mmCIF, and this avoids a dependency that
    may not be importable in every venv on this VM.
    """
    text = path.read_text()
    cols: list[str] = []
    rows: list[list[str]] = []
    in_loop = False
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("_atom_site."):
            cols.append(s.split(".", 1)[1].split()[0])
            in_loop = True
            continue
        if in_loop and (s.startswith("ATOM") or s.startswith("HETATM")):
            f = s.split()
            if len(f) >= len(cols):
                rows.append(f)
            continue
        if in_loop and rows and (s.startswith("#") or s.startswith("loop_") or s.startswith("data_")):
            break
    if not cols:
        raise SystemExit(f"{path}: no _atom_site loop found")
    return cols, rows


def write_plain_cif(out: Path, entries: list[dict], entry_id: str) -> None:
    """Write a minimal, valid plain mmCIF with just the atoms we kept."""
    lines = [f"data_{entry_id}", "#", f"_entry.id  {entry_id}", "#", "loop_"]
    lines += [f"_atom_site.{c}" for c in OUT_COLS]
    for i, e in enumerate(entries, start=1):
        e = dict(e)
        e["id"] = str(i)
        lines.append(" ".join(str(e.get(c, ".")) for c in OUT_COLS))
    lines.append("#")
    out.write_text("\n".join(lines) + "\n")


def cmd_framework(a) -> int:
    """Extract one chain's REAL framework residues into a plain mmCIF.

    Drops expandable-segment sentinel residues (whose label_comp_id is a length
    range like `12-21`, not an amino acid) and any non-standard residue, keeping
    original numbering so the natural gaps mark the CDR insertion points.
    Reports the sentinel lengths it found, because those are the bundle's own
    CDR length priors and we want them in the ledger rather than guessed.
    """
    src = Path(a.src)
    cols, rows = read_atom_site(src)
    ix = {c: i for i, c in enumerate(cols)}
    need = ["label_asym_id", "label_seq_id", "label_comp_id", "label_atom_id",
            "Cartn_x", "Cartn_y", "Cartn_z", "type_symbol"]
    missing = [c for c in need if c not in ix]
    if missing:
        raise SystemExit(f"{src}: atom_site lacks {missing}")

    kept: list[dict] = []
    sentinels: dict[int, str] = {}
    seen: dict[int, str] = {}
    for f in rows:
        if f[ix["label_asym_id"]] != a.chain:
            continue
        try:
            n = int(f[ix["label_seq_id"]])
        except ValueError:
            continue
        comp = f[ix["label_comp_id"]]
        if comp not in AA3:
            # a length-range sentinel (`13-13`) or an unknown residue
            sentinels[n] = comp
            continue
        if a.keep_backbone_only and f[ix["label_atom_id"]] not in ("N", "CA", "C", "O"):
            continue
        seen[n] = comp
        kept.append({
            "group_PDB": "ATOM",
            "type_symbol": f[ix["type_symbol"]],
            "label_atom_id": f[ix["label_atom_id"]],
            "label_alt_id": ".",
            "label_comp_id": comp,
            "label_asym_id": a.out_chain or a.chain,
            "label_entity_id": "1",
            "label_seq_id": str(n),
            "pdbx_PDB_ins_code": "?",
            "auth_seq_id": str(n),
            "auth_comp_id": comp,
            "auth_asym_id": a.out_chain or a.chain,
            "auth_atom_id": f[ix["label_atom_id"]],
            "Cartn_x": f[ix["Cartn_x"]],
            "Cartn_y": f[ix["Cartn_y"]],
            "Cartn_z": f[ix["Cartn_z"]],
            "pdbx_PDB_model_num": "1",
        })

    if not kept:
        raise SystemExit(f"{src}: no standard residues found in chain {a.chain}")

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_plain_cif(out, kept, a.entry_id or out.stem)

    ks = sorted(seen)
    seq = "".join(AA3[seen[i]] for i in ks)
    segs = contiguous_segments(ks)
    report = {
        "src": str(src), "out": str(out),
        "chain_in": a.chain, "chain_out": a.out_chain or a.chain,
        "n_residues": len(ks), "n_atoms": len(kept),
        "numbering": [ks[0], ks[-1]],
        "kept_segments": [f"{s}-{e}" for s, e in segs],
        "gaps": [i for i in range(ks[0], ks[-1]) if i not in seen],
        "sentinels_found": {str(k): v for k, v in sorted(sentinels.items())},
        "sequence": seq,
    }
    print(json.dumps(report, indent=2))
    return 0


def contiguous_segments(nums: list[int]) -> list[tuple[int, int]]:
    segs, start, prev = [], nums[0], nums[0]
    for n in nums[1:]:
        if n != prev + 1:
            segs.append((start, prev))
            start = n
        prev = n
    segs.append((start, prev))
    return segs


def build_contig(antigen_segments: list[str], fw_segments: list[tuple[str, list[tuple[int, int]]]],
                 cdr_lengths: list[str]) -> str:
    """Compose a contig string.

    antigen_segments: ["C489-580", "D591-607"] -- each becomes its own chain
    fw_segments: [("H", [(1,22),(24,37),(39,75),(77,86)]), ...]
    cdr_lengths: ["13-13", "9-10", "12-21"] interleaved between the kept framework
                 segments of each chain (len(segments)-1 entries per chain)
    """
    parts: list[str] = []
    for seg in antigen_segments:
        if parts:
            parts.append("/0")
        parts.append(seg)
    li = 0
    for chain, segs in fw_segments:
        parts.append("/0")
        for k, (s, e) in enumerate(segs):
            parts.append(f"{chain}{s}-{e}")
            if k < len(segs) - 1:
                if li >= len(cdr_lengths):
                    raise SystemExit("not enough cdr_lengths for the framework segments")
                parts.append(cdr_lengths[li])
                li += 1
    # "/0" is a chain break token; join with ", " exactly as the bundled spec does
    out = []
    for p in parts:
        out.append(p)
    return ", ".join(out)


def build_spec(design_id: str, inputs: list[str], contig: str,
               antigen_sel: list[str], fw_sel: list[str], hotspots: list[str],
               ab_chains: list[str], framework_mode: str) -> dict:
    """Assemble one RFD4 design spec.

    C_CRD   absolute coordinates on the antigen backbone. Needs >= 3 absolute
            atoms or `pass` and target-aligned RMSD come back EMPTY -- that is
            the exact failure already seen in runs/her2/scored-boltz2, so the
            antigen backbone is always pinned here.
    C_CRDGRP 1 on framework backbones: internal geometry preserved, rigid-body
            placement free. One group so a paired VH/VL stays paired.
    C_SEQ   antigen always. Framework only in `seqfix` mode; omitting it in
            `bbonly` mode is what lets the model lower identity by redesigning
            framework sequence.
    C_HOT   the epitope atoms we want contacted.
    C_CATH  2.60.40 immunoglobulin fold on the antibody chains.
    """
    conds: dict = {
        "C_CRD": {"true": [{"select": f"sel('{s}') and not has_nan_coord()"} for s in antigen_sel]},
        "C_CRDGRP": {"1": [{"select": f"sel('{s}') and not has_nan_coord()"} for s in fw_sel]},
    }
    seq_sel = list(antigen_sel_plain(antigen_sel))
    if framework_mode == "seqfix":
        seq_sel += list(fw_sel_plain(fw_sel))
    conds["C_SEQ"] = {"true": [{"select": f"sel('{s}')"} for s in seq_sel]}
    if hotspots:
        conds["C_HOT"] = {"true": [{"select": f"sel('{h}')"} for h in hotspots]}
    conds["C_CATH"] = {"2.60.40": [{"select": f"sel('{c}')"} for c in ab_chains]}
    return {design_id: {"input": inputs, "contig": contig, "conditions": conds}}


def antigen_sel_plain(sels: list[str]):
    """C_CRD selectors carry an atom filter (/[N,CA,C,O]); C_SEQ must not."""
    for s in sels:
        yield s.split("/[")[0].rstrip("/") if "/[" in s else s


def fw_sel_plain(sels: list[str]):
    for s in sels:
        if "/[N,CA,C,O]" in s:
            yield s.replace("/[N,CA,C,O]", "")
        else:
            yield s


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("framework", help="extract a framework chain to a plain mmCIF")
    s.add_argument("--src", required=True)
    s.add_argument("--chain", required=True)
    s.add_argument("--out-chain")
    s.add_argument("--out", required=True)
    s.add_argument("--entry-id")
    s.add_argument("--keep-backbone-only", action="store_true")
    s.set_defaults(fn=cmd_framework)

    a = ap.parse_args()
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
