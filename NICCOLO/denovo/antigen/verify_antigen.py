"""Re-parse the antigen-only CIF and assert it is a valid RFD4 design input."""
import numpy as np, sys
sys.path.insert(0, "/home/hackathon/hackathon-team-3/NICCOLO/denovo/antigen")
ANT = "/home/hackathon/hackathon-team-3/NICCOLO/denovo/antigen/her2_dIV_CD.cif"
SRC = "/home/hackathon/hackathon-team-3/apheris/rfd4/Data/12-her2/her2_fv_1n8z.cif"

def raw_atom_site(path):
    cols, in_loop, rows = [], False, []
    for line in open(path):
        if line.startswith("_atom_site."):
            cols.append(line.strip().split(".", 1)[1]); in_loop = True; continue
        if in_loop and not cols: continue
        if in_loop and (line.startswith("#") or line.startswith("loop_")): break
        if not (in_loop and line.startswith(("ATOM", "HETATM"))): continue
        t = line.split()
        assert len(t) == len(cols), f"{len(t)} vs {len(cols)}: {line}"
        rows.append(dict(zip(cols, t)))
    return cols, rows

cols, rows = raw_atom_site(ANT)
print(f"atom_site columns ({len(cols)}): {cols}")
print(f"atom_site rows: {len(rows)}  (incl. any hydrogens/HETATM)")

# 1) no A/B atoms survive, by every chain-bearing column
for col in ("label_asym_id", "auth_asym_id"):
    bad = sorted({r[col] for r in rows} - {"C", "D"})
    print(f"  chains in {col}: {sorted({r[col] for r in rows})}  unexpected={bad}")
    assert not bad
badent = sorted({r["label_entity_id"] for r in rows} - {"3", "4"})
print(f"  label_entity_id present: {sorted({r['label_entity_id'] for r in rows})} unexpected={badent}")
assert not badent
assert "A" not in open(ANT).read().split("loop_\n_atom_site")[0] or True  # header chains checked below

# header must not advertise trastuzumab any more
hdr = open(ANT).read().split("_atom_site.group_PDB")[0]
for s in ("DIQMTQSPSS", "EVQLVESGGG", "1 polymer", "2 polymer", "1 polypeptide", "2 polypeptide"):
    assert s not in hdr, f"trastuzumab metadata survived in header: {s}"
print("  header: no trastuzumab entity / sequence rows remain  OK")
assert " A ?" not in hdr and " B ?" not in hdr
print("  header: entity_poly strand ids A/B gone  OK")

# 2) per-chain residue inventory
inv = {}
for r in rows:
    inv.setdefault(r["label_asym_id"], {})[int(r["label_seq_id"])] = r["label_comp_id"]
for c in sorted(inv):
    n = sorted(inv[c])
    gaps = [(a, b) for a, b in zip(n, n[1:]) if b - a != 1]
    print(f"  chain {c}: {len(n)} residues, range {n[0]}-{n[-1]}, "
          f"contiguous={not gaps} gaps={gaps}")
assert len(inv["C"]) == 92 and min(inv["C"]) == 489 and max(inv["C"]) == 580
assert len(inv["D"]) == 17 and min(inv["D"]) == 591 and max(inv["D"]) == 607

# 3) label_seq_id == auth_seq_id everywhere
assert all(r["label_seq_id"] == r["auth_seq_id"] for r in rows)
print("  label_seq_id == auth_seq_id for all rows  OK")

# 4) coordinates: finite, no NaN, no sentinel
X = np.array([[float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"])] for r in rows])
print(f"  coords finite={np.isfinite(X).all()}  nan={int(np.isnan(X).sum())}  "
      f"bbox min={X.min(0).round(2).tolist()} max={X.max(0).round(2).tolist()}")
assert np.isfinite(X).all()

# 5) ids contiguous
ids = [int(r["id"]) for r in rows]
aai = [int(r["atom_array_index"]) for r in rows]
assert ids == list(range(1, len(rows) + 1)) and aai == list(range(len(rows)))
print("  _atom_site.id 1..N and atom_array_index 0..N-1 contiguous  OK")

# 6) coordinates identical to the source for every kept atom (no silent perturbation)
_, srows = raw_atom_site(SRC)
skey = {(r["label_asym_id"], r["label_seq_id"], r["label_atom_id"]):
        (r["Cartn_x"], r["Cartn_y"], r["Cartn_z"]) for r in srows}
nmis = sum(1 for r in rows if skey[(r["label_asym_id"], r["label_seq_id"], r["label_atom_id"])]
           != (r["Cartn_x"], r["Cartn_y"], r["Cartn_z"]))
print(f"  coordinate strings differing from source: {nmis}")
assert nmis == 0
nsrc_cd = sum(1 for r in srows if r["label_asym_id"] in ("C", "D"))
print(f"  source C/D atom rows={nsrc_cd}  kept={len(rows)}  all kept={nsrc_cd == len(rows)}")
assert nsrc_cd == len(rows)

# 7) independent parse with biotite (the kind of reader a validator uses)
import biotite.structure.io.pdbx as pdbx
import biotite.structure as struc
cif = pdbx.CIFFile.read(ANT)
arr = pdbx.get_structure(cif, model=1, use_author_fields=False)
print(f"  biotite parse: {arr.array_length()} atoms, chains={sorted(set(arr.chain_id))}, "
      f"{struc.get_residue_count(arr)} residues, nan_coord={int(np.isnan(arr.coord).sum())}")
assert sorted(set(arr.chain_id)) == ["C", "D"]
assert struc.get_residue_count(arr) == 109

# 8) hotspot atom existence check
present = {(r["label_asym_id"], int(r["label_seq_id"]), r["label_atom_id"]) for r in rows}
SETS = {
 # A: control - reproduce Apheris's four CDRH3-proximal atoms verbatim
 "A_control_cdrh3":    ["D/593/NZ", "C/573/CD1", "D/591/CG2", "C/572/CD"],
 # B: burial-ranked O-ring core of the VH-side sub-patch
 "B_buried_core_VH":   ["C/558/OE1", "C/560/OD2", "C/571/CB", "C/572/CB", "C/573/CZ"],
 # C: VL-side juxtamembrane sub-patch, zero CDRH3 contact
 "C_juxtamembrane_VL": ["D/596/OD2", "D/598/OE1", "D/600/CB", "D/602/OE1", "D/603/CG"],
}
print("\nhotspot atom existence in her2_dIV_CD.cif:")
ok = True
for name, sel in SETS.items():
    for s in sel:
        c, n, a = s.split("/")
        hit = (c, int(n), a) in present
        ok &= hit
        print(f"  {name:20s} {s:14s} present={hit}")
print(f"\nALL CHECKS PASSED={ok}")
