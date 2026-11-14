#!/usr/bin/env python3
"""Check a design spec against its structures BEFORE spending GPU time.

RFD4 selections fail silently when they name the wrong atoms, so this resolves every selector
in the spec against the CIFs it references and verifies what the design intends:

  1. referenced files exist, chain IDs are distinct across files
  2. contig: retained ranges are exactly the residues present, CDR lengths parse
  3. every selector resolves to atoms; fixed-geometry selectors are backbone-only
  4. C_CRD covers all of HER2 domain IV, C_CRDGRP covers the whole framework backbone
  5. C_SEQ keeps both framework cysteines; freed positions are never Cys/Gly/Pro
  6. every hotspot is a single real atom, in contact with trastuzumab in PDB 1N8Z (within 4.5 A of an
     antibody heavy atom, and of its CDRH3 specifically: the organisers' definition), with unchanged coordinates
  7. hotspots are compact enough for one VHH paratope
Then prints the identity arithmetic for the track's 70% line (informational).

    python3 check_spec.py spec/design_fixedfw.json [spec/design_freefw.json ...]

Exit status is non-zero if any check fails.
"""
import json
import math
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from build_inputs import HER2_CIF, parse_atom_site  # noqa: E402  (file parsing only, not the logic under test)

BACKBONE = {"N", "CA", "C", "O"}
CONTACT_CUTOFF = 4.5  # A, same cutoff as the interface measurement
CDRH3 = {("B", k) for k in range(97, 110)}  # trastuzumab CDRH3 (SRWGGDGFYAMDY) in the organisers' file
MAX_HOTSPOT_SPAN = 20.0  # A; a VHH paratope spans roughly 18-22 A, so hotspots must fit inside one
SEL_RE = re.compile(r"^sel\('([^']+)'\)$")

failures = []


def check(ok, label, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f"  -- {detail}" if detail else ""))
    if not ok:
        failures.append(label)


def load_structure(path):
    """{(chain, resnum): {"comp": str, "atoms": {name: (x, y, z)}}}"""
    out = {}
    for a in parse_atom_site(Path(path).read_text(encoding="utf8")):
        k = (a["label_asym_id"], int(a["label_seq_id"]))
        r = out.setdefault(k, {"comp": a["label_comp_id"], "atoms": {}})
        r["atoms"][a["label_atom_id"]] = (float(a["Cartn_x"]), float(a["Cartn_y"]), float(a["Cartn_z"]))
    return out


def parse_res(spec):
    if spec == "*":
        return None
    out = set()
    for part in spec.strip("[]").split(","):
        part = part.strip()
        if "-" in part:
            a, b = part.split("-")
            out.update(range(int(a), int(b) + 1))
        else:
            out.add(int(part))
    return out


def resolve(select, table):
    """Atoms (chain, res, name) named by one sel('...') string; unknown syntax raises."""
    m = SEL_RE.match(select)
    if not m:
        raise ValueError(f"unrecognised selector {select!r}")
    parts = m.group(1).split("/")
    chain = parts[0]
    res = parse_res(parts[2]) if len(parts) > 2 else None
    atoms = parts[3].strip("[]").replace(" ", "").split(",") if len(parts) > 3 and parts[3] != "*" else None
    hits = []
    for (c, n), r in table.items():
        if c != chain or (res is not None and n not in res):
            continue
        for name in r["atoms"]:
            if atoms is None or name in atoms:
                hits.append((c, n, name))
    return hits


def selects(cond, key):
    """All selector strings under one condition, across its value groups."""
    return [e["select"] for group in cond.get(key, {}).values() for e in group]


def epitope_contacts(complex_, atom_key, only=None):
    """Closest trastuzumab heavy atom to one HER2 atom in the 1N8Z complex: (distance, 'B:TYR33').
    `only` restricts the trastuzumab residues considered, as a set of (chain, resnum)."""
    p = complex_[atom_key[:2]]["atoms"][atom_key[2]]
    best = (1e9, "")
    for (c, n), r in complex_.items():
        if c not in ("A", "B") or (only is not None and (c, n) not in only):
            continue
        for xyz in r["atoms"].values():
            d = math.dist(p, xyz)
            if d < best[0]:
                best = (d, f"{c}:{r['comp']}{n}")
    return best


def check_spec(path, complex_):
    path = Path(path)
    spec = json.loads(path.read_text(encoding="utf8"))
    assert len(spec) == 1, "one design ID per file keeps one run = one variant"
    (design_id, d), = spec.items()
    print(f"\n=== {path.name}  (design id {design_id}) ===")

    # 1. inputs
    files = [path.parent / f for f in d["input"]]
    check(all(f.is_file() for f in files), "referenced structure files exist", ", ".join(d["input"]))
    tables = [load_structure(f) for f in files]
    chains = [{c for c, _ in t} for t in tables]
    check(sum(len(c) for c in chains) == len(set().union(*chains)), "chain IDs distinct across files",
          " | ".join(",".join(sorted(c)) for c in chains))
    table = {}
    for t in tables:
        table.update(t)

    # 2. contig
    retained, designed = {}, []
    for tok in (t.strip() for t in d["contig"].split(",")):
        if tok == "/0":
            continue
        m = re.fullmatch(r"([A-Z])(\d+)-(\d+)", tok)
        if m:
            retained.setdefault(m.group(1), set()).update(range(int(m.group(2)), int(m.group(3)) + 1))
        elif re.fullmatch(r"\d+-\d+", tok):
            designed.append(tuple(int(x) for x in tok.split("-")))
        else:
            check(False, f"contig token {tok!r} parses")
    for c, got in sorted(retained.items()):
        have = {n for (cc, n) in table if cc == c}
        check(got == have, f"contig chain {c} = residues in file", f"{len(got)} vs {len(have)}")
    check(len(designed) == 3, "three designed CDR segments", str(designed))

    cond = d["conditions"]

    # 3. selectors resolve
    for key in ("C_CRD", "C_CRDGRP", "C_SEQ", "C_HOT", "C_CATH"):
        for s in selects(cond, key):
            n = len(resolve(s, table))
            check(n > 0, f"{key}: {s}", f"{n} atoms")

    # 3b/4. fixed-geometry selectors are backbone-only and complete
    for key in ("C_CRD", "C_CRDGRP"):
        names = {a[2] for s in selects(cond, key) for a in resolve(s, table)}
        check(names <= BACKBONE, f"{key} selects backbone atoms only", ",".join(sorted(names)))
    crd = {a[:2] for s in selects(cond, "C_CRD") for a in resolve(s, table)}
    her2 = {k for k in table if k[0] in ("C", "D")}
    check(crd == her2, "C_CRD covers all HER2 residues", f"{len(crd)}/{len(her2)}")
    grp = {a[:2] for s in selects(cond, "C_CRDGRP") for a in resolve(s, table)}
    fw = {k for k in table if k[0] == "H"}
    check(grp == fw, "C_CRDGRP covers every framework residue", f"{len(grp)}/{len(fw)}")

    # 5. sequence conditions
    seq = {a[:2] for s in selects(cond, "C_SEQ") for a in resolve(s, table)}
    check(her2 <= seq, "C_SEQ fixes HER2 sequence")
    fixed_fw = {k for k in seq if k[0] == "H"}
    freed = sorted(n for (c, n) in fw - fixed_fw)
    cys = sorted(n for (c, n) in fw if table[(c, n)]["comp"] == "CYS")
    check(all(("H", n) in fixed_fw for n in cys), "framework cysteines stay fixed", f"Cys at {cys}")
    bad = [n for n in freed if table[("H", n)]["comp"] in ("CYS", "GLY", "PRO")]
    check(not bad, "freed positions exclude Cys/Gly/Pro", f"freed {len(freed)}: {freed}" if freed else "none freed")

    # 6. hotspots
    hot = selects(cond, "C_HOT")
    pts = []
    for s in hot:
        hits = resolve(s, table)
        check(len(hits) == 1, f"hotspot {s} is exactly one atom")
        if len(hits) == 1:
            c, n, name = hits[0]
            dist, near = epitope_contacts(complex_, (c, n, name))
            same = complex_[(c, n)]["atoms"][name] == table[(c, n)]["atoms"][name]
            check(dist <= CONTACT_CUTOFF and same,
                  f"hotspot {table[(c, n)]['comp']}{n}/{name} touches trastuzumab in 1N8Z",
                  f"{dist:.2f} A to {near}; coordinates {'unchanged' if same else 'CHANGED'}")
            d3, near3 = epitope_contacts(complex_, (c, n, name), only=CDRH3)
            check(d3 <= CONTACT_CUTOFF,
                  f"hotspot {table[(c, n)]['comp']}{n}/{name} is within {CONTACT_CUTOFF} A of trastuzumab CDRH3",
                  f"{d3:.2f} A to {near3} (RFD4-Proteina's hotspot definition, with CDRH3 as the binder)")
            pts.append(table[(c, n)]["atoms"][name])

    # 7. compactness
    span = max((math.dist(p, q) for p in pts for q in pts), default=0.0)
    check(span <= MAX_HOTSPOT_SPAN, "hotspots fit one VHH paratope", f"max separation {span:.1f} A")

    # identity arithmetic (informational)
    n_fix = len(fixed_fw)
    lo, hi = sum(a for a, _ in designed), sum(b for _, b in designed)
    n_all = len(fw)
    print(f"  -- identity floor vs the source framework (zero chance matches in CDRs):")
    print(f"     fixed framework residues {n_fix} of {n_all}; designed CDR residues {lo}-{hi}")
    print(f"     total length {n_fix + len(freed) + lo}-{n_fix + len(freed) + hi}")
    for p in ((0.0,) if not freed else (0.0, 0.3, 0.5, 1.0)):
        a = (n_fix + p * len(freed)) / (n_all + lo)
        b = (n_fix + p * len(freed)) / (n_all + hi)
        tag = f" if freed positions recover the source residue {int(p * 100)}% of the time" if freed else ""
        print(f"     identity {b:.1%} (longest CDRs) .. {a:.1%} (shortest){tag}")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    complex_ = load_structure(HER2_CIF)
    for p in sys.argv[1:]:
        check_spec(p, complex_)
    print()
    if failures:
        print(f"{len(failures)} check(s) FAILED:")
        for f in failures:
            print("  -", f)
        sys.exit(1)
    print("all checks passed")


if __name__ == "__main__":
    main()
