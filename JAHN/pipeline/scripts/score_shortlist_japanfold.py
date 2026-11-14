"""Score the loop's shortlist against HER2 with JapanFold `af2ig` — provenance, not a screen.

AF2-IG re-predicts the binder-target complex from an initial guess: PDB 1N8Z (the trastuzumab
Fab-HER2 crystal) supplies the HER2 chain and the VH backbone, and each design's VH sequence is the
binder. Every design shares the parent framework, so the binder chain length always matches the
backbone and the only difference between targets is the CDRH3 window.

Read this score the way NICCOLO/eda/AGENT_GUIDE.md section 6 says to read every fold metric: a
plausibility and provenance check on what we actually submit — not a screen, not a ranker. The VL
chain is absent (af2ig takes one binder chain), so absolute ipTM values are not comparable to
full-Fab predictions; the comparison between designs on the same scaffold is the useful signal.

    JAPANFOLD_API_KEY=... python scripts/score_shortlist_japanfold.py SHORTLIST.csv -k 5 -o runs/japanfold
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.request
from pathlib import Path

import polars as pl

BASE = "https://api.japanfold.aiand.com"
PDB_1N8Z = "https://files.rcsb.org/download/1N8Z.pdb"
VH_LEN = 120          # residues of 1N8Z chain B kept as the binder backbone (the Fv's VH)


def _req(url: str, data=None, token: str | None = None, timeout: int = 120):
    headers = {"Content-Type": "application/json"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    r = urllib.request.Request(url, data=(json.dumps(data).encode() if data else None), headers=headers,
                                method="POST" if data else "GET")
    with urllib.request.urlopen(r, timeout=timeout) as resp:
        return json.loads(resp.read().decode())


def structure_lines() -> tuple[str, str]:
    """The HER2 chain (C) and the VH backbone (chain B, first VH_LEN residues) as PDB ATOM lines."""
    cache = Path(os.environ.get("JF_CACHE", Path.home() / ".cache/japanfold"))
    cache.mkdir(parents=True, exist_ok=True)
    pdb = cache / "1n8z.pdb"
    if not pdb.exists():
        urllib.request.urlretrieve(PDB_1N8Z, pdb)
    her2, vh, seen = [], [], set()
    for ln in pdb.read_text().splitlines():
        if not ln.startswith("ATOM"):
            continue
        ch, resnum, atom = ln[21], int(ln[22:26]), ln[12:16].strip()
        if atom != "CA" or (ch, resnum) in seen:
            continue
        seen.add((ch, resnum))
        if ch == "C":
            her2.append(ln)
        elif ch == "B" and len(vh) < VH_LEN:
            vh.append(ln)
    return "\n".join(her2), "\n".join(vh)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("shortlist")
    ap.add_argument("-k", "--top", type=int, default=5)
    ap.add_argument("-o", "--out", default="runs/japanfold")
    ap.add_argument("--rank-col", default="rank_score")
    a = ap.parse_args()

    token = os.environ.get("JAPANFOLD_API_KEY")
    if not token:
        sys.exit("JAPANFOLD_API_KEY is not set (never write the key into a file)")
    df = pl.read_csv(a.shortlist)
    rank_col = a.rank_col if a.rank_col in df.columns else ("rank" if "rank" in df.columns else None)
    if rank_col:
        df = df.sort(rank_col)
    top = df.head(a.top)
    if "heavy_sequence" not in df.columns:
        sys.exit("the shortlist needs a heavy_sequence column")
    her2, vh_backbone = structure_lines()
    print(f"initial guess: {len(her2.splitlines())} HER2 + {len(vh_backbone.splitlines())} VH backbone residues")

    targets = []
    for r in top.iter_rows(named=True):
        vh_seq = r["heavy_sequence"]
        if len(vh_seq) != VH_LEN:
            print(f"skip {r['sequence_id']}: heavy is {len(vh_seq)} aa, backbone is {VH_LEN}")
            continue
        content = (f"target:\n  structure: |\n" +
                   "".join(f"    {ln}\n" for ln in (her2 + "\n" + vh_backbone).splitlines()) +
                   f"  chain: C\nbinder:\n  sequence: {vh_seq}\n  chain: B\n")
        targets.append({"name": r["sequence_id"], "content": content})
    if not targets:
        sys.exit("no scorable designs")

    job = _req(f"{BASE}/v1/predictions", {"model": "af2ig", "name": "loop-shortlist", "targets": targets}, token)
    jid = job.get("job") or job.get("id")
    print(f"job {jid}: {len(targets)} designs", flush=True)
    for _ in range(60):
        s = _req(f"{BASE}/v1/jobs/{jid}?wait=30", token=token)
        if s.get("status") in ("done", "failed", "succeeded"):
            break
        print(f"  status {s.get('status')}", flush=True)
    if s.get("status") not in ("done", "succeeded"):
        sys.exit(f"job ended as {s.get('status')}: {json.dumps(s)[:400]}")
    results = _req(f"{BASE}/v1/jobs/{jid}/results", token=token)
    rows = results.get("results") or results.get("rows") or results
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "af2ig_raw.json").write_text(json.dumps(results, indent=2)[:200000])
    keep = pl.DataFrame([{**{k: r.get(k) for k in ("name", "iptm", "plddt", "ptm", "interface_pae", "pae")},
                          "sequence_id": r.get("name")} for r in rows])
    keep.write_csv(out / "af2ig_scores.csv")
    print(keep)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
