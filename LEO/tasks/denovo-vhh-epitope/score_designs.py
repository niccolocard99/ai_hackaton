#!/usr/bin/env python3
"""Join the per-design evidence into one scored table and pick a diverse shortlist.

Inputs (all keyed by sequence_id):
  designs.csv           from make_designs_table.py (heavy_sequence)
  designs_metrics.csv   from make_designs_table.py (RFD4 metrics)
  denovo report         from NICCOLO/scripts/check_denovo.py -o ... (track, identity, CDRH3, passes)
  epitope_check.csv     from epitope_check.py (geometry against 1N8Z)

    ./score_designs.py --designs D --metrics M --report R --epitope E -o scored.csv [--shortlist 12]

What the score is, and is not. These are filters for DESIGN INTENT: is the design valid, does it land on
trastuzumab's epitope, does it touch the hotspots, does it overlap trastuzumab's Fv. None of them has been
shown to predict SPR affinity, competition in the lab or thermostability: NICCOLO/eda/AGENT_GUIDE.md section 6
measured fold-confidence metrics at chance for binders, and there is no labelled de novo data. Every
threshold below is the author's choice, set before looking at the full batch, not derived from data.

Rules (in the order they are applied):
  eligible   base validation passes for the assigned track (de_novo), RFD4 hotspot_obedience == 1 and no
             clashing atoms
  on_target  >= 3 of 4 hotspot atoms within 4.5 A of the VHH, >= 60% of the VHH's HER2 contacts inside
             trastuzumab's epitope, and VHH/Fv steric overlap (> 0 atom pairs under 3 A)
  track_risk an ESTIMATE of how likely the official abguard is to call the design 70%+ identical to the
             caplacizumab framework it was built on. The local checker measures identity to trastuzumab
             only, so it cannot see this. fw_identity_est = (83 fixed framework residues + 2 expected chance
             matches in the CDRs) / length. HIGH >= 0.70, MEDIUM >= 0.67, else LOW. It is arithmetic, not
             abguard's number; fixedfw designs only (freefw framework identity is lower and not modelled)
<<<<<<< HEAD
  strong_epitope  epitope_cov >= 0.60 and contacts_in_ep >= 0.70: the team's definition (EXP-011 and
             NICCOLO/SCORED_DENOVO_VHH.md, which found epitope coverage the only axis with a usable signal for
             these VHHs; binding is unpredicted and thermostability has no baseline)
  rank       eligible, then strong_epitope, then on_target, then not-HIGH track_risk, then larger epitope
             coverage, then larger share of contacts inside the epitope, then more hotspots
  shortlist  walk the ranking, skip HIGH track_risk, and keep a design only if its CDRH3 is under 70%
             identical to every design already kept (a few distinct molecules beat several near-copies)
=======
  rank       eligible, then on_target, then lower track_risk, then more hotspots, then a larger share of
             contacts inside the epitope, then larger epitope coverage
  shortlist  walk the ranking and keep a design only if its CDRH3 is under 70% identical to every design
             already kept (a few distinct molecules beat several near-copies)
>>>>>>> origin/main
"""
import argparse
import csv
import sys

FIXED_FRAMEWORK = 83
CHANCE_CDR_MATCHES = 2
RISK_ORDER = {"LOW": 0, "MEDIUM": 1, "HIGH": 2, "n/a": 1}


def read(path, key="sequence_id"):
    with open(path, newline="", encoding="utf8") as f:
        return {r[key]: r for r in csv.DictReader(f)}


def identity(a, b):
    """1 - edit distance / longer length. Gaps cost as much as mismatches: free gaps would let the conserved
    CDRH3 flanks (AA...DY) inflate the similarity of unrelated loops."""
    n, m = len(a), len(b)
    if not n or not m:
        return 0.0
    prev = list(range(m + 1))
    for i in range(1, n + 1):
        cur = [i] + [0] * m
        for j in range(1, m + 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] != b[j - 1]))
        prev = cur
    return 1 - prev[m] / max(n, m)


def to_float(x, default=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--designs", required=True)
    ap.add_argument("--metrics", required=True)
    ap.add_argument("--report", required=True, nargs="+", help="check_denovo.py report(s)")
    ap.add_argument("--epitope", required=True, nargs="+", help="epitope_check.py output(s)")
    ap.add_argument("-o", "--out", default="scored.csv")
    ap.add_argument("--shortlist", type=int, default=12)
    a = ap.parse_args()

    designs, metrics = read(a.designs), read(a.metrics)
    report, epitope = {}, {}
    for p in a.report:
        report.update(read(p))
    for p in a.epitope:
        epitope.update(read(p))

    rows, missing = [], []
    for sid, d in designs.items():
        m, r, e = metrics.get(sid), report.get(sid), epitope.get(sid)
        if not (m and r and e):
            missing.append(sid)
            continue
        heavy = d["heavy_sequence"]
        est = (FIXED_FRAMEWORK + CHANCE_CDR_MATCHES) / len(heavy)
        variant = m["variant"]
        risk = "n/a" if variant != "fixedfw" else ("HIGH" if est >= 0.70 else "MEDIUM" if est >= 0.67 else "LOW")
        eligible = (r["passes"].lower() == "true" and r["track_assigned"] == "de_novo"
                    and to_float(m["hotspot_obedience"]) == 1.0 and to_float(m["clash_frac_1p5A"]) == 0.0)
        n_hot = int(e["n_hotspots_4p5"])
        in_ep = to_float(e["contacts_in_ep"])
        on_target = n_hot >= 3 and in_ep >= 0.60 and int(e["fv_pairs_3p0"]) > 0
<<<<<<< HEAD
        strong = to_float(e["epitope_cov"]) >= 0.60 and in_ep >= 0.70
=======
>>>>>>> origin/main
        rows.append({
            "sequence_id": sid, "variant": variant, "run": m["run"], "heavy_sequence": heavy, "length": len(heavy),
            "cdrh3": r["cdrh3"], "track_assigned": r["track_assigned"], "base_validation": "PASS" if r["passes"].lower() == "true" else "FAIL",
            "failures": r["failures"], "identity_to_trastuzumab": round(to_float(r["identity_to_trastuzumab"]), 3),
            "fw_identity_est": round(est, 3), "track_risk": risk,
            "hotspot_obedience": m["hotspot_obedience"], "clash_frac": m["clash_frac_1p5A"],
            "hbond_pairs": m["hbond_pairs"], "saltbridge_pairs": m["saltbridge_pairs"], "hydrophobic_pairs": m["hydrophobic_pairs"],
            "n_hotspots_4p5": n_hot, "n_contacts": e["n_contacts"], "epitope_cov": e["epitope_cov"], "contacts_in_ep": e["contacts_in_ep"],
<<<<<<< HEAD
            "fv_pairs_3p0": e["fv_pairs_3p0"], "eligible": eligible, "on_target": on_target, "strong_epitope": strong})

    rows.sort(key=lambda x: (not x["eligible"], not x["strong_epitope"], not x["on_target"], x["track_risk"] == "HIGH",
                             -to_float(x["epitope_cov"]), -to_float(x["contacts_in_ep"]), -x["n_hotspots_4p5"]))
=======
            "fv_pairs_3p0": e["fv_pairs_3p0"], "eligible": eligible, "on_target": on_target})

    rows.sort(key=lambda x: (not x["eligible"], not x["on_target"], RISK_ORDER[x["track_risk"]], -x["n_hotspots_4p5"],
                             -to_float(x["contacts_in_ep"]), -to_float(x["epitope_cov"])))
>>>>>>> origin/main
    for i, x in enumerate(rows, 1):
        x["rank"] = i

    keep = []
    for x in rows:
<<<<<<< HEAD
        if not (x["eligible"] and x["on_target"]) or x["track_risk"] == "HIGH":
=======
        if not (x["eligible"] and x["on_target"]):
>>>>>>> origin/main
            continue
        if all(identity(x["cdrh3"], k["cdrh3"]) < 0.70 for k in keep):
            keep.append(x)
        if len(keep) == a.shortlist:
            break
    short = {k["sequence_id"] for k in keep}
    for x in rows:
        x["shortlist"] = x["sequence_id"] in short

    cols = ["rank", "shortlist", "sequence_id", "variant", "length", "cdrh3", "base_validation", "track_assigned",
<<<<<<< HEAD
            "identity_to_trastuzumab", "fw_identity_est", "track_risk", "eligible", "strong_epitope", "on_target", "n_hotspots_4p5",
=======
            "identity_to_trastuzumab", "fw_identity_est", "track_risk", "eligible", "on_target", "n_hotspots_4p5",
>>>>>>> origin/main
            "n_contacts", "epitope_cov", "contacts_in_ep", "fv_pairs_3p0", "hbond_pairs", "saltbridge_pairs",
            "hydrophobic_pairs", "hotspot_obedience", "clash_frac", "failures", "run", "heavy_sequence"]
    with open(a.out, "w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    from collections import Counter
    n = len(rows)
    print(f"{n} designs scored ({len(missing)} missing evidence)  ->  {a.out}")
    print(f"base validation PASS (de novo): {sum(x['base_validation'] == 'PASS' for x in rows)}   "
          f"eligible: {sum(x['eligible'] for x in rows)}   on_target: {sum(x['on_target'] for x in rows)}   "
<<<<<<< HEAD
          f"eligible AND on_target: {sum(x['eligible'] and x['on_target'] for x in rows)}   "
          f"eligible AND strong_epitope: {sum(x['eligible'] and x['strong_epitope'] for x in rows)}")
=======
          f"eligible AND on_target: {sum(x['eligible'] and x['on_target'] for x in rows)}")
>>>>>>> origin/main
    print("track_risk among eligible+on_target:", dict(Counter(x["track_risk"] for x in rows if x["eligible"] and x["on_target"])))
    fails = Counter(x["failures"][:90] for x in rows if x["base_validation"] == "FAIL")
    if fails:
        print("base-validation failures:", dict(fails))
    print(f"\nshortlist ({len(keep)}, CDRH3 < 70% identical to every earlier pick):")
    print(f"{'#':>3} {'sequence_id':<44} {'len':>3} {'CDRH3':<22} {'hot':>3} {'inEp':>5} {'cov':>5} {'risk':>6} {'id_trz':>6}")
    for x in keep:
        print(f"{x['rank']:>3} {x['sequence_id'][-44:]:<44} {x['length']:>3} {x['cdrh3']:<22} {x['n_hotspots_4p5']:>3} "
              f"{x['contacts_in_ep']:>5} {x['epitope_cov']:>5} {x['track_risk']:>6} {x['identity_to_trastuzumab']:>6}")
    if missing:
        print("\nmissing evidence for:", ", ".join(missing[:8]), "..." if len(missing) > 8 else "", file=sys.stderr)


if __name__ == "__main__":
    main()
