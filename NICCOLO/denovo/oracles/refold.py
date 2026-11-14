#!/usr/bin/env python3
"""
Refolding oracle harness: score antibody designs against HER2 domain IV.

Design notes
------------
* submit / poll / collect are SEPARATE functions. `SUBMITTER` is a module-level
  hook; to route through NICCOLO/denovo/broker.py later, set
  `refold.SUBMITTER = broker_submit` (same signature) and nothing else changes.
* Every submission carries an `Idempotency-Key`, which is DERIVED from the batch
  content (sha256), so a re-run of an interrupted submit returns the SAME job
  instead of paying twice.
* A durable jsonl ledger (`runs.jsonl`) records every job and every per-design
  row. Anything already in the ledger is never resubmitted.
* Only ONE job is ever waited on per call (JapanFold allows 2 concurrent
  wait-holds per account).
* A failure of one design/target degrades that row only; the run continues.
* `pair_chains_iptm` is keyed by chain POSITION ("0","1",...). The engine orders
  chains by SORTED chain id, NOT input order. We therefore (a) build the mapping
  as sorted(chain_ids), and (b) VERIFY it against the returned structure's real
  per-chain residue counts, which are a unique fingerprint because every chain in
  our panel has a distinct length. See `verify_chain_indexing`.

Usage
-----
  python3 refold.py submit   --csv controls.csv --model rf3
  python3 refold.py run      --csv controls.csv --model rf3     # submit+poll+collect
  python3 refold.py collect  --job <id>
  python3 refold.py table
"""
import argparse, csv, hashlib, json, os, sys, time, urllib.error, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.environ.get("JAPANFOLD_BASE_URL", "https://api.japanfold.aiand.com")
LEDGER = os.path.join(HERE, "runs.jsonl")
STRUCT_DIR = os.path.join(HERE, "structures")

# ---------------------------------------------------------------- antigen
HER2_C = "CHQLCARGHCWGPGPTQCVNCSQFLRGQECVEECRVLQGLPREYVNARHCLPCHPECQPQNGSVTCFGPEADQCVACAHYKDPPFCVARCPS"  # 489-580
HER2_GAP = "GVKPDLSYMP"   # 581-590, real sequence recovered from PDB 1S78
HER2_D = "IWKFPDEEGACQPCPIN"  # 591-607
HER2_DIV_SINGLE = HER2_C + HER2_GAP + HER2_D   # 119 aa, one contiguous chain

# chain ids. Chosen so INPUT order (H,L,A,...) differs from SORTED order
# (A,...,H,L): that difference is what makes the indexing check meaningful.
CH_HEAVY, CH_LIGHT, CH_AG1, CH_AG2 = "H", "L", "A", "B"
AG_ROLE = "antigen"


def _key(*parts):
    return hashlib.sha256("|".join(parts).encode()).hexdigest()


# ---------------------------------------------------------------- HTTP
def _req(method, path, body=None, timeout=120):
    url = path if path.startswith("http") else BASE + path
    data = json.dumps(body).encode() if body is not None else None
    r = urllib.request.Request(url, data=data, method=method)
    r.add_header("Authorization", "Bearer " + os.environ["JAPANFOLD_API_KEY"])
    r.add_header("X-JapanFold-Client", "niccolo-refold")
    if data:
        r.add_header("Content-Type", "application/json")
    if isinstance(body, dict) and body.pop("_idem", None):
        pass
    try:
        with urllib.request.urlopen(r, timeout=timeout) as fh:
            return json.loads(fh.read().decode())
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"{method} {url} -> {e.code}: {e.read().decode()[:600]}") from None


def _post_idem(path, body, idem, timeout=120):
    url = BASE + path
    r = urllib.request.Request(url, data=json.dumps(body).encode(), method="POST")
    r.add_header("Authorization", "Bearer " + os.environ["JAPANFOLD_API_KEY"])
    r.add_header("X-JapanFold-Client", "niccolo-refold")
    r.add_header("Content-Type", "application/json")
    r.add_header("Idempotency-Key", idem)
    try:
        with urllib.request.urlopen(r, timeout=timeout) as fh:
            return json.loads(fh.read().decode())
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"POST {url} -> {e.code}: {e.read().decode()[:600]}") from None


# ---------------------------------------------------------------- ledger
def ledger_append(rec):
    rec["_at"] = time.time()
    with open(LEDGER, "a") as fh:
        fh.write(json.dumps(rec) + "\n")


def ledger_read():
    if not os.path.exists(LEDGER):
        return []
    out = []
    for line in open(LEDGER):
        line = line.strip()
        if line:
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                pass
    return out


# ---------------------------------------------------------------- targets
def build_target(row):
    """-> (target_name, yaml_content, chain_spec{chain_id: (role, length)})"""
    sid = row["sequence_id"]
    heavy = (row.get("heavy_sequence") or "").strip()
    light = (row.get("light_sequence") or "").strip()
    mode = (row.get("antigen_mode") or "single").strip()
    if not heavy:
        raise ValueError(f"{sid}: heavy_sequence is required")

    chains = []                      # (chain_id, sequence, role) in INPUT order
    chains.append((CH_HEAVY, heavy, "heavy"))
    if light:
        chains.append((CH_LIGHT, light, "light"))
    if mode == "split":
        chains.append((CH_AG1, HER2_C, AG_ROLE))
        chains.append((CH_AG2, HER2_D, AG_ROLE))
    else:
        chains.append((CH_AG1, HER2_DIV_SINGLE, AG_ROLE))

    lines = ["sequences:"]
    for cid, seq, _role in chains:
        lines.append(f"  - protein: {{id: {cid}, sequence: {seq}}}")
    spec = {cid: (role, len(seq)) for cid, seq, role in chains}
    return sid, "\n".join(lines) + "\n", spec


# ---------------------------------------------------------------- submit (swappable)
def japanfold_submit(targets, model, name, params=None):
    """targets: [(name, yaml)] -> job_id. Replaceable by a broker submitter."""
    body = {
        "model": model,
        "name": name,
        "targets": [{"name": n, "content": c} for n, c in targets],
    }
    if params:
        body["params"] = params
    idem = _key("submit", model, name, *(n + c for n, c in targets),
                json.dumps(params or {}, sort_keys=True))
    j = _post_idem("/v1/predictions", body, idem)
    return j["id"], idem


SUBMITTER = japanfold_submit     # <-- set to broker_submit later


def submit(targets, model, name, params=None):
    job_id, idem = SUBMITTER(targets, model, name, params)
    ledger_append({"rec": "job", "job": job_id, "model": model, "name": name,
                   "idempotency_key": idem, "targets": [n for n, _ in targets]})
    return job_id


# ---------------------------------------------------------------- poll (ONE job per call)
def poll(job_id, wait=60):
    return _req("GET", f"/v1/jobs/{job_id}?wait={int(wait)}", timeout=int(wait) + 40)


def wait_for(job_id, wait=60, max_s=1800, verbose=True):
    t0 = time.time()
    while time.time() - t0 < max_s:
        j = poll(job_id, wait)
        st = j.get("status")
        if verbose:
            q = j.get("queue") or {}
            print(f"  [{int(time.time()-t0):>4}s] {job_id[:12]} {st}"
                  + (f" queue={q.get('reason')} ahead={q.get('accounts_ahead')}" if q else ""),
                  flush=True)
        if st not in ("queued", "running"):
            return j
    return {"status": "timeout", "id": job_id}


# ---------------------------------------------------------------- chain indexing
def chain_order(chain_spec):
    """The engine's chain POSITION order. Documented/observed: SORTED chain id."""
    return sorted(chain_spec.keys())


def struct_chain_lengths(path):
    """Real per-chain residue count + mean pLDDT from a returned mmCIF/PDB.
    This is the ground truth we check the position mapping against."""
    per = {}
    if path.endswith((".cif", ".mmcif")):
        hdr, inat = [], False
        for line in open(path):
            if line.startswith("_atom_site."):
                hdr.append(line.strip().split(".", 1)[1]); inat = True; continue
            if inat:
                if not line.startswith(("ATOM", "HETATM")):
                    if per: break
                    continue
                f = line.split()
                if len(f) < len(hdr): continue
                d = dict(zip(hdr, f))
                cid = d.get("auth_asym_id") or d.get("label_asym_id")
                rid = d.get("auth_seq_id") or d.get("label_seq_id")
                try: b = float(d.get("B_iso_or_equiv", "nan"))
                except ValueError: b = float("nan")
                per.setdefault(cid, {})[rid] = b
    else:
        for line in open(path):
            if line.startswith("ATOM"):
                cid, rid = line[21], line[22:27].strip()
                try: b = float(line[60:66])
                except ValueError: b = float("nan")
                per.setdefault(cid, {})[rid] = b
    out = {}
    for cid, res in per.items():
        vals = [v for v in res.values() if v == v]
        out[cid] = {"n_res": len(res),
                    "mean_plddt": round(sum(vals) / len(vals), 2) if vals else None}
    return out


def verify_chain_indexing(chain_spec, observed):
    """Compare submitted chain lengths to the structure's real per-chain lengths.
    Every chain in the control panel has a UNIQUE length, so matching lengths
    pins chain identity without trusting any ordering convention."""
    checks, ok = [], True
    for cid, (role, n) in sorted(chain_spec.items()):
        got = observed.get(cid)
        match = bool(got) and got["n_res"] == n
        ok &= match
        checks.append({"chain": cid, "role": role, "submitted_len": n,
                       "structure_len": got["n_res"] if got else None, "match": match})
    return ok, checks


def pair_iptm_lookup(pair_matrix, order, a, b):
    """pair_matrix is keyed by chain POSITION as a string."""
    if not pair_matrix:
        return None
    try:
        i, j = order.index(a), order.index(b)
    except ValueError:
        return None
    for ii, jj in ((i, j), (j, i)):
        row = pair_matrix.get(str(ii))
        if isinstance(row, dict) and str(jj) in row:
            return row[str(jj)]
        if isinstance(row, list) and jj < len(row):
            return row[jj]
    return None


# ---------------------------------------------------------------- collect
def fetch_structure(job_id, artifacts, target_name):
    """Artifacts carry `path` + `target` (not `name`), and `url` is BASE-relative."""
    os.makedirs(STRUCT_DIR, exist_ok=True)
    cands = [a for a in artifacts if a.get("type") in (None, "structure", "mmcif", "pdb", "cif")]
    pick = ([a for a in cands if a.get("target") == target_name]
            or [a for a in cands if target_name in str(a.get("path") or a.get("name") or "")])
    if not pick:
        return None
    a = pick[0]
    fname = a.get("path") or a.get("name") or f"{target_name}.cif"
    dest = os.path.join(STRUCT_DIR, f"{job_id[:8]}_{os.path.basename(fname)}")
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    url = a["url"]
    r = urllib.request.Request(url if url.startswith("http") else BASE + url, method="GET")
    r.add_header("Authorization", "Bearer " + os.environ["JAPANFOLD_API_KEY"])
    try:
        with urllib.request.urlopen(r, timeout=180) as fh, open(dest, "wb") as out:
            out.write(fh.read())
        return dest
    except Exception as e:
        print(f"    ! structure fetch failed for {target_name}: {e}", flush=True)
        return None


def job_cost(job_id):
    for p in (f"/v1/jobs/{job_id}/usage", f"/v1/jobs/{job_id}"):
        try:
            d = _req("GET", p, timeout=60)
        except RuntimeError:
            continue
        for k in ("amount", "usd", "cost"):
            if isinstance(d.get(k), (int, float)):
                return float(d[k])
        for holder in ("usage", "billing", "detail"):
            h = d.get(holder) or {}
            if isinstance(h, dict):
                if isinstance(h.get("amount"), (int, float)):
                    return float(h["amount"])
                if isinstance(h.get("amount_micros"), (int, float)):
                    return h["amount_micros"] / 1e6
        if isinstance(d.get("amount_micros"), (int, float)):
            return d["amount_micros"] / 1e6
    return None


def collect(job_id, specs, model=None):
    """specs: {target_name: chain_spec}. Returns list of per-design dict rows."""
    job = _req("GET", f"/v1/jobs/{job_id}", timeout=60)
    res = _req("GET", f"/v1/jobs/{job_id}/results", timeout=120)
    cost = job_cost(job_id)
    arts = res.get("artifacts") or job.get("artifacts") or []
    rows_in = res.get("rows") or res.get("targets") or []
    out = []
    for r in rows_in:
        name = r.get("name") or r.get("target") or r.get("id")
        spec = specs.get(name, {})
        rec = {"rec": "design", "job": job_id, "model": model or job.get("model"),
               "sequence_id": name, "status": "ok", "cost_usd": cost,
               "iptm": r.get("iptm"), "ptm": r.get("ptm"),
               # rf3/of3 rows call the complex pLDDT `plddt` (0-1 scale)
               "complex_plddt": r.get("complex_plddt", r.get("plddt")),
               "confidence_score": r.get("confidence_score", r.get("ranking_score")),
               "has_clash": r.get("has_clash"), "n_tokens": r.get("n_tokens"),
               "runtime_s": r.get("runtime_s"), "row_status": r.get("status"),
               "chains_ptm": r.get("chains_ptm"),
               "pair_chains_iptm_raw": r.get("pair_chains_iptm")}
        try:
            order = chain_order(spec) if spec else []
            rec["position_to_chain"] = {str(i): c for i, c in enumerate(order)}
            pm = r.get("pair_chains_iptm")
            ab = [c for c, (role, _n) in spec.items() if role in ("heavy", "light")]
            ag = [c for c, (role, _n) in spec.items() if role == AG_ROLE]
            pairs = {}
            for a in ab:
                for g in ag:
                    pairs[f"{a}-{g}"] = pair_iptm_lookup(pm, order, a, g)
            rec["ab_ag_pairs"] = pairs
            vals = [v for v in pairs.values() if isinstance(v, (int, float))]
            rec["pair_iptm_ab_ag_max"] = max(vals) if vals else None
            rec["pair_iptm_ab_ag_mean"] = round(sum(vals) / len(vals), 4) if vals else None
            # internal (non-interface) pairs, used as the indexing anchors
            internal = {}
            if CH_HEAVY in spec and CH_LIGHT in spec:
                internal["H-L_obligate_Fv_pair"] = pair_iptm_lookup(pm, order, CH_HEAVY, CH_LIGHT)
            if len(ag) == 2:
                internal["A-B_same_HER2_domain"] = pair_iptm_lookup(pm, order, ag[0], ag[1])
            rec["internal_pairs"] = internal

            path = fetch_structure(job_id, arts, name)
            if path:
                obs = struct_chain_lengths(path)
                ok, checks = verify_chain_indexing(spec, obs)
                rec["chain_id_verified"] = ok
                rec["chain_checks"] = checks
                rec["per_chain_plddt"] = {c: v["mean_plddt"] for c, v in obs.items()}
                abp = [obs[c]["mean_plddt"] for c in ab
                       if c in obs and obs[c]["mean_plddt"] is not None]
                rec["ab_mean_plddt"] = round(sum(abp) / len(abp), 2) if abp else None
                rec["structure"] = path
        except Exception as e:
            rec["status"] = "partial"
            rec["error"] = f"{type(e).__name__}: {e}"
        ledger_append(rec)
        out.append(rec)
    missing = set(specs) - {r["sequence_id"] for r in out}
    for m in sorted(missing):
        rec = {"rec": "design", "job": job_id, "sequence_id": m,
               "status": "failed", "error": "no result row returned", "cost_usd": cost}
        ledger_append(rec)
        out.append(rec)
    return out


# ---------------------------------------------------------------- driver
def load_csv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def cmd_run(args):
    rows = load_csv(args.csv)
    done = {r["sequence_id"] for r in ledger_read()
            if r.get("rec") == "design" and r.get("status") == "ok"
            and r.get("model") == args.model}
    todo, specs = [], {}
    for row in rows:
        sid = row["sequence_id"]
        if sid in done and not args.force:
            print(f"  skip {sid} (already in ledger for {args.model})")
            continue
        try:
            name, yaml, spec = build_target(row)
        except Exception as e:
            print(f"  ! {sid}: {e}")
            ledger_append({"rec": "design", "sequence_id": sid, "status": "failed",
                           "error": str(e), "model": args.model})
            continue
        todo.append((name, yaml)); specs[name] = spec
    if not todo:
        print("nothing to do"); return
    print(f"submitting {len(todo)} targets to {args.model} as one job")
    job_id = submit(todo, args.model, args.name or f"her2-controls-{args.model}",
                    params={"seed": args.seed} if args.seed is not None else None)
    print("job:", job_id)
    ledger_append({"rec": "specs", "job": job_id, "specs":
                   {k: {c: list(v) for c, v in s.items()} for k, s in specs.items()}})
    j = wait_for(job_id, wait=args.wait, max_s=args.max_s)
    print("final status:", j.get("status"))
    if j.get("status") != "succeeded":
        print(json.dumps({k: j.get(k) for k in ("status", "error", "partial")}, indent=1))
        if j.get("status") == "timeout":
            return
    for r in collect(job_id, specs, model=args.model):
        print(f"  {r['sequence_id']:<22} {r['status']:<8} iptm={r.get('iptm')} "
              f"pair_max={r.get('pair_iptm_ab_ag_max')} verified={r.get('chain_id_verified')}")


def cmd_collect(args):
    specs = {}
    for r in ledger_read():
        if r.get("rec") == "specs" and r.get("job") == args.job:
            specs = {k: {c: tuple(v) for c, v in s.items()} for k, s in r["specs"].items()}
    for r in collect(args.job, specs, model=args.model):
        print(json.dumps(r, indent=1))


def cmd_table(args):
    recs = [r for r in ledger_read() if r.get("rec") == "design"]
    best = {}
    for r in recs:
        best[(r.get("model"), r.get("sequence_id"))] = r
    hdr = f"{'model':<8} {'sequence_id':<22} {'iptm':>6} {'ptm':>6} {'pair_max':>9} {'pair_mean':>9} {'ab_plddt':>9} {'cplddt':>7} {'ok':>4}"
    print(hdr); print("-" * len(hdr))
    for (m, sid), r in sorted(best.items(), key=lambda kv: (str(kv[0][0]), str(kv[0][1]))):
        f = lambda v, w=6: (f"{v:>{w}.3f}" if isinstance(v, (int, float)) else f"{'-':>{w}}")
        print(f"{str(m):<8} {str(sid):<22} {f(r.get('iptm'))} {f(r.get('ptm'))} "
              f"{f(r.get('pair_iptm_ab_ag_max'),9)} {f(r.get('pair_iptm_ab_ag_mean'),9)} "
              f"{f(r.get('ab_mean_plddt'),9)} {f(r.get('complex_plddt'),7)} "
              f"{str(r.get('chain_id_verified')):>4}")
    tot = {}
    for r in recs:
        if r.get("job") and isinstance(r.get("cost_usd"), (int, float)):
            tot[r["job"]] = r["cost_usd"]
    print(f"\ntotal spend across {len(tot)} jobs: ${sum(tot.values()):.4f}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run"); r.set_defaults(fn=cmd_run)
    r.add_argument("--csv", default=os.path.join(HERE, "controls.csv"))
    r.add_argument("--model", default="rf3")
    r.add_argument("--name"); r.add_argument("--seed", type=int)
    r.add_argument("--wait", type=int, default=90)
    r.add_argument("--max-s", dest="max_s", type=int, default=1800)
    r.add_argument("--force", action="store_true")
    c = sub.add_parser("collect"); c.set_defaults(fn=cmd_collect)
    c.add_argument("--job", required=True); c.add_argument("--model")
    t = sub.add_parser("table"); t.set_defaults(fn=cmd_table)
    a = ap.parse_args()
    if "JAPANFOLD_API_KEY" not in os.environ:
        sys.exit("JAPANFOLD_API_KEY is not set")
    a.fn(a)


if __name__ == "__main__":
    main()
