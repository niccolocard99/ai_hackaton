#!/usr/bin/env python3
"""Fetch full text for cached AMASS records that have it.

Usage: AMASS_API_KEY=amass_... python fetch_fulltext.py [--limit N] [--dry-run]

Cost: $0.01 per get-by-id. Only records with hasFulltext=True are fetched for
biomedcore; patents are fetched with claims+description. Everything is cached
under fulltext/<amassId>.json and never re-fetched, so re-runs are free.

Rate limit is 60 requests / 60 s, so we sleep 1.1 s between calls and back off
on HTTP 429 using the Retry-After header.
"""
import argparse
import glob
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

# Paper titles contain Greek letters and dashes that a Windows cp1252 console
# cannot encode; without this, printing a title crashes the run mid-way.
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "cache")
OUT = os.path.join(HERE, "fulltext")
BASE = "https://api.amass.tech/api/v1"

KEY = os.environ.get("AMASS_API_KEY")
if not KEY:
    sys.exit("set AMASS_API_KEY in the environment (never commit it)")

INCLUDES = {
    "biomedcore": ["fulltext"],
    "patentcore": ["claims", "description"],
}


def targets():
    """Yield (core, amassId, label) for every record worth fetching, deduped."""
    seen = set()
    for path in sorted(glob.glob(os.path.join(CACHE, "*.json"))):
        core = os.path.basename(path).partition("__")[0]
        if core not in INCLUDES:
            continue
        try:
            data = json.load(open(path, encoding="utf-8"))["data"]
        except Exception as exc:  # noqa: BLE001
            print(f"skip {path}: {exc}")
            continue
        for rec in data:
            aid = rec.get("amassId")
            if not aid or aid in seen:
                continue
            if core == "biomedcore" and not rec.get("hasFulltext"):
                continue
            seen.add(aid)
            label = (rec.get("title") or "")[:70]
            yield core, aid, label


def fetch(core, aid):
    qs = urllib.parse.urlencode([("include", i) for i in INCLUDES[core]])
    url = f"{BASE}/cores/{core}/records/{aid}?{qs}"
    req = urllib.request.Request(url, headers={"Authorization": f"Bearer {KEY}"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            if exc.code == 429:
                wait = int(exc.headers.get("Retry-After") or 2 ** (attempt + 2))
                print(f"  429, waiting {wait}s")
                time.sleep(wait)
                continue
            return {"_error": f"HTTP {exc.code}: {exc.read()[:200].decode('utf-8','replace')}"}
        except Exception as exc:  # noqa: BLE001
            if attempt == 3:
                return {"_error": str(exc)}
            time.sleep(2 ** attempt)
    return {"_error": "retries exhausted"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="stop after N fetches")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    todo = [t for t in targets() if not os.path.exists(os.path.join(OUT, f"{t[1]}.json"))]
    have = len(glob.glob(os.path.join(OUT, "*.json")))
    print(f"cached already: {have}   to fetch: {len(todo)}   est cost: ${len(todo)*0.01:.2f}")
    if args.dry_run:
        return

    done = failed = 0
    for core, aid, label in todo:
        if args.limit and done >= args.limit:
            break
        res = fetch(core, aid)
        if "_error" in res:
            failed += 1
            print(f"FAIL {aid} {label}: {res['_error']}")
        else:
            json.dump(res, open(os.path.join(OUT, f"{aid}.json"), "w", encoding="utf-8"))
            done += 1
            body = res.get("data", {})
            n = len(body.get("fulltext") or body.get("description") or "")
            print(f"ok   {aid} [{n:>6} chars] {label}")
        time.sleep(1.1)

    print(f"\nfetched={done} failed={failed} cost=~${done*0.01:.2f}")


if __name__ == "__main__":
    main()
