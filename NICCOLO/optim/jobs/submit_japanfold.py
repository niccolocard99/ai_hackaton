#!/usr/bin/env python3
"""Submit a prepared JapanFold payload. DRY-RUN unless --go is given.

    python3 submit_japanfold.py rf3_payload.json            # show the plan and the cost estimate; spends nothing
    python3 submit_japanfold.py rf3_payload.json --go       # POST, wait for the job, download results to results/<payload name>/

Safe to re-run: the Idempotency-Key is derived from the payload hash, so a retry returns the first job instead of billing twice.
The API key comes from $JAPANFOLD_API_KEY (or ~/.config/japanfold/api_key) and is never printed. The key is also stripped
from any redirect to a different host (the result archive is served from object storage).
JapanFold holds two waits per account: run one payload at a time.
"""
import argparse, hashlib, json, os, pathlib, sys, time, urllib.error, urllib.parse, urllib.request

BASE = os.environ.get("JAPANFOLD_BASE_URL", "https://api.japanfold.aiand.com")
# USD per target, from this account's own billing history (273-target calibration jobs at $0.26/processor-hour). History, not a quote.
PER_TARGET = {"rf3": 0.0034, "boltz2": 0.0018}


def api_key():
    k = os.environ.get("JAPANFOLD_API_KEY")
    p = pathlib.Path("~/.config/japanfold/api_key").expanduser()
    if not k and p.exists(): k = p.read_text().strip()
    if not k: sys.exit("no JapanFold key: export JAPANFOLD_API_KEY in your shell (never paste it into a chat)")
    return k


class DropAuthOnForeignRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        new = super().redirect_request(req, fp, code, msg, headers, newurl)
        if new is not None and urllib.parse.urlparse(newurl).netloc != urllib.parse.urlparse(req.full_url).netloc:
            new.headers.pop("Authorization", None); new.unredirected_hdrs.pop("Authorization", None)
        return new


OPENER = urllib.request.build_opener(DropAuthOnForeignRedirect)


def call(method, path, body=None, extra=None, raw=False, timeout=130):
    h = {"Authorization": f"Bearer {api_key()}", "X-JapanFold-Client": "skill", "Content-Type": "application/json"}
    h.update(extra or {})
    req = urllib.request.Request(BASE + path, data=None if body is None else json.dumps(body).encode(), headers=h, method=method)
    try:
        with OPENER.open(req, timeout=timeout) as r:
            data = r.read()
            return data if raw else json.loads(data)
    except urllib.error.HTTPError as e:
        return {"_http": e.code, "_body": e.read().decode()[:800]}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("payload"); ap.add_argument("--go", action="store_true")
    a = ap.parse_args()
    payload = json.load(open(a.payload)); model, n = payload["model"], len(payload["targets"])
    idem = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()[:32]
    est = PER_TARGET.get(model, 0.01) * n
    print(f"payload {a.payload}: model={model} targets={n} params={payload.get('params')} idempotency-key={idem}")
    print(f"estimated cost ~${est:.2f} (history, not a quote)")
    if not a.go:
        print("DRY RUN: nothing submitted. Re-run with --go to submit."); return
    t_submit = time.time()
    r = call("POST", "/v1/predictions", payload, {"Idempotency-Key": idem})
    if "id" not in r: sys.exit(f"submit failed: {r}")
    jid = r["id"]; print(time.strftime("%H:%M:%S"), "SUBMITTED job", jid, r.get("status"), "warnings:", r.get("warnings"), flush=True)
    pathlib.Path(a.payload).with_suffix(".job_id").write_text(f"{jid}\n{t_submit:.0f}\n")
    while True:
        j = call("GET", f"/v1/jobs/{jid}?wait=60")
        st = j.get("status"); print(time.strftime("%H:%M:%S"), f"+{time.time() - t_submit:.0f}s", st, j.get("queue") or "", flush=True)
        if st not in ("queued", "running"): break
        if j.get("wait_refused"): time.sleep(int(j.get("retry_after", 10)))
    t_done = time.time()
    out = pathlib.Path("results") / pathlib.Path(a.payload).stem; out.mkdir(parents=True, exist_ok=True)
    json.dump(j, open(out / "job.json", "w"), indent=1)
    if st != "succeeded": sys.exit(f"job ended {st}: {j.get('error')}")
    if j.get("partial"): print("PARTIAL:", j["partial"])
    json.dump(call("GET", f"/v1/jobs/{jid}/results"), open(out / "results.json", "w"), indent=1)
    (out / "archive.zip").write_bytes(call("GET", f"/v1/jobs/{jid}/archive", raw=True, timeout=900))
    print(f"saved {out}\nTIMING: submit->job finished {t_done - t_submit:.0f}s, download {time.time() - t_done:.0f}s, total {time.time() - t_submit:.0f}s", flush=True)


if __name__ == "__main__":
    main()
