#!/usr/bin/env python3
"""GPU stages on the local H100: NanoMelt (predicted Tm) and Boltz-2 (structure plausibility).

What these are for, stated plainly, because the ledger measured it:

  * Boltz-2 confidence does NOT rank CDRH3 variants by affinity. On 60 measured SPR
    variants every metric sits at chance (ipTM AUROC 0.502, min_ipae 0.497, EXP-004),
    on 180 binders it is uncorrelated with KD (Spearman -0.065, CI spans zero), and
    ipTM >= 0.7 passes 100% of true non-binders. Nothing here ranks designs.
  * What a fold can say that the sequence models cannot is whether the CDRH3 stays in
    the trastuzumab pose. We report the H3 displacement from a trastuzumab control
    folded in the same batch, after superposing on HER2, with a noise floor measured
    from the control at a second seed. The H3 loop is only predicted to about 2 A
    (TH-011), so this flags outliers and nothing finer.
  * NanoMelt is nanobody-trained and out of domain for this VH (EXP-005). It found
    CDRH3 thermostability-inert, so a CDRH3-only design should land within its 4.1 C
    error of the control. A design that does not deserves a second look; that is all.

The H100 is shared. Before each GPU job we wait for the card to be free, and Boltz
runs are cut into small chunks so another job can slot in between them.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from her2lib import LOOP_SLICE, TZ_LOOP, TZ_VH, TZ_VL, vh_with_loop  # noqa: E402

GPU = HERE / "gpu"
NANOMELT = "/home/hackathon/nick-work/nmenv/bin/nanomelt"
RFD = Path("/home/hackathon/hackathon-team-3/apheris/rfd4")
REF_Q = RFD / "runs/nick-ref/prepared/holo/queries.json"   # carries HER2 chains C and D (1N8Z numbering)
REF_T = RFD / "runs/nick-ref/prepared/holo/template.cif"
BOLTZ_IMG = "quay.io/apheris/oslo-hackathon:apheris-boltz-v1.3.0-boltz2"
LOOP_START = LOOP_SLICE.start          # 0-based index of the 10-mer in the VH


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# --------------------------------------------------------------------------- GPU etiquette
def gpu_used_mb() -> float:
    out = subprocess.run(["nvidia-smi", "--query-gpu=memory.used", "--format=csv,noheader,nounits"],
                         capture_output=True, text=True).stdout.strip().splitlines()
    return max(float(x) for x in out) if out else 0.0


def wait_for_gpu(max_wait_s: int = 1500, poll_s: int = 15, free_below_mb: float = 3000) -> bool:
    """Block until nobody else is using the card. False if we gave up waiting."""
    t0, said = time.time(), False
    while gpu_used_mb() > free_below_mb:
        if time.time() - t0 > max_wait_s:
            log(f"GPU still busy after {max_wait_s}s ({gpu_used_mb():.0f} MiB used); skipping this stage")
            return False
        if not said:
            log(f"GPU busy ({gpu_used_mb():.0f} MiB used by another job); waiting for it to free up")
            said = True
        time.sleep(poll_s)
    return True


# --------------------------------------------------------------------------- NanoMelt
def run_nanomelt(loops: list[str], tag: str, timeout: int = 1800) -> pd.DataFrame:
    """Predicted Tm of each design VH, and the change against the trastuzumab VH."""
    d = GPU / tag
    d.mkdir(parents=True, exist_ok=True)
    fa, out = d / "in.fa", d / "nanomelt.csv"
    ids = {f"d{i:05d}": l for i, l in enumerate(loops)}
    with fa.open("w") as fh:
        fh.write(f">WT\n{TZ_VH}\n")
        for k, l in ids.items():
            fh.write(f">{k}\n{vh_with_loop(l)}\n")
    t0 = time.time()
    env = dict(os.environ, PATH=f"{Path(NANOMELT).parent}:{os.environ.get('PATH', '')}")
    p = subprocess.run([NANOMELT, "predict", "-i", str(fa), "-o", str(out), "-align", "-ncpu", "4"],
                       capture_output=True, text=True, timeout=timeout, cwd=d, env=env)
    if not out.exists():
        raise RuntimeError(f"nanomelt produced no output:\n{p.stdout[-800:]}\n{p.stderr[-800:]}")
    df = pd.read_csv(out)
    tm = dict(zip(df["ID"].astype(str), df["NanoMelt Tm (C)"].astype(float)))
    wt = tm["WT"]
    res = pd.DataFrame({"cdrh3": list(ids.values()), "tm_pred": [tm[k] for k in ids]})
    res["dtm_vs_wt"] = res.tm_pred - wt
    res.attrs["wt_tm"] = wt
    log(f"NanoMelt: {len(res)} designs + control in {time.time() - t0:.0f}s; control Tm {wt:.2f} C")
    return res


# --------------------------------------------------------------------------- structure helpers
def read_ca(path: Path) -> dict[str, np.ndarray]:
    """CA coordinates per chain from a Boltz mmCIF: columns are (label_seq_id, x, y, z)."""
    chains: dict[str, list] = {}
    with path.open() as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            f = line.split()
            if f[3] != "CA":
                continue
            chains.setdefault(f[6], []).append((int(f[8]), float(f[10]), float(f[11]), float(f[12])))
    return {c: np.array(v, float) for c, v in chains.items()}


def kabsch(mobile: np.ndarray, ref: np.ndarray):
    """R, t with  mobile @ R.T + t  ~  ref  (least-squares superposition)."""
    mc, rc = mobile.mean(0), ref.mean(0)
    U, _, Vt = np.linalg.svd((mobile - mc).T @ (ref - rc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, rc - mc @ R.T


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(((a - b) ** 2).sum(1).mean()))


def pose_metrics(model: Path, wt_model: Path) -> dict:
    """H3 displacement from the control, in the HER2 frame (pose) and after local fit (shape)."""
    m, w = read_ca(model), read_ca(wt_model)
    her_m = np.vstack([m["C"][:, 1:], m["D"][:, 1:]])
    her_w = np.vstack([w["C"][:, 1:], w["D"][:, 1:]])
    R, t = kabsch(her_m, her_w)
    sl = slice(LOOP_START, LOOP_START + len(TZ_LOOP))
    vh_m, vh_w = m["B"][:, 1:] @ R.T + t, w["B"][:, 1:]
    vl_m, vl_w = m["A"][:, 1:] @ R.T + t, w["A"][:, 1:]
    R2, t2 = kabsch(vh_m[sl], vh_w[sl])
    return dict(vh_rmsd=rmsd(vh_m, vh_w), vl_rmsd=rmsd(vl_m, vl_w),
                h3_rmsd_pose=rmsd(vh_m[sl], vh_w[sl]),
                h3_rmsd_shape=rmsd(vh_m[sl] @ R2.T + t2, vh_w[sl]))


# --------------------------------------------------------------------------- Boltz-2
def _queries(loops: list[str]) -> dict:
    ref = json.loads(REF_Q.read_text())["queries"][0]
    her2 = [dict(c) for c in ref["chains"] if c["chain_id"] in ("C", "D")]
    assert len(her2) == 2, "reference queries.json does not carry HER2 chains C and D"
    qs = []
    for l in loops:
        qs.append({"query_id": l, "chains": [
            {"chain_id": "A", "sequence": TZ_VL, "msa_path": None},
            {"chain_id": "B", "sequence": vh_with_loop(l), "msa_path": None},
            *[dict(c) for c in her2]]})
    return {"queries": qs}


def _boltz_call(d: Path, seed: int, timeout: int) -> None:
    uid, gid = os.getuid(), os.getgid()
    cmd = ["docker", "run", "--rm", "--gpus", "all", "--shm-size", "16g", "--user", f"{uid}:{gid}",
           "-e", "HOME=/tmp", "-v", "/etc/passwd:/etc/passwd:ro", "-v", f"{d}:/work", "-w", "/work",
           BOLTZ_IMG, "python", "-m", "apheris_boltz.cli", "predict",
           "--input", "in/queries.json", "--output", "out", "--template-paths", "in/template.cif",
           "--write-detailed-scores", "false", "--num-diffusion-samples", "1",
           "--seeds", str(seed), "--query-failures", "skip"]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    (d / "boltz.log").write_text((p.stdout or "") + "\n--- stderr ---\n" + (p.stderr or ""), encoding="utf-8")
    if p.returncode != 0:
        raise RuntimeError(f"boltz exited {p.returncode}; see {d / 'boltz.log'}\n{(p.stderr or '')[-600:]}")


def run_boltz(loops: list[str], tag: str, seed: int = 42, chunk: int = 20,
              timeout: int = 2400) -> pd.DataFrame:
    """Fold each design (plus the trastuzumab control) against HER2 domain IV, template-guided."""
    loops = list(dict.fromkeys(loops))
    rows = []
    for i0 in range(0, len(loops), chunk):
        part = loops[i0 : i0 + chunk]
        d = GPU / f"{tag}_c{i0 // chunk:02d}"
        if (d / "out").exists():
            shutil.rmtree(d / "out")
        (d / "in").mkdir(parents=True, exist_ok=True)
        (d / "in" / "queries.json").write_text(json.dumps(_queries([TZ_LOOP] + part), indent=1))
        shutil.copy(REF_T, d / "in" / "template.cif")
        if not wait_for_gpu():
            break
        t0 = time.time()
        _boltz_call(d, seed, timeout)
        log(f"Boltz-2 chunk {i0 // chunk}: {len(part) + 1} folds in {time.time() - t0:.0f}s")
        wt_model = d / "out" / TZ_LOOP / f"seed_{seed}" / f"{TZ_LOOP}_seed_{seed}_sample_1_model.cif"
        for l in part:
            base = d / "out" / l / f"seed_{seed}"
            sj = base / f"{l}_seed_{seed}_sample_1_scores.json"
            cif = base / f"{l}_seed_{seed}_sample_1_model.cif"
            if not (sj.exists() and cif.exists() and wt_model.exists()):
                continue
            s = json.loads(sj.read_text())
            cc = np.array(s["chain_chain_interface_predicted_tm_score"])
            r = dict(cdrh3=l, iptm=s["interface_predicted_tm_score"], ptm=s["predicted_tm_score"],
                     plddt=s["average_predicted_local_distance_difference_test"],
                     gpde=s["global_predicted_distance_error"],
                     ab_her2_iptm=float(np.mean([cc[0, 2], cc[0, 3], cc[1, 2], cc[1, 3]])),
                     model_cif=str(cif))
            r.update(pose_metrics(cif, wt_model))
            rows.append(r)
    return pd.DataFrame(rows)


def wt_noise_floor(seed_ref: int = 42, seed_alt: int = 43) -> dict:
    """How far does the control H3 move between two seeds? Cached; a one-off 2-fold job."""
    f = GPU / "wt_noise.json"
    if f.exists():
        return json.loads(f.read_text())
    cifs = {}
    for seed in (seed_ref, seed_alt):
        dd = GPU / "wtnoise" / f"s{seed}"
        if (dd / "out").exists():
            shutil.rmtree(dd / "out")
        (dd / "in").mkdir(parents=True, exist_ok=True)
        (dd / "in" / "queries.json").write_text(json.dumps(_queries([TZ_LOOP]), indent=1))
        shutil.copy(REF_T, dd / "in" / "template.cif")
        if not wait_for_gpu():
            return {}
        _boltz_call(dd, seed, 1200)
        cifs[seed] = dd / "out" / TZ_LOOP / f"seed_{seed}" / f"{TZ_LOOP}_seed_{seed}_sample_1_model.cif"
    res = pose_metrics(cifs[seed_alt], cifs[seed_ref])
    res["seeds"] = [seed_ref, seed_alt]
    f.write_text(json.dumps(res, indent=1))
    log(f"control noise floor between seeds: H3 pose {res['h3_rmsd_pose']:.2f} A, "
        f"shape {res['h3_rmsd_shape']:.2f} A, VH {res['vh_rmsd']:.2f} A")
    return res


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    pool = pd.read_csv(HERE / "out/round_00/shortlist.csv").head(n)
    print(run_nanomelt(pool.cdrh3.tolist(), "selftest_nm").round(2).to_string(index=False))


def run_nanomelt_vh(seqs: dict[str, str], tag: str, timeout: int = 1800) -> dict[str, float]:
    """Predicted Tm for arbitrary VH sequences (id -> VH), with the trastuzumab control as id WT."""
    d = GPU / tag
    d.mkdir(parents=True, exist_ok=True)
    fa, out = d / "in.fa", d / "nanomelt.csv"
    ids = {f"d{i:05d}": k for i, k in enumerate(seqs)}
    with fa.open("w") as fh:
        fh.write(f">WT\n{TZ_VH}\n")
        for new, k in ids.items():
            fh.write(f">{new}\n{seqs[k]}\n")
    env = dict(os.environ, PATH=f"{Path(NANOMELT).parent}:{os.environ.get('PATH', '')}")
    p = subprocess.run([NANOMELT, "predict", "-i", str(fa), "-o", str(out), "-align", "-ncpu", "4"],
                       capture_output=True, text=True, timeout=timeout, cwd=d, env=env)
    if not out.exists():
        raise RuntimeError(f"nanomelt produced no output:\n{p.stdout[-600:]}\n{p.stderr[-600:]}")
    df = pd.read_csv(out)
    tm = dict(zip(df["ID"].astype(str), df["NanoMelt Tm (C)"].astype(float)))
    res = {k: tm[new] for new, k in ids.items()}
    res["WT"] = tm["WT"]
    return res


# --------------------------------------------------------------------------- Boltz-2 on exact sequences
def _safe(qid: str) -> str:
    import re
    return re.sub(r"[^A-Za-z0-9_-]", "_", qid)


def _queries_seqs(items: dict) -> dict:
    ref = json.loads(REF_Q.read_text())["queries"][0]
    her2 = [dict(c) for c in ref["chains"] if c["chain_id"] in ("C", "D")]
    qs = []
    for qid, (vh, vl) in items.items():
        qs.append({"query_id": _safe(qid), "chains": [
            {"chain_id": "A", "sequence": vl, "msa_path": None},
            {"chain_id": "B", "sequence": vh, "msa_path": None},
            *[dict(c) for c in her2]]})
    return {"queries": qs}


def run_boltz_seqs(items: dict, tag: str, seed: int = 42, chunk: int = 10, timeout: int = 2400) -> pd.DataFrame:
    """Fold exact designs {id: (VH, VL)} against HER2 domain IV, with a trastuzumab control in every chunk."""
    ids, rows = list(items), []
    for i0 in range(0, len(ids), chunk):
        part = ids[i0 : i0 + chunk]
        d = GPU / f"{tag}_c{i0 // chunk:02d}"
        if (d / "out").exists():
            shutil.rmtree(d / "out")
        (d / "in").mkdir(parents=True, exist_ok=True)
        batch = {"WT": (TZ_VH, TZ_VL), **{k: items[k] for k in part}}
        (d / "in" / "queries.json").write_text(json.dumps(_queries_seqs(batch), indent=1))
        shutil.copy(REF_T, d / "in" / "template.cif")
        if not wait_for_gpu():
            break
        t0 = time.time()
        _boltz_call(d, seed, timeout)
        log(f"Boltz-2 (exact sequences) chunk {i0 // chunk}: {len(part) + 1} folds in {time.time() - t0:.0f}s")
        wt_model = d / "out" / "WT" / f"seed_{seed}" / f"WT_seed_{seed}_sample_1_model.cif"
        for qid in part:
            q = _safe(qid)
            base = d / "out" / q / f"seed_{seed}"
            sj, cif = base / f"{q}_seed_{seed}_sample_1_scores.json", base / f"{q}_seed_{seed}_sample_1_model.cif"
            if not (sj.exists() and cif.exists() and wt_model.exists()):
                continue
            s = json.loads(sj.read_text())
            cc = np.array(s["chain_chain_interface_predicted_tm_score"])
            r = dict(id=qid, iptm=s["interface_predicted_tm_score"], ptm=s["predicted_tm_score"],
                     plddt=s["average_predicted_local_distance_difference_test"],
                     ab_her2_iptm=float(np.mean([cc[0, 2], cc[0, 3], cc[1, 2], cc[1, 3]])), model_cif=str(cif))
            r.update(pose_metrics(cif, wt_model))
            rows.append(r)
    return pd.DataFrame(rows)
