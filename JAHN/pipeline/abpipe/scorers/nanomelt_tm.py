"""NanoMelt VH Tm as a scorer (LEO handoff item 2): the first non-additive score in the loop.

NanoMelt predicts apparent melting temperatures (semi-supervised ensemble, trained on nanobody
Tm data). The loop's pre-filter is additive over one-hot CDRH3 positions; Tm is not, which is
exactly the regime where extra rounds earn their keep. NICCOLO/thermo measured its error at
about 4.1 degC and it is trained on nanobodies, so the parent VH is out of domain: it is a
TIE-BREAK between comparable designs, never a gate.

score__tm_pred: predicted Tm of the design's VH (higher is better).
score__dtm_vs_wt: the same minus the parent VH's predicted Tm (higher is better).

The binary lives in ~/nick-work/nmenv/bin (NICCOLO/thermo/vm_install_nanomelt.sh). The parent is
folded in as a control row on every call, so the delta is always against the same baseline.
"""
from __future__ import annotations

import subprocess
import time
from pathlib import Path

import polars as pl

from abpipe.base import Context, Scorer
from abpipe.registry import register

DEFAULT_BIN = "~/nick-work/nmenv/bin/nanomelt"


@register("scorer", "nanomelt_tm")
class NanoMeltTm(Scorer):
    """Predicted VH melting temperature from NanoMelt — tie-break only, never a gate."""

    version = "1"
    outputs = {"score__tm_pred": "higher", "score__dtm_vs_wt": "higher"}

    def __init__(self, nanomelt_bin: str = DEFAULT_BIN, ncpu: int = 4, timeout_s: int = 1800,
                 do_align: bool = False):
        self.bin = Path(nanomelt_bin).expanduser()
        if not self.bin.exists():
            raise FileNotFoundError(f"{self.bin} not found; install with NICCOLO/thermo/vm_install_nanomelt.sh")
        self.ncpu, self.timeout_s, self.do_align = ncpu, timeout_s, do_align

    def score(self, df: pl.DataFrame, ctx: Context) -> pl.DataFrame:
        parent = ctx.require_parent()
        d = ctx.run_dir / "_nanomelt"
        d.mkdir(parents=True, exist_ok=True)
        fa, out = d / "in.fa", d / "nanomelt.csv"
        ids = {}
        with fa.open("w") as fh:
            fh.write(f">__wt__\n{parent.heavy}\n")
            for sid, h in zip(df["sequence_id"].to_list(), df["heavy_sequence"].to_list()):
                if h and h != parent.heavy and sid not in ids:
                    ids[sid] = h
            for sid, h in ids.items():
                fh.write(f">{sid}\n{h}\n")
        if not ids:
            return df.select("sequence_id",
                             pl.lit(None, dtype=pl.Float64).alias("score__tm_pred"),
                             pl.lit(None, dtype=pl.Float64).alias("score__dtm_vs_wt"))
        t0 = time.time()
        cmd = [str(self.bin), "predict", "-i", str(fa), "-o", str(out), "-ncpu", str(self.ncpu)]
        if self.do_align:
            cmd.append("-align")
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=self.timeout_s, cwd=d)
        if not out.exists():
            raise RuntimeError(f"nanomelt produced no output (exit {p.returncode}):\n"
                               f"{(p.stdout or '')[-600:]}\n{(p.stderr or '')[-600:]}")
        res = pl.read_csv(out)
        tm = dict(zip(res["ID"].astype(str).to_list(), res["NanoMelt Tm (C)"].astype(float).to_list()))
        wt = tm.get("__wt__")
        rows = [{"sequence_id": sid, "score__tm_pred": tm.get(sid),
                 "score__dtm_vs_wt": (tm[sid] - wt) if (tm.get(sid) is not None and wt is not None) else None}
                for sid in df["sequence_id"].to_list()]
        print(f"nanomelt: {len(ids)} designs + control in {time.time() - t0:.0f}s; "
              f"control Tm {wt:.2f} C" if wt is not None else f"nanomelt: {len(ids)} designs")
        return pl.DataFrame(rows, schema={"sequence_id": pl.String, "score__tm_pred": pl.Float64,
                                          "score__dtm_vs_wt": pl.Float64})
