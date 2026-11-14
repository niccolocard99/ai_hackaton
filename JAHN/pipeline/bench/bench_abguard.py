"""Measure where abguard time goes: numbering vs novelty cross-join. Read-only."""
import sys, time
import polars as pl

GEN = "/home/hackathon/nick-1m/runs/pilot50k/pilot50k/01_generated.parquet"
df = pl.read_parquet(GEN)
print("generated cols:", df.columns, "rows:", df.height, flush=True)

from abguard.annotate import annotate
from abguard.config import PipelineConfig
from abguard import reference as refmod

need = ["sequence_id", "kind", "heavy_sequence", "light_sequence", "design_description"]
sub = df.select([c for c in need if c in df.columns])
if "light_sequence" in sub.columns:
    sub = sub.with_columns(pl.col("light_sequence").fill_null(""))
print("submission cols:", sub.columns, flush=True)

def bench(n, reference=None, gate=False, label=""):
    frame = sub.head(n)
    t = time.perf_counter()
    out = annotate(frame, PipelineConfig(), include_liabilities=False,
                   reference=reference, gate_novelty=gate)
    dt = time.perf_counter() - t
    print(f"{label:28s} n={n:7d}  {dt:8.2f}s  {n/dt:10.1f} designs/s", flush=True)
    return dt

# warm up (imports, lazy init)
bench(50, label="warmup numbering")
print("--- NUMBERING ONLY (abguard_hard path) ---", flush=True)
for n in (1000, 10000, 50000):
    bench(n, label="numbering only")

print("--- WITH NOVELTY GATE (track_check path) ---", flush=True)
bp = refmod.bundled_reference()
print("bundled reference:", bp, flush=True)
ref = refmod.cdrh3_reference(bp) if bp else None
if ref is not None:
    nref = ref.select(pl.len()).collect().item()
    print("reference rows:", nref, flush=True)
    for n in (100, 500, 2000):
        bench(n, reference=ref, gate=True, label="numbering+novelty")
