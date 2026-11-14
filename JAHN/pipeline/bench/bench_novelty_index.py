"""Index build + 1M-query throughput, against the measured abguard baseline."""
import os, random, sys, time
import numpy as np
import polars as pl
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from abpipe.novelty_index import NoveltyIndex

ref_path = sys.argv[1]
n_query = int(sys.argv[2]) if len(sys.argv) > 2 else 1_000_000
d = int(sys.argv[3]) if len(sys.argv) > 3 else 1

refs = pl.read_parquet(ref_path).get_column("cdrh3").to_list()
print(f"reference loops: {len(refs):,}")

TZ = "SRWGGDGFYAMDY"          # trastuzumab CDRH3; the mutable window is positions 2..11
AA = "ACDEFGHIKLMNPQRSTVWY"
rng = random.Random(11)
queries = []
for _ in range(n_query):
    chars = list(TZ)
    for p in rng.sample(range(2, 12), rng.randint(4, 8)):
        chars[p] = rng.choice(AA)
    queries.append("".join(chars))
print(f"queries: {len(queries):,} ({len(set(queries)):,} distinct)")

t = time.perf_counter()
index = NoveltyIndex.build(refs, max_distance=d)
index.min_distance([TZ])                      # force the lazy bucket build
build = time.perf_counter() - t
print(f"build (max_distance={d}): {build:.1f}s")

t = time.perf_counter()
dist = index.min_distance(queries)
query = time.perf_counter() - t
rate = len(queries) / query
print(f"query: {query:.1f}s  ->  {rate:,.0f} designs/s")
print(f"novel (min_distance >= 2): {int((dist >= 2).sum()):,} / {len(queries):,}")
print(f"distance histogram: {dict(zip(*[x.tolist() for x in np.unique(dist, return_counts=True)]))}")

ABGUARD = 37.2   # measured: numbering+novelty, 4 polars threads, niced
print(f"\nabguard novelty baseline: {ABGUARD} designs/s -> {len(queries)/ABGUARD/3600:.2f} h for {len(queries):,}")
print(f"this index:               {rate:,.0f} designs/s -> {query/60:.2f} min (+{build:.0f}s build)")
print(f"speedup on the novelty stage: {rate/ABGUARD:,.0f}x")
