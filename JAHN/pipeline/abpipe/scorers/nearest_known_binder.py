from __future__ import annotations

import numpy as np
import polars as pl
from rapidfuzz import process
from rapidfuzz.distance import Hamming, Levenshtein

from abpipe import data
from abpipe.base import Context, Scorer
from abpipe.distance import check_metric
from abpipe.registry import register


@register("scorer", "nearest_known_binder")
class NearestKnownBinder(Scorer):
    """Closeness of the CDRH3 window to FACS-measured variants of the chosen class(es).

    score__binder_distance: distance to the nearest known binder (lower is better).
    score__neighbor_binder_frac: share of binders among the k nearest FACS variants, counting every
    variant tied with the k-th so the result does not depend on the order of the reference table.
    Binders are FACS class high or medium — the team's rule: {high, medium} are binders, everything
    else is a non-binder. The classes are not an affinity ranking, only a binder/non-binder split.
    metric auto uses Hamming (substituted positions, as the FACS edit_distance counts them) because
    every candidate window and every FACS variant has the same length; it falls back to Levenshtein
    if the FACS lengths differ from the window's.
    """

    version = "4"
    outputs = {"score__binder_distance": "lower", "score__neighbor_binder_frac": "higher"}

    def __init__(self, binder_classes: tuple = ("high", "medium"), k: int = 25, chunk: int = 256, metric: str = "auto"):
        self.binder_classes, self.k, self.chunk = tuple(binder_classes), k, chunk
        self.metric = check_metric(metric)

    def score(self, df: pl.DataFrame, ctx: Context) -> pl.DataFrame:
        p = ctx.require_parent()
        s, length = p.window_start, len(p.window)
        win = df.select(pl.when(pl.col("heavy_sequence").str.len_chars() == len(p.heavy))
                        .then(pl.col("heavy_sequence").str.slice(s, length)).otherwise(None).alias("w"))["w"].to_list()
        facs = data.facs_table(ctx)
        refs = facs["cdrh3"].to_list()
        same_len = set(facs["cdrh3"].str.len_chars().unique().to_list()) == {length}
        if self.metric == "hamming" and not same_len:
            raise ValueError("metric 'hamming' needs every FACS variant to be as long as the window")
        scorer = Hamming.distance if (same_len and self.metric != "levenshtein") else Levenshtein.distance
        is_b = facs["binding_class"].is_in(list(self.binder_classes)).to_numpy()
        ok = [i for i, w in enumerate(win) if w is not None]
        dist, frac = np.full(len(win), np.nan), np.full(len(win), np.nan)
        for lo in range(0, len(ok), self.chunk):
            idx = ok[lo:lo + self.chunk]
            d = process.cdist([win[i] for i in idx], refs, scorer=scorer, dtype=np.uint8, workers=-1)
            dist[idx] = d[:, is_b].min(axis=1)
            kth = np.partition(d, self.k - 1, axis=1)[:, self.k - 1]       # distance of the k-th nearest
            near = d <= kth[:, None]                                       # the k nearest plus every tie at that distance
            frac[idx] = (near & is_b[None, :]).sum(axis=1) / near.sum(axis=1)
        return pl.DataFrame({"sequence_id": df["sequence_id"],
                             "score__binder_distance": pl.Series(dist).fill_nan(None),
                             "score__neighbor_binder_frac": pl.Series(frac).fill_nan(None)})
