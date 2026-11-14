from __future__ import annotations

import polars as pl

from abpipe.base import Context, Scorer
from abpipe.registry import register


@register("scorer", "pass_through")
class PassThrough(Scorer):
    """Promote an existing column to a score without recomputing anything.

    score__parent_score (from info__parent_score, higher is better) carries the ranked-list score a
    generator inherited for its variants — a prior for ordering, never evidence of binding. It is
    what the de novo VHH lane ranks on until a validated scorer exists for that track
    (NICCOLO/eda/AGENT_GUIDE.md: the pre-filter does not transfer to de novo designs).
    """

    version = "1"
    outputs = {"score__value": "higher"}

    def __init__(self, column: str = "info__parent_score", output: str = "score__parent_score",
                 direction: str = "higher"):
        if direction not in ("higher", "lower"):
            raise ValueError("direction must be 'higher' or 'lower'")
        if not output.startswith("score__"):
            raise ValueError("output must be a score__* column")
        self.column, self.outputs = column, {output: direction}

    def score(self, df: pl.DataFrame, ctx: Context) -> pl.DataFrame:
        if self.column not in df.columns:
            raise ValueError(f"pass_through: column {self.column!r} is not in the table "
                             f"(generator must provide it); found {sorted(df.columns)}")
        return df.select("sequence_id", pl.col(self.column).cast(pl.Float64).alias(
            next(iter(self.outputs))))
