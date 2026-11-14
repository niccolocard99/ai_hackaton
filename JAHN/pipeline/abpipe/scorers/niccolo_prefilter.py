from __future__ import annotations

import numpy as np
import polars as pl

from abpipe import _niccolo
from abpipe.base import Context, Scorer
from abpipe.registry import register


@register("scorer", "niccolo_prefilter")
class NiccoloPrefilter(Scorer):
    """Niccolo's HER2 CDRH3 pre-filter (NICCOLO/eda/prefilter.py): P(measurable KD) and a predicted-KD ranking.

    score__p_binder: logistic regression on the one-hot CDRH3 window, trained on 620 measured SPR outcomes
    (higher is better). Gate on it; it is a ranking signal, not a calibrated probability for novel designs.
    score__pred_log10_kd: ridge regression to log10 KD in nM, trained on the 179 variants with a measured KD
    (lower is tighter). A RANKING, not a predicted KD: it compresses at the tight end.
    Scope: trastuzumab's 10-residue CDRH3 window on the trastuzumab framework (optimisation track only). It says
    nothing about thermostability or epitope preservation. Rows it cannot score (framework changed, non-standard
    residue) get null scores. The FACS score is context only and measured to add nothing (+0.671 with it, +0.673 without), so prefilter.py
    now defaults to use_facs=False and emits no `p_facs_high`; emit_facs then raises rather than returning nulls.
    """

    version = "1"
    outputs = {"score__p_binder": "higher", "score__pred_log10_kd": "lower"}

    def __init__(self, eda_dir: str = _niccolo.EDA_REL, splits_dir: str | None = None,
                 spr_csv: str = _niccolo.SPR_REL, emit_facs: bool = False):
        self.eda_dir, self.splits_dir, self.spr_csv, self.emit_facs = eda_dir, splits_dir, spr_csv, emit_facs
        self.outputs = dict(type(self).outputs)
        if emit_facs:
            self.outputs["score__p_facs_high"] = "higher"
        self.version = f"1+{_niccolo.quick_tag(eda_dir)}"      # cache key follows the model's source

    def score(self, df: pl.DataFrame, ctx: Context) -> pl.DataFrame:
        fit = _niccolo.get_prefilter(data_dir=ctx.data_dir, cache_dir=ctx.cache_dir, repo_root=ctx.repo_root,
                                     eda_dir=self.eda_dir, splits_dir=self.splits_dir, spr_csv=self.spr_csv)
        parent = ctx.require_parent()
        _niccolo.check_parent(parent, fit)
        if self.emit_facs and not _niccolo.has_facs(fit):
            raise _niccolo.PrefilterUnavailable(
                "emit_facs=True but this pre-filter carries no FACS model. prefilter.py defaults to "
                "use_facs=False, because the FACS feature was measured to add nothing to the ranker "
                "(+0.671 with it, +0.673 without; NICCOLO/eda/AGENT_GUIDE.md section 3). Drop emit_facs, or "
                "fit with PreFilter(use_facs=True) and the clonotype splits.")
        res = _niccolo.score_windows(fit, _niccolo.windows(df["heavy_sequence"].to_list(), parent))
        cols = {"sequence_id": df["sequence_id"],
                "score__p_binder": pl.Series(res["p_binder"]).fill_nan(None),
                "score__pred_log10_kd": pl.Series(res["pred_log10_kd"]).fill_nan(None)}
        if self.emit_facs:
            cols["score__p_facs_high"] = pl.Series(res["p_facs_high"]).fill_nan(None)
        return pl.DataFrame(cols)
