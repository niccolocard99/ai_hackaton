"""Shared loader for Niccolo's HER2 CDRH3 pre-filter (NICCOLO/eda/prefilter.py).

The scorer, the gate filter and the generator all need the same fitted model, so it is fitted once per process.
Nothing here re-implements the model: prefilter.py is imported from the repo and used as NICCOLO/eda/AGENT_GUIDE.md
documents it. This file only supplies what that script expects from its working directory (it reads `data/...`
relative to the cwd) and turns abpipe candidates into the 10-residue windows it scores.
"""
from __future__ import annotations

import hashlib
import importlib.util
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"
L = 10                                        # the window length prefilter.py is built for
EDA_REL = "NICCOLO/eda"
SPR_REL = "her2-trastuzumab/trastuzumab_her2_spr_all_768.csv"
_SPLITS_SEARCH = ("~/nick-work/splits", "~/nick-work/optim/eda/data/splits")
_AA_SET = frozenset(AA)
_FITTED: dict = {}


class PrefilterUnavailable(RuntimeError):
    """The pre-filter, its data or a dependency is missing. The message says what and how to fix it."""


def _sha(path: Path, n: int = 10) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:n]


def repo_guess() -> Path:
    """<repo>/JAHN/pipeline/abpipe/_niccolo.py -> <repo>."""
    return Path(__file__).resolve().parents[3]


def resolve_eda(eda_dir: str, repo_root=None) -> Path:
    p = Path(os.path.expanduser(str(eda_dir)))
    if p.is_absolute():
        return p
    for root in (repo_root, repo_guess(), Path.cwd()):
        if root is not None and (Path(root) / p / "prefilter.py").exists():
            return Path(root) / p
    return Path(repo_root or repo_guess()) / p


def find_splits(explicit=None, required: bool = True):
    """The FACS clonotype splits (train.csv, val.csv) behind prefilter.py's optional FACS feature.

    Only needed for `use_facs=True`. prefilter.py defaults to False, so with `required=False` a missing
    splits directory returns None instead of raising, and the pre-filter runs off the SPR table alone.
    """
    for cand in (explicit, os.environ.get("PREFILTER_SPLITS_DIR"), *_SPLITS_SEARCH):
        if not cand:
            continue
        p = Path(os.path.expanduser(str(cand)))
        if (p / "train.csv").exists() and (p / "val.csv").exists():
            return p.resolve()
    if not required:
        return None
    raise PrefilterUnavailable(
        "FACS clonotype splits (train.csv, val.csv) not found. Pass `splits_dir` to the module or set "
        "PREFILTER_SPLITS_DIR; on the team VM they are in ~/nick-work/splits (NICCOLO/eda/AGENT_GUIDE.md, section 2).")


def quick_tag(eda_dir: str, repo_root=None) -> str:
    """Fingerprint of prefilter.py, put in module versions so a changed model never serves stale cached scores."""
    f = resolve_eda(eda_dir, repo_root) / "prefilter.py"
    return _sha(f, 8) if f.exists() else "missing"


@dataclass
class Fitted:
    pf: object              # prefilter.PreFilter, exactly as prefilter.py fits it
    mod: object             # the imported prefilter module
    tag: str
    n_facs_dropped: int     # FACS rows removed so the FACS model never saw an SPR sequence
    clf_w: np.ndarray       # (L, 20) logistic weights; row = window position, column = residue in AA order
    rank_w: np.ndarray      # (L, 20) ridge weights on log10 KD (negative = tighter binding)
    wt_window: str          # trastuzumab's window as prefilter.py defines it


def _link(src: Path, dst: Path) -> None:
    if dst.is_symlink() or dst.exists():
        if dst.is_symlink() and os.readlink(dst) == str(src):
            return
        dst.unlink()
    dst.symlink_to(src)


def _import(src: Path):
    spec = importlib.util.spec_from_file_location("niccolo_prefilter", src)
    mod = importlib.util.module_from_spec(spec)
    cwd = os.getcwd()
    os.chdir(src.parent)        # prefilter.py makes ./models at import; keep that inside NICCOLO/eda (git-ignored)
    try:
        spec.loader.exec_module(mod)
    finally:
        os.chdir(cwd)
    return mod


def get_prefilter(*, data_dir, cache_dir, repo_root=None, eda_dir: str = EDA_REL, splits_dir=None,
                  spr_csv: str = SPR_REL, use_facs: bool = False) -> Fitted:
    src = resolve_eda(eda_dir, repo_root) / "prefilter.py"
    if not src.exists():
        raise PrefilterUnavailable(f"{src} not found. It lives on origin/main under NICCOLO/eda; `git pull` first.")
    spr = (Path(data_dir) / spr_csv).resolve()
    if not spr.exists():
        raise PrefilterUnavailable(f"SPR table not found: {spr}")
    # the splits are only needed for use_facs=True, which prefilter.py no longer defaults to
    splits = find_splits(splits_dir, required=use_facs)
    sizes = ".".join(str((splits / f).stat().st_size) for f in ("train.csv", "val.csv")) if splits else "nofacs"
    tag = f"{_sha(src, 8)}.{_sha(spr, 8)}.{sizes}"
    if tag in _FITTED:
        return _FITTED[tag]
    stage = Path(cache_dir) / "niccolo_prefilter_data"
    stage.mkdir(parents=True, exist_ok=True)
    _link(spr, stage / "trastuzumab_her2_spr_all_768.csv")
    if splits:
        _link(splits, stage / "splits")
    try:
        mod = _import(src)
        mod.DATA = stage                              # prefilter.load() reads DATA/<file>
        spr_df, facs_df = mod.load(with_facs=use_facs) if use_facs else mod.load()
        pf = mod.PreFilter(use_facs=use_facs).fit(spr_df, facs_df) if use_facs else mod.PreFilter().fit(spr_df)
    except ImportError as e:
        raise PrefilterUnavailable(f"prefilter.py needs pandas and scikit-learn in this environment ({e}); "
                                   "install with `uv pip install pandas scikit-learn`.") from e
    if mod.AA != AA or mod.L != L:
        raise PrefilterUnavailable("prefilter.py changed its alphabet or window length; update abpipe/_niccolo.py")
    fit = Fitted(pf=pf, mod=mod, tag=tag, n_facs_dropped=int(pf.n_facs_dropped),
                 clf_w=np.asarray(pf.clf.coef_, float).reshape(L, 20),
                 rank_w=np.asarray(pf.rank.coef_, float)[: L * 20].reshape(L, 20), wt_window=mod.TZ)
    _FITTED[tag] = fit
    return fit


def check_parent(parent, fit: Fitted) -> None:
    if parent.window != fit.wt_window:
        raise PrefilterUnavailable(
            f"the pre-filter is specific to trastuzumab's window {fit.wt_window!r}; this run's parent window is "
            f"{parent.window!r}. It does not transfer to other antibodies or to de novo designs.")


def windows(heavy: list, parent) -> list:
    """The parent's mutable window out of each heavy chain; None when the chain cannot be scored by prefilter.py
    (different length from the parent, or a residue outside the 20 standard amino acids)."""
    s, e = parent.window_start, parent.window_start + len(parent.window)
    out = []
    for h in heavy:
        w = h[s:e] if h is not None and len(h) == len(parent.heavy) else None
        out.append(w if w is not None and len(w) == L and _AA_SET.issuperset(w) else None)
    return out


SCORE_KEYS = ("p_binder", "pred_log10_kd", "p_facs_high")


def score_windows(fit: Fitted, wins: list, chunk: int = 200_000) -> dict:
    """prefilter.PreFilter.score on every scorable window; NaN where the window is None.

    The keys follow what the fitted model actually emits, which is not fixed: `p_facs_high` is present only
    when prefilter.py was fitted with `use_facs=True`. It defaults to False upstream, because the FACS
    feature was measured to add nothing to the ranker (+0.671 with it, +0.673 without) while costing an
    11 MB data dependency. Callers must treat that key as optional -- use `has_facs(fit)`.
    """
    n = len(wins)
    ok = np.array([i for i, w in enumerate(wins) if w is not None], dtype=int)
    out = {k: np.full(n, np.nan) for k in ("p_binder", "pred_log10_kd")}
    keys = None
    for lo in range(0, len(ok), chunk):
        idx = ok[lo:lo + chunk]
        res = fit.pf.score([wins[i] for i in idx])
        if keys is None:
            keys = [k for k in SCORE_KEYS if k in res.columns]
            out = {k: np.full(n, np.nan) for k in keys}
        for k in keys:
            out[k][idx] = res[k].to_numpy()
    return out


def has_facs(fit: Fitted) -> bool:
    """Whether this fitted model carries the FACS feature, i.e. emits `p_facs_high`."""
    return getattr(fit.pf, "facs", None) is not None
