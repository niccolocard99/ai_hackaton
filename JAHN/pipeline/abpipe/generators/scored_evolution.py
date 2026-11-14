from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import polars as pl

from abpipe import _niccolo, data, novelty
from abpipe.base import Context, Generator
from abpipe.registry import register

AA, L, NA = _niccolo.AA, _niccolo.L, 20
_W20 = 20 ** np.arange(L - 1, -1, -1, dtype=np.int64)
_N, _P, _S, _T = (AA.index(c) for c in "NPST")
_OPS = ("sample", "mutate", "cross")


def _z(t: np.ndarray) -> np.ndarray:
    s = float(t.std())
    return (t - t.mean()) / (s if s > 0 else 1.0)


def _softmax(logits: np.ndarray) -> np.ndarray:
    e = np.exp(logits - logits.max(axis=1, keepdims=True))
    return e / e.sum(axis=1, keepdims=True)


def file_sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sha_pinned(path, expected: str | None) -> bool:
    """True when `expected` is absent, or is the file's full sha256, or a >=8-hex-char prefix of it."""
    if not expected:
        return True
    full = file_sha256(path).lower()
    e = str(expected).lower()
    return full == e or (len(e) >= 8 and full.startswith(e))


def history_digest(paths) -> str:
    return hashlib.sha256("".join(file_sha256(p) for p in paths).encode()).hexdigest()


def read_table(path) -> pl.DataFrame:
    path = Path(path)
    return pl.read_parquet(path) if path.suffix == ".parquet" else pl.read_csv(path)


def _topk_mask(rng, rows: int, k, mu: np.ndarray) -> np.ndarray:
    """Row i gets exactly k[i] True entries, drawn without replacement with probability proportional to mu."""
    g = np.log(mu)[None, :] + rng.gumbel(size=(rows, L))
    return np.argsort(np.argsort(-g, axis=1), axis=1) < np.asarray(k)[:, None]


@register("generator", "scored_evolution")
class ScoredEvolution(Generator):
    """Propose CDRH3 windows from the pre-filter, the known datasets and the previous round's best.

    Prior, every round: per-position residue preferences = z-scored classifier weights + z-scored (negated) ranker
    weights from NICCOLO/eda/prefilter.py + z-scored FACS 'high' enrichment (organiser data), softmaxed at
    `temperature`. Feedback, when `parents` is a previous round's ranked table: its top `n_elite` windows pull the
    residue, mutated-position and substitution-count distributions towards themselves with weight `eta`
    (a cross-entropy-method update). Children come from three operators: `sample` (a fresh window from the
    distributions), `mutate` (1-2 edits to an elite) and `cross` (uniform crossover of two elites).
    Every child must have min_subs..max_subs substitutions, no excluded residue, no new N-glycosylation sequon,
    be at least min_dist_provided substitutions from every provided window, and be new (not in `history`).
    """

    version = "1"

    def __init__(self, parents: str | None = None, parents_sha256: str | None = None, history: tuple = (),
                 history_sha256: str | None = None, n_elite: int = 100, frac_sample: float = 0.4,
                 frac_mutate: float = 0.4, frac_cross: float = 0.2, min_subs: int = 4, max_subs: int = 8,
                 w_clf: float = 1.0, w_rank: float = 1.0, w_facs: float = 0.5, temperature: float = 1.0,
                 eta: float = 0.5, min_dist_provided: int = 3, exclude: str = "C",
                 eda_dir: str = _niccolo.EDA_REL, splits_dir: str | None = None, spr_csv: str = _niccolo.SPR_REL):
        if not 1 <= min_subs <= max_subs <= L:
            raise ValueError(f"need 1 <= min_subs <= max_subs <= {L}")
        if min(frac_sample, frac_mutate, frac_cross) < 0 or frac_sample + frac_mutate + frac_cross <= 0:
            raise ValueError("operator fractions must be non-negative and not all zero")
        if not 0.0 <= eta <= 1.0 or temperature <= 0 or n_elite < 1:
            raise ValueError("need 0 <= eta <= 1, temperature > 0, n_elite >= 1")
        if min_dist_provided not in (0, 1, 2, 3):
            raise ValueError("min_dist_provided must be 0, 1, 2 or 3")
        self.parents, self.parents_sha256, self.history = parents, parents_sha256, tuple(history or ())
        self.history_sha256, self.n_elite = history_sha256, n_elite
        self.fracs = np.array([frac_sample, frac_mutate, frac_cross], float)
        self.min_subs, self.max_subs = min_subs, max_subs
        self.w_clf, self.w_rank, self.w_facs, self.temperature, self.eta = w_clf, w_rank, w_facs, temperature, eta
        self.min_dist_provided, self.exclude = min_dist_provided, exclude
        self.eda_dir, self.splits_dir, self.spr_csv = eda_dir, splits_dir, spr_csv
        self.version = f"1+{_niccolo.quick_tag(eda_dir)}"

    # ------------------------------------------------------------------ inputs
    def _elites(self, parent) -> np.ndarray:
        if not self.parents:
            return np.zeros((0, L), np.uint8)
        if self.parents_sha256 and not sha_pinned(self.parents, self.parents_sha256):
            raise ValueError(f"{self.parents} changed since the config was written (sha256 mismatch)")
        df = read_table(self.parents)
        if "rank" in df.columns:
            df = df.sort("rank")
        wins = [w for w in _niccolo.windows(df["heavy_sequence"].to_list(), parent) if w][: self.n_elite]
        return novelty.encode(wins) if wins else np.zeros((0, L), np.uint8)

    def _history(self, parent) -> np.ndarray:
        if self.history_sha256:
            files = [str(p) for p in self.history]
            e = str(self.history_sha256).lower()
            combined = history_digest(files)
            if not (combined == e or (len(e) >= 8 and combined.startswith(e))
                    or all(sha_pinned(p, e) for p in files)):
                raise ValueError("a history file changed since the config was written (sha256 mismatch)")
        codes = []
        for path in self.history:
            wins = [w for w in _niccolo.windows(read_table(path)["heavy_sequence"].to_list(), parent) if w]
            if wins:
                codes.append(novelty.encode(wins).astype(np.int64) @ _W20)
        return np.unique(np.concatenate(codes)) if codes else np.zeros(0, np.int64)

    def _facs_enrichment(self, ctx: Context) -> np.ndarray:
        """log fold-enrichment of each residue at each position among FACS 'high' variants (smoothed)."""
        facs = data.facs_table(ctx)
        S = novelty.encode([w for w in facs["cdrh3"].to_list() if len(w) == L and set(w) <= set(AA)])
        keep = [i for i, w in enumerate(facs["cdrh3"].to_list()) if len(w) == L and set(w) <= set(AA)]
        hi = facs["binding_class"].is_in(["high"]).to_numpy()[keep]
        base = hi.mean()
        E = np.zeros((L, NA))
        for p in range(L):
            a_all = np.bincount(S[:, p], minlength=NA).astype(float)
            a_hi = np.bincount(S[hi, p], minlength=NA).astype(float)
            E[p] = np.log(((a_hi + 2 * base) / (a_all + 2)) / base)
        return E

    # ------------------------------------------------------------------ proposal distributions
    def _state(self, ctx: Context, parent, fit) -> SimpleNamespace:
        wt = novelty.encode([parent.window])[0]
        allowed = np.array([a not in self.exclude for a in AA])
        pref = (self.w_clf * _z(fit.clf_w) + self.w_rank * _z(-fit.rank_w) + self.w_facs * _z(self._facs_enrichment(ctx)))
        logits = np.where(allowed[None, :], pref / self.temperature, -np.inf)
        logits[np.arange(L), wt] = -np.inf                              # proposals are substitutions
        Q0 = _softmax(logits)
        ks = np.arange(self.min_subs, self.max_subs + 1)
        kp0 = 0.6 ** (ks - self.min_subs)
        kp0 = kp0 / kp0.sum()
        elites = self._elites(parent)
        E = len(elites)
        if E:
            f = np.stack([np.bincount(elites[:, p], minlength=NA) for p in range(L)]).astype(float)
            f[np.arange(L), wt] = 0.0
            f = (f + 0.5) * allowed[None, :]
            f[np.arange(L), wt] = 0.0
            f /= f.sum(axis=1, keepdims=True)
            Q = (1 - self.eta) * Q0 + self.eta * f
            mu = (1 - self.eta) * np.full(L, 0.5) + self.eta * (((elites != wt).sum(0) + 1) / (E + 2))
            hist = np.bincount((elites != wt).sum(1), minlength=L + 1)[ks].astype(float) + 0.5
            kp = (1 - self.eta) * kp0 + self.eta * hist / hist.sum()
        else:
            Q, mu, kp = Q0, np.full(L, 0.5), kp0
        cdf = np.cumsum(Q, axis=1)
        cdf /= cdf[:, -1:]
        s, e = parent.window_start, parent.window_start + L
        flank = lambda t: np.array([AA.index(c) for c in t], np.uint8)
        left, right = flank(parent.heavy[s - 2:s]), flank(parent.heavy[e:e + 3])
        return SimpleNamespace(wt=wt, Q=Q, cdf=cdf, last=np.array([np.flatnonzero(Q[p] > 0).max() for p in range(L)]),
                               mu=mu, ks=ks, kp=kp / kp.sum(), elites=elites, left=left, right=right,
                               parent_seq=self._sequons(np.concatenate([left, wt, right])[None, :])[0],
                               excl=np.array([AA.index(c) for c in self.exclude], np.uint8),
                               index=novelty.provided_index(ctx) if self.min_dist_provided else None)

    @staticmethod
    def _sequons(ext: np.ndarray) -> np.ndarray:
        """N-X-S/T (X not P) starting at each position of the extended window."""
        return (ext[:, :-2] == _N) & (ext[:, 1:-1] != _P) & ((ext[:, 2:] == _S) | (ext[:, 2:] == _T))

    # ------------------------------------------------------------------ operators
    @staticmethod
    def _draw(rng, st, p: int, size: int) -> np.ndarray:
        idx = np.searchsorted(st.cdf[p], rng.random(size), side="right")
        return np.minimum(idx, st.last[p]).astype(np.uint8)

    def _fill(self, rng, st, S: np.ndarray, mask: np.ndarray) -> np.ndarray:
        for p in range(L):
            idx = np.flatnonzero(mask[:, p])
            if idx.size:
                S[idx, p] = self._draw(rng, st, p, idx.size)
        return S

    def _op_sample(self, rng, st, rows: int) -> np.ndarray:
        k = rng.choice(st.ks, size=rows, p=st.kp)
        return self._fill(rng, st, np.tile(st.wt, (rows, 1)), _topk_mask(rng, rows, k, st.mu))

    def _op_mutate(self, rng, st, rows: int) -> np.ndarray:
        S = st.elites[rng.integers(0, len(st.elites), rows)].copy()
        m = rng.choice([1, 2], size=rows, p=[0.7, 0.3])
        return self._fill(rng, st, S, _topk_mask(rng, rows, m, st.mu))

    def _op_cross(self, rng, st, rows: int) -> np.ndarray:
        a = st.elites[rng.integers(0, len(st.elites), rows)]
        b = st.elites[rng.integers(0, len(st.elites), rows)]
        return np.where(rng.random((rows, L)) < 0.5, a, b)

    def _accept(self, S: np.ndarray, st) -> np.ndarray:
        nm = (S != st.wt).sum(1)
        ok = (nm >= self.min_subs) & (nm <= self.max_subs)
        if len(st.excl):
            ok &= ~np.isin(S, st.excl).any(1)
        ext = np.concatenate([np.tile(st.left, (len(S), 1)), S, np.tile(st.right, (len(S), 1))], axis=1)
        ok &= ~(self._sequons(ext) & ~st.parent_seq[None, :]).any(1)
        if st.index is not None and ok.any():
            idx = np.flatnonzero(ok)
            ok[idx] &= ~st.index.within(S[idx], self.min_dist_provided - 1)
        return ok

    # ------------------------------------------------------------------ generate
    def generate(self, n: int, seed: int, ctx: Context) -> pl.DataFrame:
        parent = ctx.require_parent()
        fit = _niccolo.get_prefilter(data_dir=ctx.data_dir, cache_dir=ctx.cache_dir, repo_root=ctx.repo_root,
                                     eda_dir=self.eda_dir, splits_dir=self.splits_dir, spr_csv=self.spr_csv)
        _niccolo.check_parent(parent, fit)
        rng = np.random.default_rng(seed)
        st = self._state(ctx, parent, fit)
        fr = self.fracs.copy()
        if not len(st.elites):
            fr[1:] = 0.0                                   # nothing to mutate or cross yet: cold start
        fr = fr / fr.sum()
        ops = (self._op_sample, self._op_mutate, self._op_cross)
        taken = self._history(parent)
        quota = np.floor(n * fr).astype(int)
        quota[0] += n - quota.sum()                        # the rounding remainder goes to fresh samples
        have = np.zeros(3, int)
        kept_S, kept_op = [], []

        def draw(i: int, want: int) -> None:
            nonlocal taken
            S = ops[i](rng, st, int(max(5_000, 6 * want)))
            S = S[self._accept(S, st)]
            codes = S.astype(np.int64) @ _W20
            _, first = np.unique(codes, return_index=True)
            first.sort()
            S, codes = S[first], codes[first]
            fresh = ~np.isin(codes, taken)
            S, codes = S[fresh][:want], codes[fresh][:want]
            taken = np.union1d(taken, codes)
            kept_S.append(S)
            kept_op.append(np.full(len(S), i))
            have[i] += len(S)

        for _ in range(30):                                # every operator fills its own quota
            for i in range(3):
                if quota[i] > have[i]:
                    draw(i, int(quota[i] - have[i]))
            if (have >= quota).all():
                break
        for _ in range(30):                                # a shortfall is topped up with fresh samples
            if have.sum() >= n:
                break
            draw(0, int(n - have.sum()))
        got = int(have.sum())
        if got < n:
            print(f"warning: {self.name} produced {got} of {n} requested unique windows")
        S = np.vstack(kept_S) if kept_S else np.zeros((0, L), np.uint8)
        lab = np.concatenate(kept_op) if kept_op else np.zeros(0, int)
        perm = rng.permutation(len(S))                     # do not leave the operators in blocks
        S, lab = S[perm], lab[perm]
        op_names = [_OPS[i] for i in lab]
        s, e = parent.window_start, parent.window_start + L
        rows = []
        for i, (w, op) in enumerate(zip(novelty.decode(S), op_names)):
            subs = [f"{parent.window[p]}{p + 1}{w[p]}" for p in range(L) if w[p] != parent.window[p]]
            rows.append({"sequence_id": f"{self.name}_s{seed}_{i:06d}", "kind": "paired",
                         "heavy_sequence": parent.heavy[:s] + w + parent.heavy[e:], "light_sequence": parent.light,
                         "design_description": f"{len(subs)} substitutions in the CDRH3 window of {parent.id} ({op})",
                         "parent_id": parent.id, "origin": f"{self.name}:{op}", "seed": seed,
                         "meta": json.dumps({"window": w, "subs": subs, "op": op})})
        return pl.DataFrame(rows, schema={"sequence_id": pl.String, "kind": pl.String, "heavy_sequence": pl.String,
                                          "light_sequence": pl.String, "design_description": pl.String,
                                          "parent_id": pl.String, "origin": pl.String, "seed": pl.Int64,
                                          "meta": pl.String})
