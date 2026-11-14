#!/usr/bin/env python3
"""De novo validity gate wrapper: the official checker plus our own safety margin.

Nothing in here re-implements the competition rules. The verdict that counts is
produced by shelling out to the official track script

    NICCOLO/scripts/check_denovo.py <designs.csv> -o <report.csv> --strict-exit

and every official column is carried through verbatim. What this module adds is
headroom: the official bars are identity < 0.70 and CDRH3 >= 2 edits from
trastuzumab and from every CDRH3 in abguard's bundled reference, which is exactly
the line a design must not sit on. ``internal_pass`` applies a stricter margin
(identity <= 0.62, CDRH3 >= 4 edits on both counts) so a candidate that clears it
still clears the official gate after a refold, a re-numbering, or a reference
bundle refresh.

CLI::

    python gate.py designs.csv -o out.csv

Library::

    from gate import run_gate
    df = run_gate("designs.csv")          # -> pandas.DataFrame

``designs.csv`` is the submission format: ``sequence_id,kind,heavy_sequence,
light_sequence,design_description`` with ``kind`` in ``vhh|scfv|paired``. A VHH or
an scFv goes whole into ``heavy_sequence`` with ``light_sequence`` empty.

Two interpreters are involved and that is deliberate: the official checker needs
the environment abguard is installed in, while this module needs pandas, and no
single venv on this box has both. The abguard interpreter is only ever reached
through ``subprocess``; the CDRH3 reference is cached locally as a compact npz the
first time it is needed, after which this module needs nothing but numpy + pandas.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from glob import glob
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
OFFICIAL_CHECKER = REPO / "NICCOLO" / "scripts" / "check_denovo.py"

#: Trastuzumab's IMGT CDRH3, as abguard's immunum segmentation reads it.
TRASTUZUMAB_CDRH3 = "SRWGGDGFYAMDY"

#: The official de novo identity bar: a design is `de_novo` only below this.
IDENTITY_CAP = 0.70
#: The official minimum CDRH3 edit distance for a `de_novo` design.
OFFICIAL_CDRH3_DIST_MIN = 2

#: Our stricter internal bars, well inside the official ones.
INTERNAL_IDENTITY_MAX = 0.62
INTERNAL_CDRH3_DIST_MIN = 4

_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


# --------------------------------------------------------------------------- #
# interpreter plumbing
# --------------------------------------------------------------------------- #
def _candidate_pythons(*, needs_abguard: bool) -> list[str]:
    """Interpreters worth trying, best first."""
    out: list[str] = []
    if needs_abguard:
        out.append(os.environ.get("ABGUARD_PY", ""))
        exe = shutil.which("abguard")
        if exe:
            out.append(str(Path(exe).resolve().parent / "python"))
        out.append("/home/hackathon/.venvs/abguard/bin/python")
        # newest uv-cache env carrying the novelty gate (older builds lack --reference-db)
        pattern = "~/.cache/uv/archive-v0/*/lib/python3*/site-packages/abguard/novelty.py"
        found = glob(os.path.expanduser(pattern))
        for novelty in sorted(found, key=os.path.getmtime, reverse=True):
            out.append(str(Path(novelty).parents[4] / "bin" / "python"))
    else:
        out.append(os.environ.get("PANDAS_PY", ""))
        out.append(sys.executable)
        out.extend(sorted(glob("/home/hackathon/.venvs/*/bin/python")))
    seen: set[str] = set()
    return [p for p in out if p and p not in seen and not seen.add(p) and os.path.exists(p)]


def abguard_python() -> str:
    """The interpreter that can run the official checker.

    Resolved rather than hard-coded so a different box, or ``ABGUARD_PY``, still
    works. Verified by actually importing abguard, not by guessing from a path.
    """
    for py in _candidate_pythons(needs_abguard=True):
        probe = subprocess.run(
            [py, "-c", "import abguard"], capture_output=True, text=True, timeout=120
        )
        if probe.returncode == 0:
            return py
    raise RuntimeError(
        "no interpreter with abguard found; set ABGUARD_PY=/path/to/python that has it"
    )


def _reexec_under_pandas() -> None:
    """CLI only: restart under an interpreter that has pandas, if this one has not."""
    if os.environ.get("_DENOVO_GATE_REEXEC") == "1":
        sys.exit(
            "gate.py needs pandas and numpy. Run it with an interpreter that has them, "
            "e.g. /home/hackathon/.venvs/abmodels/bin/python, or set PANDAS_PY."
        )
    for py in _candidate_pythons(needs_abguard=False):
        probe = subprocess.run(
            [py, "-c", "import pandas, numpy"], capture_output=True, text=True, timeout=120
        )
        if probe.returncode == 0:
            os.environ["_DENOVO_GATE_REEXEC"] = "1"
            os.execv(py, [py, *sys.argv])
    sys.exit("no interpreter with pandas + numpy found; set PANDAS_PY=/path/to/python")


if __name__ == "__main__":  # pragma: no cover - import-time interpreter check
    try:
        import pandas  # noqa: F401
        import numpy  # noqa: F401
    except ImportError:
        _reexec_under_pandas()

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


# --------------------------------------------------------------------------- #
# Levenshtein
# --------------------------------------------------------------------------- #
def levenshtein(a: str, b: str) -> int:
    """Plain Levenshtein, same recurrence the official scripts use."""
    if a == b:
        return 0
    if not a:
        return len(b)
    if not b:
        return len(a)
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


# --------------------------------------------------------------------------- #
# the CDRH3 reference, cached
# --------------------------------------------------------------------------- #
def cache_dir() -> Path:
    return Path(os.environ.get("DENOVO_GATE_CACHE", str(HERE / "cache")))


_CACHE_VERSION = 1
_REFERENCE: "_Reference | None" = None
_DISTANCE_CACHE: dict[str, int] = {}

#: ASCII -> amino-acid index, with everything unknown folded into one extra slot.
_CODE_TABLE = np.full(256, len(_ALPHABET), dtype=np.uint8)
for _i, _letter in enumerate(_ALPHABET):
    _CODE_TABLE[ord(_letter)] = _i


class _Reference:
    """abguard's bundled CDRH3 set, in a form numpy can screen in one pass.

    ``counts_by_letter`` is the composition matrix transposed, so the lower-bound
    screen walks 21 contiguous columns instead of materialising a 610k x 21
    difference. ``by_length`` holds one dense uint8 matrix per distinct loop length,
    built once, so aligning a bucket is a fancy-index instead of a Python loop over
    half a million strings.
    """

    __slots__ = ("lengths", "counts_by_letter", "by_length", "slot", "path", "n")

    def __init__(self, lengths: np.ndarray, counts: np.ndarray, codes: np.ndarray,
                 offsets: np.ndarray, path: str):
        self.lengths = lengths.astype(np.int32, copy=False)
        self.counts_by_letter = np.ascontiguousarray(counts.T)
        self.path = path
        self.n = int(lengths.size)
        self.by_length = {}
        self.slot = np.zeros(self.n, dtype=np.int64)
        order = np.argsort(self.lengths, kind="stable")
        boundaries = np.flatnonzero(np.diff(self.lengths[order])) + 1
        for group in np.split(order, boundaries):
            if group.size == 0:
                continue
            length = int(self.lengths[group[0]])
            if length:
                window = offsets[group][:, None] + np.arange(length, dtype=np.int64)
                matrix = codes[window]
            else:
                matrix = np.empty((group.size, 0), dtype=np.uint8)
            self.by_length[length] = (group, matrix)
            self.slot[group] = np.arange(group.size, dtype=np.int64)


def _bundled_reference_path(abguard_py: str | None = None) -> str:
    """Ask abguard itself where its bundled reference lives."""
    py = abguard_py or abguard_python()
    code = (
        "from abguard.reference import bundled_reference, read_manifest;"
        "b = bundled_reference();"
        "print(str(b / read_manifest(b).index_table) + '.parquet')"
    )
    out = subprocess.run([py, "-c", code], capture_output=True, text=True, timeout=300)
    if out.returncode != 0:
        raise RuntimeError(f"could not locate abguard's bundled reference:\n{out.stderr}")
    return out.stdout.strip().splitlines()[-1]


def _dump_reference(abguard_py: str, destination: Path) -> str:
    """Write the bundled reference's CDRH3 column out as one loop per line."""
    parquet = _bundled_reference_path(abguard_py)
    code = (
        "import sys, polars as pl;"
        "df = pl.read_parquet(sys.argv[1]).select('cdrh3').drop_nulls().unique();"
        "open(sys.argv[2], 'w').write('\\n'.join(df['cdrh3'].to_list()))"
    )
    out = subprocess.run(
        [abguard_py, "-c", code, parquet, str(destination)],
        capture_output=True,
        text=True,
        timeout=900,
    )
    if out.returncode != 0:
        raise RuntimeError(f"could not export the CDRH3 reference:\n{out.stderr}")
    return parquet


def load_reference(abguard_py: str | None = None, *, refresh: bool = False) -> _Reference:
    """The bundled CDRH3 reference, memoised in-process and on disk.

    The npz cache exists because the only interpreter that can read abguard's
    Parquet bundle is the one that has polars, and this module runs under the one
    that has pandas. It is derived data: delete it and the next call rebuilds it.
    """
    global _REFERENCE
    if _REFERENCE is not None and not refresh:
        return _REFERENCE
    npz = cache_dir() / f"reference_cdrh3_v{_CACHE_VERSION}.npz"
    if refresh or not npz.is_file():
        npz.parent.mkdir(parents=True, exist_ok=True)
        py = abguard_py or abguard_python()
        with tempfile.TemporaryDirectory() as tmp:
            raw = Path(tmp) / "cdrh3.txt"
            parquet = _dump_reference(py, raw)
            strings = raw.read_text().split("\n")
        strings = [s for s in strings if s]
        lengths, counts, codes, offsets = _encode(strings)
        np.savez_compressed(
            npz,
            codes=codes,
            offsets=offsets,
            lengths=lengths,
            counts=counts,
            source=np.array([parquet]),
        )
    with np.load(npz, allow_pickle=False) as blob:
        _REFERENCE = _Reference(
            blob["lengths"],
            blob["counts"],
            blob["codes"],
            blob["offsets"],
            str(blob["source"][0]),
        )
    _DISTANCE_CACHE.clear()
    return _REFERENCE


def _encode(strings: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Lengths, composition matrix, flat codes and offsets for `strings`."""
    lengths = np.fromiter((len(s) for s in strings), dtype=np.int32, count=len(strings))
    offsets = np.zeros(len(strings) + 1, dtype=np.int64)
    np.cumsum(lengths, out=offsets[1:])
    codes = _CODE_TABLE[np.frombuffer("".join(strings).encode(), dtype=np.uint8)]
    counts = np.zeros((len(strings), len(_ALPHABET) + 1), dtype=np.uint8)
    rows = np.repeat(np.arange(len(strings), dtype=np.int64), lengths)
    np.add.at(counts, (rows, codes), 1)
    return lengths, counts, codes, offsets


def _code_of(sequence: str) -> np.ndarray:
    return _CODE_TABLE[np.frombuffer(sequence.encode(), dtype=np.uint8)]


def _batch_levenshtein(
    query: np.ndarray, reference: _Reference, index: np.ndarray, best: int
) -> int:
    """Smallest Levenshtein between `query` and the reference rows in `index`.

    One DP per distinct reference length, with the whole bucket advanced in
    lockstep: the recurrence is sequential in i and j but embarrassingly parallel
    across entries, so each cell is a single numpy op over the bucket. A row's
    minimum is non-decreasing in i, so once the whole bucket's best partial score
    reaches `best` the bucket cannot improve on it and the DP is abandoned.
    """
    lengths = reference.lengths[index]
    for raw in np.unique(lengths):
        n = int(raw)
        rows = index[lengths == n]
        if n == 0:
            best = min(best, len(query))
            continue
        matrix = reference.by_length[n][1][reference.slot[rows]]
        prev = np.tile(np.arange(n + 1, dtype=np.int16), (rows.size, 1))
        cur = np.empty_like(prev)
        for i in range(1, len(query) + 1):
            cur[:, 0] = i
            differs = matrix != query[i - 1]
            for j in range(1, n + 1):
                np.minimum(prev[:, j] + 1, cur[:, j - 1] + 1, out=cur[:, j])
                np.minimum(cur[:, j], prev[:, j - 1] + differs[:, j - 1], out=cur[:, j])
            prev, cur = cur, prev
            if int(prev.min()) >= best:
                break
        else:
            best = min(best, int(prev[:, n].min()))
    return best


def min_distance_to_reference(cdrh3: str, reference: _Reference) -> int:
    """Exact minimum Levenshtein from `cdrh3` to every reference CDRH3.

    610k exact DPs per candidate would be far too slow to run on every design, so
    a composition lower bound is screened first: an edit changes the amino-acid
    multiset by at most one, hence ``lev >= max(excess_a, excess_b)``, which is
    also ``>= |len(a) - len(b)|``. Buckets of equal bound are then evaluated in
    ascending order and the scan stops as soon as the bound exceeds the best exact
    distance found, so only the genuinely near entries are ever aligned. The answer
    is the same as the brute-force minimum, just reached in milliseconds.
    """
    if not cdrh3:
        return -1
    hit = _DISTANCE_CACHE.get(cdrh3)
    if hit is not None:
        return hit
    codes = _code_of(cdrh3)
    composition = np.zeros(len(_ALPHABET) + 1, dtype=np.int16)
    np.add.at(composition, codes, 1)
    excess_reference = np.zeros(reference.n, dtype=np.int16)
    excess_query = np.zeros(reference.n, dtype=np.int16)
    for letter in range(reference.counts_by_letter.shape[0]):
        difference = reference.counts_by_letter[letter].astype(np.int16)
        difference -= composition[letter]
        excess_reference += np.maximum(difference, 0)
        np.negative(difference, out=difference)
        excess_query += np.maximum(difference, 0)
    bound = np.maximum(excess_reference, excess_query)
    np.maximum(
        bound, np.abs(reference.lengths - len(cdrh3)).astype(np.int16), out=bound
    )
    order = np.argsort(bound, kind="stable")
    sorted_bound = bound[order]
    best = 1 << 30
    start = 0
    while start < order.size:
        level = int(sorted_bound[start])
        if level >= best:
            break
        stop = int(np.searchsorted(sorted_bound, level, side="right"))
        best = _batch_levenshtein(codes, reference, order[start:stop], best)
        start = stop
    _DISTANCE_CACHE[cdrh3] = int(best)
    return int(best)


# --------------------------------------------------------------------------- #
# the official gate
# --------------------------------------------------------------------------- #
#: Columns check_denovo.py writes, in order.
OFFICIAL_COLUMNS = (
    "sequence_id",
    "kind",
    "track_assigned",
    "identity_to_trastuzumab",
    "cdrh3",
    "cdrh3_min_distance",
    "cdrh3_identity_to_trastuzumab",
    "passes",
    "failures",
)


def run_official(
    designs: Path | str,
    *,
    reference_db: Path | str | None = None,
    abguard_py: str | None = None,
    checker: Path | str | None = None,
) -> tuple[pd.DataFrame, int, str]:
    """Run check_denovo.py on `designs` and return (report, exit code, stdout).

    ``--strict-exit`` is always passed: the script writes its report before it
    returns, so a non-zero exit is information rather than a lost run.
    """
    script = Path(checker or OFFICIAL_CHECKER)
    if not script.is_file():
        raise FileNotFoundError(f"official checker not found at {script}")
    py = abguard_py or abguard_python()
    with tempfile.TemporaryDirectory() as tmp:
        report = Path(tmp) / "official_report.csv"
        command = [py, str(script), str(designs), "-o", str(report), "--strict-exit"]
        if reference_db is not None:
            command += ["--reference-db", str(reference_db)]
        done = subprocess.run(command, capture_output=True, text=True, timeout=3600)
        if not report.is_file():
            raise RuntimeError(
                f"official checker wrote no report (exit {done.returncode}):\n"
                f"{done.stdout}\n{done.stderr}"
            )
        frame = pd.read_csv(report, dtype=str, keep_default_na=False)
    return frame, done.returncode, done.stdout


def _as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def _as_float(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.replace("", None), errors="coerce")


# --------------------------------------------------------------------------- #
# public entry point
# --------------------------------------------------------------------------- #
def run_gate(
    designs: Path | str | pd.DataFrame,
    *,
    output: Path | str | None = None,
    reference_db: Path | str | None = None,
    abguard_py: str | None = None,
    checker: Path | str | None = None,
    identity_max: float = INTERNAL_IDENTITY_MAX,
    cdrh3_distance_min: int = INTERNAL_CDRH3_DIST_MIN,
) -> pd.DataFrame:
    """Gate `designs` and return one enriched row per design.

    Args:
        designs: a submission CSV/TSV/Parquet path, or a DataFrame in that shape.
        output: write the enriched report here as CSV (suffix ``.parquet``/``.json``
            also honoured).
        reference_db: an alternative CDRH3 reference for the *official* call.
            ``cdrh3_dist_min_reference`` always measures against abguard's bundled
            bundle, so pass this only when deliberately overriding the official run.
        abguard_py: interpreter that has abguard. Autodetected when omitted.
        checker: path to check_denovo.py. Defaults to the repo's copy.
        identity_max: our internal identity ceiling (official cap is 0.70).
        cdrh3_distance_min: our internal CDRH3 edit floor (official floor is 2).

    Returns:
        A DataFrame with the official verdict verbatim plus our margin columns. The
        official verdict is authoritative; ``internal_pass`` is a strict subset of it.
    """
    temporary: tempfile.TemporaryDirectory | None = None
    if isinstance(designs, pd.DataFrame):
        temporary = tempfile.TemporaryDirectory()
        path = Path(temporary.name) / "designs.csv"
        designs.to_csv(path, index=False)
    else:
        path = Path(designs)
        if not path.is_file():
            raise FileNotFoundError(path)
    try:
        py = abguard_py or abguard_python()
        official, exit_code, stdout = run_official(
            path, reference_db=reference_db, abguard_py=py, checker=checker
        )
        submitted = pd.read_csv(path, dtype=str, keep_default_na=False)
    finally:
        if temporary is not None:
            temporary.cleanup()

    missing = [c for c in OFFICIAL_COLUMNS if c not in official.columns]
    if missing:
        raise RuntimeError(f"official report is missing {missing}; got {list(official.columns)}")

    reference = load_reference(py)
    cdrh3 = official["cdrh3"].fillna("")
    identity = _as_float(official["identity_to_trastuzumab"])

    unique = sorted({c for c in cdrh3 if c})
    to_trastuzumab = {c: levenshtein(c, TRASTUZUMAB_CDRH3) for c in unique}
    to_reference = {c: min_distance_to_reference(c, reference) for c in unique}

    out = pd.DataFrame(
        {
            "sequence_id": official["sequence_id"],
            "kind": official["kind"],
            "official_pass": _as_bool(official["passes"]),
            "official_track_assigned": official["track_assigned"],
            "official_reason": official["failures"],
            "identity_to_trastuzumab": identity,
            "margin_vs_identity_cap": IDENTITY_CAP - identity,
            "cdrh3": cdrh3,
            "cdrh3_len": cdrh3.str.len().astype("Int64"),
            "cdrh3_dist_trastuzumab": pd.array(
                [to_trastuzumab.get(c) for c in cdrh3], dtype="Int64"
            ),
            "cdrh3_dist_min_reference": pd.array(
                [to_reference.get(c) for c in cdrh3], dtype="Int64"
            ),
            "official_cdrh3_min_distance": pd.to_numeric(
                official["cdrh3_min_distance"].replace("", None), errors="coerce"
            ).astype("Int64"),
        }
    )

    reasons: list[str] = []
    passes: list[bool] = []
    for row in out.itertuples(index=False):
        why: list[str] = []
        if not row.official_pass:
            why.append("official gate failed")
        if pd.isna(row.identity_to_trastuzumab) or row.identity_to_trastuzumab > identity_max:
            why.append(f"identity {row.identity_to_trastuzumab} > internal max {identity_max}")
        if pd.isna(row.cdrh3_dist_trastuzumab) or row.cdrh3_dist_trastuzumab < cdrh3_distance_min:
            why.append(
                f"CDRH3 only {row.cdrh3_dist_trastuzumab} edit(s) from trastuzumab; "
                f"internal min {cdrh3_distance_min}"
            )
        if (
            pd.isna(row.cdrh3_dist_min_reference)
            or row.cdrh3_dist_min_reference < cdrh3_distance_min
        ):
            why.append(
                f"CDRH3 only {row.cdrh3_dist_min_reference} edit(s) from the reference set; "
                f"internal min {cdrh3_distance_min}"
            )
        reasons.append("; ".join(why))
        passes.append(not why)
    out["internal_pass"] = passes
    out["internal_reason"] = reasons

    descriptions = dict(
        zip(submitted.get("sequence_id", pd.Series(dtype=str)),
            submitted.get("design_description", pd.Series(dtype=str)))
    )
    out["design_description"] = [descriptions.get(i, "") for i in out["sequence_id"]]

    out.attrs["official_exit_code"] = exit_code
    out.attrs["official_stdout"] = stdout
    out.attrs["reference_path"] = reference.path
    out.attrs["reference_entries"] = reference.n
    out.attrs["abguard_python"] = py

    if output is not None:
        destination = Path(output)
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.suffix == ".parquet":
            out.to_parquet(destination, index=False)
        elif destination.suffix == ".json":
            out.to_json(destination, orient="records", indent=2)
        else:
            out.to_csv(destination, index=False)
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run the official de novo gate and add our stricter margin columns."
    )
    parser.add_argument("designs", help="submission table (sequence_id,kind,heavy_sequence,...)")
    parser.add_argument("-o", "--output", default=None, help="write the enriched report here")
    parser.add_argument("--reference-db", default=None, help="override the OFFICIAL reference")
    parser.add_argument("--abguard-py", default=None, help="interpreter that has abguard")
    parser.add_argument("--checker", default=None, help="path to check_denovo.py")
    parser.add_argument(
        "--identity-max", type=float, default=INTERNAL_IDENTITY_MAX, help="internal identity cap"
    )
    parser.add_argument(
        "--cdrh3-dist-min",
        type=int,
        default=INTERNAL_CDRH3_DIST_MIN,
        help="internal CDRH3 edit floor",
    )
    parser.add_argument(
        "--refresh-reference", action="store_true", help="rebuild the local reference cache"
    )
    parser.add_argument(
        "--strict-exit",
        action="store_true",
        help="exit 1 unless every design passes BOTH the official gate and our margin",
    )
    args = parser.parse_args(argv)

    if args.refresh_reference:
        load_reference(args.abguard_py, refresh=True)

    frame = run_gate(
        args.designs,
        output=args.output,
        reference_db=args.reference_db,
        abguard_py=args.abguard_py,
        checker=args.checker,
        identity_max=args.identity_max,
        cdrh3_distance_min=args.cdrh3_dist_min,
    )
    width = max((len(str(i)) for i in frame["sequence_id"]), default=12)
    for row in frame.itertuples(index=False):
        official = "PASS" if row.official_pass else "FAIL"
        internal = "PASS" if row.internal_pass else "FAIL"
        print(
            f"official={official} internal={internal}  {str(row.sequence_id):<{width}}  "
            f"id={row.identity_to_trastuzumab:.4f} margin={row.margin_vs_identity_cap:+.4f} "
            f"cdrh3={row.cdrh3 or '-'}({row.cdrh3_len}) "
            f"d_tz={row.cdrh3_dist_trastuzumab} d_ref={row.cdrh3_dist_min_reference}"
        )
        if not row.internal_pass:
            print(f"    <- {row.internal_reason}")
    n_official = int(frame["official_pass"].sum())
    n_internal = int(frame["internal_pass"].sum())
    print(
        f"\nofficial de novo: {n_official}/{len(frame)} eligible "
        f"(exit {frame.attrs['official_exit_code']})"
    )
    print(f"internal margin:  {n_internal}/{len(frame)} eligible")
    print(
        f"reference: {frame.attrs['reference_entries']} CDRH3 entries from "
        f"{frame.attrs['reference_path']}"
    )
    if args.output:
        print(f"wrote {args.output}")
    return 1 if (args.strict_exit and n_internal < len(frame)) else 0


if __name__ == "__main__":
    raise SystemExit(main())
