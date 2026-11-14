"""The S97 hedge (LEO handoff item 1): +1 novelty margin from a position nothing in the reference set
ever varies.

abguard's CDRH3 is the 13-mer 97-109 (SRWGGDGFYAMDY) and the mutable window is its 10-mer core
(99-108). The 525k reference loops NEVER vary position 97, so any S97X substitution is one Levenshtein
edit further from every known variant — measured by the JapanFold-RFD4 agent: an S97 substitution
lifted 23 of 24 top parents to novelty margin >= 4 while staying legal. The effect of S97X on binding
is UNMEASURED (it makes no heavy-atom contact to HER2 in 1N8Z, which is the argument for trying it);
the pre-filter cannot see position 97 (its one-hot covers the 10-mer), so children inherit their
parent's window scores — the hedge buys margin, not predicted affinity.

Each round this generator takes the previous round's selected parents (or a cold-start table) and
emits one child per parent per residue, with the parent's heavy chain edited only at position 97
(two residues before the window start). `pin` optionally also fixes window positions in the PARENT
selection (item 7: keep Y105 fixed — the largest buried surface and an alanine-scan hotspot).
"""
from __future__ import annotations

import polars as pl

from abpipe.base import Context, Generator
from abpipe.generators.scored_evolution import file_sha256, read_table
from abpipe.registry import register

AA = "ACDEFGHIKLMNPQRSTVWY"


@register("generator", "s97_hedge")
class S97Hedge(Generator):
    """+1 CDRH3 novelty distance from every reference loop, via the never-varied S97."""

    version = "1"

    def __init__(self, parents: str | None = None, parents_sha256: str | None = None,
                 residues: str = "WY", pin: tuple = (), pin_to: str = "Y"):
        for r in residues:
            if r not in AA:
                raise ValueError(f"residue {r!r} is not a standard amino acid")
        if pin and (not all(isinstance(p, int) and 0 <= p < 10 for p in pin) or pin_to not in AA):
            raise ValueError("pin must be window indices 0-9 and pin_to a standard residue")
        self.parents, self.parents_sha256 = parents, parents_sha256
        self.residues, self.pin, self.pin_to = residues, tuple(pin), pin_to

    def generate(self, n: int, seed: int, ctx: Context) -> pl.DataFrame:
        parent = ctx.require_parent()
        if not self.parents:
            raise ValueError("s97_hedge needs `parents` (a scored/ranked table with heavy_sequence)")
        if self.parents_sha256 and not file_sha256(self.parents).startswith(self.parents_sha256.lower()):
            raise ValueError(f"{self.parents} changed since the config was written (sha256 mismatch)")
        src = read_table(self.parents)
        if "rank" in src.columns:
            src = src.sort("rank")
        rows = []
        s97 = parent.window_start - 2
        for r in src.iter_rows(named=True):
            heavy = r["heavy_sequence"]
            if heavy is None or len(heavy) != len(parent.heavy) or heavy[s97] != parent.heavy[s97]:
                continue                                    # not the parent scaffold; the hedge does not apply
            if self.pin and any(heavy[parent.window_start + p] != self.pin_to for p in self.pin):
                continue                                    # parent breaks the pin (e.g. keep Y105 fixed)
            for res in self.residues:
                rows.append({
                    "sequence_id": f"{self.name}_s{seed}_{len(rows):06d}", "kind": r.get("kind") or "paired",
                    "heavy_sequence": heavy[:s97] + res + heavy[s97 + 1:],
                    "light_sequence": r.get("light_sequence") or parent.light,
                    "design_description": f"S97{res} hedge on {r['sequence_id']}; +1 CDRH3 distance "
                                          f"from every reference loop, window unchanged",
                    "parent_id": r["sequence_id"], "origin": f"{self.name}:S97{res}", "seed": seed})
            if len(rows) >= n:
                break
        if not rows:
            raise ValueError("s97_hedge produced nothing: every parent differs from the parent scaffold "
                             "at position 97, or breaks the pin")
        return pl.DataFrame(rows, schema={"sequence_id": pl.String, "kind": pl.String, "heavy_sequence": pl.String,
                                          "light_sequence": pl.String, "design_description": pl.String,
                                          "parent_id": pl.String, "origin": pl.String, "seed": pl.Int64})
