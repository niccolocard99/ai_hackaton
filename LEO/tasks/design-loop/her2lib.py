"""Shared constants and helpers for the HER2 optimisation-track design loop.

Everything here is about one thing: building a full trastuzumab-framework
VH/VL pair from a designed 10-residue CDRH3 loop, and back again.

The 10-mer is what the SPR-768 table varies, so it is the only region with
measured affinity labels. In the VH it sits between the `...YYCSR` anchor and
the `YWGQGT` anchor:

    ...TAVYYCSR | WGGDGFYAMD | YWGQGTLVTVSS
                  ^^^^^^^^^^
                  positions 99-108, the designed window

abguard's own CDRH3 definition is wider than this window (it includes the SR
prefix and the Y suffix). Which one the identity gate uses decides how many
substitutions a legal design needs, so it is measured in probe_abguard.py
rather than assumed.
"""
from __future__ import annotations

AA = "ACDEFGHIKLMNPQRSTVWY"

# Trastuzumab, exactly as the SPR-768 table carries it.
TZ_VH = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFT"
         "ISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")
TZ_VL = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGT"
         "DFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")

TZ_LOOP = "WGGDGFYAMD"          # the 10-mer the SPR table varies
PREFIX, SUFFIX = "YYCSR", "YWGQGT"  # anchors around it in the VH
TZ_KD_NM = 0.152                 # trastuzumab's own measured KD, the bar to beat

# Position labels for the 10-mer, Kabat-style, as the ledger uses them.
POSITIONS = [99, 100, 101, 102, 103, 104, 105, 106, 107, 108]
assert len(POSITIONS) == len(TZ_LOOP)

# EXP-002, cross-checked over 8 trastuzumab-HER2 structures: in how many of the 8
# does this CDRH3 position contact HER2 (4.5 A heavy-atom cutoff). Copied from
# LEO/tasks/litreview-optimisation/consensus_4.5.txt, "CDRH3 POSITION CONSENSUS"
# table, column n_X/n_informative (all 8 structures informative for every one of
# these positions). 0 means it never touches the antigen in any structure, so it
# is the cheapest place to spend a substitution.
#
# NOTE: an earlier draft of this file assumed 8 for positions 99-101 and 103-105.
# That was a guess, not a measurement, and it was wrong (G100 is 3, G101 is 2,
# F104 is 4). Any output produced before this table was corrected is invalid.
CONTACT_FREQ = {99: 5, 100: 3, 101: 2, 102: 6, 103: 7, 104: 4, 105: 7,
                106: 0, 107: 0, 108: 0}
N_STRUCT = 8
FREE_POSITIONS = [p for p, n in CONTACT_FREQ.items() if n == 0]  # 106, 107, 108
# Position 109 (Y, the last residue of the 11-mer) contacts in 1 of 8 but lies
# outside the 10-mer the SPR table varies, so no affinity label can speak to it.

_i = TZ_VH.index(PREFIX) + len(PREFIX)
LOOP_SLICE = slice(_i, _i + len(TZ_LOOP))
assert TZ_VH[LOOP_SLICE] == TZ_LOOP, "loop anchors do not locate the 10-mer"


def vh_with_loop(loop: str) -> str:
    """Trastuzumab VH carrying `loop` in place of the wild-type 10-mer."""
    if len(loop) != len(TZ_LOOP):
        raise ValueError(f"loop must be {len(TZ_LOOP)} aa, got {len(loop)}")
    return TZ_VH[: LOOP_SLICE.start] + loop + TZ_VH[LOOP_SLICE.stop :]


def nmut(loop: str, ref: str = TZ_LOOP) -> int:
    """Hamming distance between two equal-length loops."""
    return sum(a != b for a, b in zip(loop, ref))


def subs(loop: str, ref: str = TZ_LOOP) -> list[str]:
    """Human-readable substitution list, e.g. ['D102S', 'A106L']."""
    return [f"{r}{p}{a}" for r, a, p in zip(ref, loop, POSITIONS) if r != a]


def epitope_risk(loop: str, ref: str = TZ_LOOP) -> float:
    """Expected number of HER2-contacting CDRH3 positions a design disturbs.

    Sum over mutated positions of (structures contacting / 8): each term is the
    chance, over the 8 reference structures, that the mutated position was a
    contact. 0.0 means every substitution sits on a position that never touches
    HER2 in any structure. This is geometry, not energy (a rare contact can
    still be a hot spot), so it ranks designs by how much contact they spend; it
    does not say any single substitution breaks the epitope.

    Floor for a legal design: it needs >= 4 substitutions and only 3 positions
    are free, so the minimum is the 3 free ones plus the least-contacting fourth,
    G101 at 2/8 -> 0.25.
    """
    return sum(CONTACT_FREQ[p] for r, a, p in zip(ref, loop, POSITIONS) if r != a) / N_STRUCT


def design_rows(loops, prefix="LEO-OPT"):
    """Design-table rows in the shape the track checkers want."""
    for k, loop in enumerate(loops, 1):
        yield {
            "sequence_id": f"{prefix}-{k:04d}",
            "kind": "paired",
            "heavy_sequence": vh_with_loop(loop),
            "light_sequence": TZ_VL,
            "design_description": f"CDRH3 {loop} ({'/'.join(subs(loop)) or 'wild-type'})",
        }
