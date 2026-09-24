"""
LEVEL 2 - Sequence alignment from scratch (dynamic programming).

    python -m pytest level2_alignment/test_level2.py -q

You implement global (Needleman-Wunsch) and local (Smith-Waterman) alignment
with a real substitution matrix and affine gap penalties. This is the first
algorithm in the ladder where a subtle bug still produces a plausible-looking
answer, so the tests check the score AND the traceback consistency.

Conventions used throughout:
  * gap character is '-'
  * affine gap cost for a run of length k is: gap_open + gap_extend * k
    (both given as POSITIVE numbers; they are subtracted from the score)
  * ties in the traceback are broken diagonal > up > left, so results are
    reproducible
"""

from blosum62 import BLOSUM62  # noqa: F401  (dict of (a, b) -> int)

GAP = "-"


def score_pair(a, b, matrix=BLOSUM62):
    """Substitution score for two residues. Unknown residues score 0."""
    raise NotImplementedError("TODO")


def score_alignment(aln_a, aln_b, matrix=BLOSUM62, gap_open=11, gap_extend=1):
    """Score an already-aligned pair of strings with affine gaps.

    This is your referee: any alignment your DP produces must score exactly
    what the DP said it did. Gap-vs-gap columns should not occur; raise
    ValueError if you see one.
    """
    raise NotImplementedError("TODO")


def needleman_wunsch(a, b, matrix=BLOSUM62, gap_open=11, gap_extend=1):
    """Global alignment. Return (score, aligned_a, aligned_b).

    Use three DP matrices (M, Ix, Iy) - the Gotoh formulation - so that an
    affine gap is charged open+extend once, not open per residue.
    """
    raise NotImplementedError("TODO")


def smith_waterman(a, b, matrix=BLOSUM62, gap_open=11, gap_extend=1):
    """Local alignment. Return (score, aligned_a, aligned_b, start_a, start_b)
    where start_a/start_b are 0-based starts of the aligned region in the
    ORIGINAL sequences.

    Differences from global alignment, all of which matter:
      * cells are floored at 0
      * the traceback starts at the maximum cell, not the corner
      * the traceback stops when it reaches a 0
    """
    raise NotImplementedError("TODO")


def percent_identity(aln_a, aln_b):
    """Identical columns / aligned columns (columns with a gap are counted in
    the denominator). Return 0.0 for an empty alignment."""
    raise NotImplementedError("TODO")
