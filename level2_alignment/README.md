# Level 2 — Alignment from scratch

**Time:** 1–2 days · **Prereqs:** Level 1, comfort with 2-D arrays · **New ideas:** dynamic programming, affine gaps, global vs local, tie-breaking

## Why this level exists

Alignment is the first algorithm in the ladder where a subtle bug still gives a *plausible* answer. A DP with the wrong gap bookkeeping returns a number and a pretty alignment, and nothing screams. The defence is a **referee**: `score_alignment` re-scores the alignment your DP produced, independently. If the two numbers disagree, you have a bug — and the tests enforce this on every pair.

That habit — an independent check that doesn't share code with the thing it's checking — is the single most transferable thing in this repo.

## What you build

`score_pair`, `score_alignment`, `needleman_wunsch` (global), `smith_waterman` (local), `percent_identity`, with BLOSUM62 and affine gaps (`open=11`, `extend=1`, the BLAST defaults).

Use the **Gotoh** three-matrix formulation: `M` (column ends in a match), `Ix` (gap in b), `Iy` (gap in a). A single-matrix DP cannot charge affine gaps correctly.

```bash
python -m pytest level2_alignment/test_level2.py -q
python level2_alignment/solution.py
```

## Checkpoints

1. `needleman_wunsch(s, s)` returns the sum of the diagonal of BLOSUM62 for `s` and no gaps.
2. Every alignment you produce re-scores to exactly the DP score.
3. Removing a contiguous block of residues yields **one** gap run, not several — that's the affine penalty doing its job.
4. Smith–Waterman on the `dom_host1/dom_host2` pair recovers the 35-residue shared domain; Needleman–Wunsch on the same pair smears it. Run the demo and look at both.

## Questions to answer in your notes

- Why can't one DP matrix handle affine gaps? Write the recurrence that fails.
- Your traceback breaks ties diagonal > up > left. Pick a pair where a different rule gives a *different alignment with the same score*, and say why reproducibility matters when the alignment feeds a downstream model.
- Global alignment gave 38.7% identity on the domain pair, local gave ~95%. If a paper reports "X% identity," what do you now have to ask before believing it?
- Complexity is O(nm) time and memory. At what sequence length does the memory become the problem, and what does Hirschberg's algorithm trade to fix it?

## Extensions

1. **Linear-space alignment** (Hirschberg). Same score, O(min(n,m)) memory. The most satisfying rewrite in this repo.
2. **Multiple sequence alignment**, progressive: align the closest pair, then align the profile to the next sequence. Use the four `prot_*` sequences.
3. **Statistical significance**: shuffle one sequence 1000× and compare the real local score to the null distribution. You've just derived the idea behind a BLAST E-value.
4. **Benchmark against Biopython** (`pip install biopython`, `Bio.Align.PairwiseAligner`). Same score with the same parameters? If not, find out which of you is wrong — it isn't always them, but it usually is you.

## Graduating

Tests pass, your alignments always re-score to the DP score, and you can explain out loud when local beats global.
