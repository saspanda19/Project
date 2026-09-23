"""LEVEL 2 reference solution - Gotoh affine-gap alignment.

    python level2_alignment/solution.py

Writes figures/level2_dotplot.png and prints an identity matrix over the
synthetic protein family in data/proteins.fasta.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from blosum62 import BLOSUM62  # noqa: E402

GAP = "-"
NEG = float("-inf")


def score_pair(a, b, matrix=BLOSUM62):
    return matrix.get((a, b), matrix.get((b, a), 0))


def score_alignment(aln_a, aln_b, matrix=BLOSUM62, gap_open=11, gap_extend=1):
    if len(aln_a) != len(aln_b):
        raise ValueError("aligned strings must be the same length")
    total, in_gap = 0, False
    for ca, cb in zip(aln_a, aln_b):
        if ca == GAP and cb == GAP:
            raise ValueError("gap aligned to gap")
        if ca == GAP or cb == GAP:
            total -= gap_extend + (0 if in_gap else gap_open)
            in_gap = True
        else:
            total += score_pair(ca, cb, matrix)
            in_gap = False
    return total


def _gotoh(a, b, matrix, gap_open, gap_extend, local):
    """Shared DP core. M[i][j] ends with a match/mismatch column,
    Ix ends with a gap in b (consumes a), Iy ends with a gap in a."""
    n, m = len(a), len(b)
    M = [[NEG] * (m + 1) for _ in range(n + 1)]
    Ix = [[NEG] * (m + 1) for _ in range(n + 1)]
    Iy = [[NEG] * (m + 1) for _ in range(n + 1)]
    # pointers: 0 diag(M), 1 up(Ix), 2 left(Iy); per matrix
    pM = [[None] * (m + 1) for _ in range(n + 1)]
    pX = [[None] * (m + 1) for _ in range(n + 1)]
    pY = [[None] * (m + 1) for _ in range(n + 1)]

    M[0][0] = 0
    if not local:
        for i in range(1, n + 1):
            Ix[i][0] = -(gap_open + gap_extend * i)
            pX[i][0] = 1 if i > 1 else 0
        for j in range(1, m + 1):
            Iy[0][j] = -(gap_open + gap_extend * j)
            pY[0][j] = 2 if j > 1 else 0

    best, best_cell = (0, (0, 0)) if local else (None, None)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = score_pair(a[i - 1], b[j - 1], matrix)
            cands = ((M[i - 1][j - 1], 0), (Ix[i - 1][j - 1], 1),
                     (Iy[i - 1][j - 1], 2))
            val, ptr = max(cands, key=lambda t: (t[0], -t[1]))
            val = NEG if val == NEG else val + s
            if local and (val == NEG or val < 0):
                val, ptr = 0, None
            M[i][j], pM[i][j] = val, ptr

            cands = ((M[i - 1][j] - gap_open - gap_extend, 0),
                     (Ix[i - 1][j] - gap_extend, 1))
            val, ptr = max(cands, key=lambda t: (t[0], -t[1]))
            Ix[i][j], pX[i][j] = val, ptr

            cands = ((M[i][j - 1] - gap_open - gap_extend, 0),
                     (Iy[i][j - 1] - gap_extend, 2))
            val, ptr = max(cands, key=lambda t: (t[0], -t[1]))
            Iy[i][j], pY[i][j] = val, ptr

            if local and M[i][j] > best:
                best, best_cell = M[i][j], (i, j)
    return M, Ix, Iy, pM, pX, pY, best, best_cell


def _traceback(a, b, M, Ix, Iy, pM, pX, pY, i, j, local):
    out_a, out_b = [], []
    layer = max(((M[i][j], 0), (Ix[i][j], 1), (Iy[i][j], 2)),
                key=lambda t: (t[0], -t[1]))[1] if not local else 0
    while i > 0 or j > 0:
        if local and layer == 0 and M[i][j] == 0:
            break
        if layer == 0:
            if i == 0 or j == 0:
                layer = 1 if i > 0 else 2
                continue
            out_a.append(a[i - 1])
            out_b.append(b[j - 1])
            nxt = pM[i][j]
            i, j = i - 1, j - 1
            layer = 0 if nxt is None else nxt
        elif layer == 1:
            out_a.append(a[i - 1])
            out_b.append(GAP)
            nxt = pX[i][j]
            i -= 1
            layer = 0 if nxt is None else nxt
        else:
            out_a.append(GAP)
            out_b.append(b[j - 1])
            nxt = pY[i][j]
            j -= 1
            layer = 0 if nxt is None else nxt
    return "".join(reversed(out_a)), "".join(reversed(out_b)), i, j


def needleman_wunsch(a, b, matrix=BLOSUM62, gap_open=11, gap_extend=1):
    M, Ix, Iy, pM, pX, pY, _, _ = _gotoh(a, b, matrix, gap_open, gap_extend,
                                         local=False)
    n, m = len(a), len(b)
    score = max(M[n][m], Ix[n][m], Iy[n][m])
    aa, bb, _, _ = _traceback(a, b, M, Ix, Iy, pM, pX, pY, n, m, local=False)
    return score, aa, bb


def smith_waterman(a, b, matrix=BLOSUM62, gap_open=11, gap_extend=1):
    M, Ix, Iy, pM, pX, pY, best, cell = _gotoh(a, b, matrix, gap_open,
                                               gap_extend, local=True)
    i, j = cell
    aa, bb, si, sj = _traceback(a, b, M, Ix, Iy, pM, pX, pY, i, j, local=True)
    return best, aa, bb, si, sj


def percent_identity(aln_a, aln_b):
    if not aln_a:
        return 0.0
    same = sum(1 for x, y in zip(aln_a, aln_b) if x == y and x != GAP)
    return same / len(aln_a)


def pretty(aln_a, aln_b, width=60):
    mid = "".join("|" if x == y and x != GAP else
                  (" " if GAP in (x, y) else
                   ("+" if score_pair(x, y) > 0 else "."))
                  for x, y in zip(aln_a, aln_b))
    out = []
    for k in range(0, len(aln_a), width):
        out += [aln_a[k:k + width], mid[k:k + width], aln_b[k:k + width], ""]
    return "\n".join(out)


# ------------------------------------------------------------------ demo ---
def read_fasta(path):
    """(Level 1 again - 8 lines, so the demo stays self-contained.)"""
    recs, name, desc, chunks = [], None, "", []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name:
                recs.append((name, desc, "".join(chunks)))
            h = line[1:].split(None, 1)
            name, desc, chunks = h[0], (h[1] if len(h) > 1 else ""), []
        elif line:
            chunks.append(line.upper())
    if name:
        recs.append((name, desc, "".join(chunks)))
    return recs


def main():
    recs = read_fasta(os.path.join(ROOT, "data", "proteins.fasta"))
    seqs = {n: s for n, _d, s in recs}

    print("Global alignments against the ancestor (prot_parent):\n")
    parent = seqs["prot_parent"]
    for name in ("prot_close", "prot_mid", "prot_far"):
        sc, aa, bb = needleman_wunsch(parent, seqs[name])
        assert sc == score_alignment(aa, bb), "DP score must match the referee"
        print(f"  {name:<12} score={sc:>5}  identity={percent_identity(aa, bb):.1%}"
              f"  length={len(aa)}")

    print("\nLocal alignment of two unrelated proteins that share one domain:")
    sc, aa, bb, si, sj = smith_waterman(seqs["dom_host1"], seqs["dom_host2"])
    print(f"  score={sc}  starts at {si} / {sj}  length={len(aa)}"
          f"  identity={percent_identity(aa, bb):.1%}")
    print(pretty(aa, bb))
    gsc, gaa, gbb = needleman_wunsch(seqs["dom_host1"], seqs["dom_host2"])
    print(f"  the GLOBAL alignment of the same pair: score={gsc}, "
          f"identity={percent_identity(gaa, gbb):.1%} - "
          "it smears the shared domain across unrelated flanks.")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception:
        return
    x, y = seqs["dom_host1"], seqs["dom_host2"]
    w = 4
    dot = np.zeros((len(x) - w, len(y) - w))
    for i in range(len(x) - w):
        for j in range(len(y) - w):
            dot[i, j] = sum(1 for k in range(w) if x[i + k] == y[j + k])
    figdir = os.path.join(ROOT, "figures")
    os.makedirs(figdir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.imshow(dot.T >= 3, cmap="Greys", origin="lower", interpolation="nearest")
    ax.set_xlabel("dom_host1")
    ax.set_ylabel("dom_host2")
    ax.set_title("Level 2 - dot plot (word size 4)\nthe diagonal streak is the shared domain")
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "level2_dotplot.png"), dpi=150)
    print(f"\nwrote {os.path.join(figdir, 'level2_dotplot.png')}")


if __name__ == "__main__":
    main()
