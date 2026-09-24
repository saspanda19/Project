"""LEVEL 1 reference solution. Read it AFTER you have your own working version
-- comparing two working implementations teaches more than copying one.

Run the demo (prints a report, writes figures/level1_gc.png):

    python level1_sequences/solution.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from starter import CODON_TABLE  # noqa: E402

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def read_fasta(path):
    records, name, desc, chunks = [], None, "", []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, desc, "".join(chunks)))
                header = line[1:].strip().split(None, 1)
                name = header[0]
                desc = header[1] if len(header) > 1 else ""
                chunks = []
            else:
                chunks.append(line.upper())
    if name is not None:
        records.append((name, desc, "".join(chunks)))
    return records


def gc_content(seq):
    if not seq:
        return 0.0
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq)


def reverse_complement(seq):
    return "".join(COMPLEMENT.get(b, "N") for b in seq.upper()[::-1])


def translate(seq, stop_at_stop=False):
    seq = seq.upper()
    aas = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        aa = CODON_TABLE.get(seq[i:i + 3], "X")
        if aa == "*" and stop_at_stop:
            break
        aas.append(aa)
    return "".join(aas)


def find_orfs(seq, min_aa=30):
    orfs = []
    for strand, s in (("+", seq.upper()), ("-", reverse_complement(seq))):
        for frame in range(3):
            i = frame
            while i + 3 <= len(s):
                if s[i:i + 3] == "ATG":
                    j = i
                    while j + 3 <= len(s):
                        if CODON_TABLE.get(s[j:j + 3], "X") == "*":
                            protein = translate(s[i:j])
                            if len(protein) >= min_aa:
                                orfs.append({"strand": strand, "frame": frame,
                                             "start": i, "end": j + 3,
                                             "protein": protein})
                            i = j  # continue scanning after this ORF
                            break
                        j += 3
                    else:
                        break  # ran off the end with no stop: not an ORF
                i += 3
    orfs.sort(key=lambda o: -len(o["protein"]))
    return orfs


def sliding_gc(seq, window=50, step=10):
    out = []
    for start in range(0, len(seq) - window + 1, step):
        out.append((start + window // 2, gc_content(seq[start:start + window])))
    return out


# ------------------------------------------------------------------ demo ---
def main():
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    records = read_fasta(os.path.join(root, "data", "dna_sample.fasta"))

    print(f"{'name':<10}{'len':>6}{'GC':>8}{'ORFs':>6}  longest protein")
    print("-" * 64)
    for name, _desc, seq in records:
        orfs = find_orfs(seq, min_aa=30)
        longest = orfs[0]["protein"] if orfs else ""
        preview = (longest[:24] + "...") if len(longest) > 24 else longest
        print(f"{name:<10}{len(seq):>6}{gc_content(seq):>8.3f}"
              f"{len(orfs):>6}  {len(longest):>3}aa {preview}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        print("\n(matplotlib unavailable - skipping figure)")
        return
    figdir = os.path.join(root, "figures")
    os.makedirs(figdir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4))
    for name, _d, seq in records:
        pts = sliding_gc(seq, window=40, step=5)
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], label=name, lw=1.4)
    ax.set_xlabel("position (bp)")
    ax.set_ylabel("GC fraction (40 bp window)")
    ax.set_title("Level 1 - sliding-window GC content")
    ax.legend(fontsize=7, ncol=3)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out = os.path.join(figdir, "level1_gc.png")
    fig.savefig(out, dpi=150)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
