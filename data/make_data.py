"""
Generate every dataset used by the ladder. Deterministic: same seed -> same files.

Everything here is SYNTHETIC on purpose:
  * it runs offline, in seconds, on a laptop
  * the ground truth is known, so you can tell whether your code is right
  * the planted signal is realistic in *shape* (a motif, a composition bias,
    a batch confound) even though the molecules are made up

Level READMEs tell you how to swap in real data (RCSB, UniProt, PDB) when you
want to graduate from the sandbox.

    python data/make_data.py
"""
import os
import math
import random

HERE = os.path.dirname(os.path.abspath(__file__))
AA = "ACDEFGHIKLMNPQRSTVWY"
HYDROPHOBIC = "AVILMFWY"

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
SENSE_CODONS = [c for c, a in CODON_TABLE.items() if a != "*"]


def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, desc, seq in records:
            fh.write(f">{name} {desc}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


# ---------------------------------------------------------------- level 1 ---
def make_dna(rng):
    """Six 'genes'. Each has a real ORF on the forward strand; gc_target
    controls base composition so GC-content code has something to find."""
    records = []
    specs = [
        ("gene_A", 0.35, 90), ("gene_B", 0.50, 120), ("gene_C", 0.65, 75),
        ("gene_D", 0.42, 150), ("gene_E", 0.58, 60), ("gene_F", 0.70, 105),
    ]
    for name, gc, orf_codons in specs:
        def base():
            return rng.choice("GC") if rng.random() < gc else rng.choice("AT")

        def codon():
            # pick a sense codon whose own GC is close to the target, so the
            # whole record (not just the UTRs) hits gc_target
            cands = [rng.choice(SENSE_CODONS) for _ in range(8)]
            return min(cands, key=lambda c: abs(
                (c.count("G") + c.count("C")) / 3 - gc) + rng.random() * 0.05)

        utr5 = "".join(base() for _ in range(rng.randint(20, 40)))
        body = "".join(codon() for _ in range(orf_codons))
        stop = rng.choice(["TAA", "TAG", "TGA"])
        utr3 = "".join(base() for _ in range(rng.randint(20, 40)))
        seq = utr5 + "ATG" + body + stop + utr3
        desc = f"synthetic|gc_target={gc}|orf_aa={orf_codons + 1}"
        records.append((name, desc, seq))
    write_fasta(os.path.join(HERE, "dna_sample.fasta"), records)
    return records


# ---------------------------------------------------------------- level 2 ---
def make_protein_pairs(rng):
    """Pairs with known evolutionary distance: a parent sequence plus mutated
    children. You know the true alignment, so you can sanity-check yours."""
    records = []
    parent = "".join(rng.choice(AA) for _ in range(120))
    records.append(("prot_parent", "ancestor", parent))
    for name, subs, indels in [("prot_close", 8, 1), ("prot_mid", 30, 3),
                               ("prot_far", 60, 6)]:
        s = list(parent)
        for _ in range(subs):
            i = rng.randrange(len(s))
            s[i] = rng.choice(AA)
        for _ in range(indels):
            i = rng.randrange(len(s))
            if rng.random() < 0.5:
                del s[i:i + rng.randint(1, 4)]
            else:
                s[i:i] = [rng.choice(AA) for _ in range(rng.randint(1, 4))]
        records.append((name, f"subs={subs};indels={indels}", "".join(s)))
    # a local-alignment case: shared domain inside unrelated flanks
    domain = "".join(rng.choice(AA) for _ in range(35))
    for name in ("dom_host1", "dom_host2"):
        left = "".join(rng.choice(AA) for _ in range(rng.randint(25, 45)))
        right = "".join(rng.choice(AA) for _ in range(rng.randint(25, 45)))
        records.append((name, "shares a 35aa domain with the other dom_host",
                        left + domain + right))
    write_fasta(os.path.join(HERE, "proteins.fasta"), records)
    return records


# ---------------------------------------------------------------- level 3 ---
def make_binder_dataset(rng, n=1200):
    """Peptide binder/non-binder labels with THREE deliberate traps:

      1. real signal   : an 'RGD' motif + hydrophobic enrichment -> binder
      2. near-duplicates: binders come in families of 3 mutated siblings, so a
                          random split leaks family members across train/test
      3. batch confound : column 'plate' correlates with label at 0.8

    A beginner who does train_test_split(random_state=0) and reports AUC will
    get a beautiful, wrong number. Level 3 is about finding out why.
    """
    rows = []
    n_fam = n // 6  # each family -> 3 siblings; one positive + one negative fam

    def family(parent, label, tag):
        out = []
        for sib in range(3):
            s = list(parent)
            for _ in range(rng.randint(0, 1)):  # siblings are NEARLY identical
                i = rng.randrange(len(s))
                s[i] = rng.choice(AA)
            out.append({"peptide_id": f"{tag}_sib{sib}", "family": tag,
                        "sequence": "".join(s), "label": label})
        return out

    for f in range(n_fam):
        # positives: the motif is present only ~55% of the time, and the
        # compositional signal is weak -- so the task is learnable but not
        # trivially separable. Real assay data behaves like this.
        core = "".join(rng.choice(AA) for _ in range(rng.randint(10, 17)))
        if rng.random() < 0.55:
            i = rng.randrange(1, len(core) - 1)
            core = core[:i] + "RGD" + core[i:]
        core = "".join(rng.choice(HYDROPHOBIC)
                       if rng.random() < 0.18 else c for c in core)
        rows += family(core, 1, f"pos{f:04d}")

        # negatives: same length distribution, and ~10% carry RGD by chance
        neg = "".join(rng.choice(AA) for _ in range(rng.randint(10, 17)))
        if rng.random() < 0.10:
            i = rng.randrange(1, len(neg) - 1)
            neg = neg[:i] + "RGD" + neg[i:]
        rows += family(neg, 0, f"neg{f:04d}")

    # assay noise: 7% of labels are simply wrong, as in any real screen
    for r in rows:
        if rng.random() < 0.07:
            r["label"] = 1 - r["label"]
    rng.shuffle(rows)
    for r in rows:
        # plate is 80% predictable from label -> a shortcut feature
        good = rng.random() < 0.8
        r["plate"] = ("P1" if r["label"] == 1 else "P2") if good else \
                     ("P2" if r["label"] == 1 else "P1")
    path = os.path.join(HERE, "peptide_binders.csv")
    with open(path, "w") as fh:
        fh.write("peptide_id,family,sequence,plate,label\n")
        for r in rows:
            fh.write(f"{r['peptide_id']},{r['family']},{r['sequence']},"
                     f"{r['plate']},{r['label']}\n")
    return rows


# ---------------------------------------------------------------- level 4 ---
def _helix_ca(n, origin, axis_offset, phase=0.0):
    """Ideal alpha-helix CA trace: r=2.3 A, 100 deg/residue, 1.5 A rise."""
    out = []
    for i in range(n):
        ang = math.radians(100.0 * i) + phase
        out.append((origin[0] + axis_offset[0] + 2.3 * math.cos(ang),
                    origin[1] + axis_offset[1] + 2.3 * math.sin(ang),
                    origin[2] + 1.5 * i))
    return out


def make_structure(rng):
    """A CA-only three-helix bundle (chain A) plus a peptide (chain B) docked
    against helix 1. Ideal geometry, not a real protein -- but the file format,
    the contact map and the interface logic are exactly the real thing."""
    chains = []
    helices = [
        _helix_ca(24, (0, 0, 0), (0.0, 0.0)),
        _helix_ca(24, (0, 0, 36), (10.0, 0.0)),   # antiparallel handled below
        _helix_ca(24, (0, 0, 0), (5.0, 9.0), phase=1.1),
    ]
    helices[1] = [(x, y, 36 - (z - 36)) for (x, y, z) in helices[1]]
    seqA, coordsA = [], []
    linker = "GSGSG"
    for hi, h in enumerate(helices):
        for (x, y, z) in h:
            seqA.append(rng.choice("AELKQRMILVF"))  # helix-favouring-ish
            coordsA.append((x, y, z))
        if hi < len(helices) - 1:
            last = coordsA[-1]
            for k, c in enumerate(linker, start=1):
                seqA.append(c)
                coordsA.append((last[0] + 1.2 * k, last[1] + 1.0 * k,
                                last[2] - 0.8 * k))
    # chain B: short peptide laid ~8 A off helix 1, roughly antiparallel
    seqB, coordsB = [], []
    for i in range(12):
        seqB.append("RGDWYKPTAEIL"[i])
        coordsB.append((-7.5 + 0.4 * math.sin(i), 1.2 * math.cos(i * 0.9),
                        4.0 + 1.5 * i))

    three = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
             "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
             "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
             "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR"}
    lines = ["REMARK  SYNTHETIC CA-ONLY MODEL GENERATED BY make_data.py",
             "REMARK  chain A = 3-helix bundle, chain B = docked peptide"]
    serial = 1
    for ch, seq, coords in (("A", seqA, coordsA), ("B", seqB, coordsB)):
        for i, (aa, (x, y, z)) in enumerate(zip(seq, coords), start=1):
            lines.append(
                f"ATOM  {serial:5d}  CA  {three[aa]} {ch}{i:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C")
            serial += 1
        lines.append(f"TER   {serial:5d}      {three[seq[-1]]} {ch}"
                     f"{len(seq):4d}")
        serial += 1
        chains.append((ch, "".join(seq)))
    lines.append("END")
    with open(os.path.join(HERE, "bundle_peptide.pdb"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return chains


def main():
    rng = random.Random(20260923)
    dna = make_dna(rng)
    prots = make_protein_pairs(rng)
    rows = make_binder_dataset(rng)
    chains = make_structure(rng)
    print(f"dna_sample.fasta      {len(dna)} sequences")
    print(f"proteins.fasta        {len(prots)} sequences")
    print(f"peptide_binders.csv   {len(rows)} rows "
          f"({sum(r['label'] for r in rows)} positives)")
    print(f"bundle_peptide.pdb    chains " +
          ", ".join(f"{c}={len(s)}aa" for c, s in chains))


if __name__ == "__main__":
    main()
