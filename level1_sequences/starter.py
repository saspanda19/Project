"""
LEVEL 1 - Sequence toolkit.  Fill in every TODO, then run:

    python -m pytest level1_sequences/test_level1.py -q

No Biopython. Writing these six functions by hand once is worth more than
importing them a hundred times, because every later level reuses the ideas
(indexing, frames, strandedness, off-by-one errors at sequence ends).
"""

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


def read_fasta(path):
    """Return [(name, description, sequence), ...]."""
    records = []
    name, description, sequence = None, "", []
    
    #Read the file line by line and strip each line
    
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
                
    # When a line starts with >, a new record begins. Save the record you were building, if there is one, then start a new one.
            
            if line.startswith(">"):
                if name is not None:
                    records.append((name, description, "".join(sequence)))
                parts = line[1:].split(None, 1)
                name = parts[0]
                description = parts[1] if len(parts) > 1 else ""
                sequence = []
            else:
                sequence.append(line.upper())

    # the last record never sees another '>', so it is banked here
    if name is not None:
        records.append((name, description, "".join(sequence)))

    return records
                
                
    """Rules that trip people up:
      * the header line is '>name rest of the description'
      * sequence lines wrap; join them
      * uppercase the sequence, strip whitespace
      * a record with no description gets ''
    """
    #raise NotImplementedError("TODO")


def gc_content(seq):
    """Fraction of G or C in seq, as a float in [0, 1].

    An empty sequence returns 0.0 (decide, document, and test edge cases --
    this is the habit that separates working code from demo code)."""
    if not seq:
        return 0.0
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq)
    #raise NotImplementedError("TODO")


def reverse_complement(seq):
    """Reverse complement of a DNA string. A<->T, C<->G, N->N."""
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(complement.get(base, "N") for base in seq.upper()[::-1])
    #raise NotImplementedError("TODO")


def translate(seq, stop_at_stop=False):
    """Translate seq in frame 0 into one-letter amino acids.

    '*' marks a stop codon. Trailing 1-2 leftover bases are ignored.
    If stop_at_stop, cut the protein at the first '*' (the '*' is dropped).
    """
    seq = seq.upper()
    aas = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        aa = CODON_TABLE.get(seq[i:i + 3], "X")
        if aa == "*" and stop_at_stop:
            break
        aas.append(aa)
    return "".join(aas)
    #raise NotImplementedError("TODO")


def find_orfs(seq, min_aa=30):
    """Find open reading frames on BOTH strands.

    Return a list of dicts, each:
        {"strand": "+" or "-", "frame": 0|1|2, "start": int, "end": int,
         "protein": str}

    Definition used here: ATG ... in-frame stop codon. 'start' and 'end' are
    0-based, half-open, ON THE STRAND SEARCHED (do not try to map minus-strand
    coordinates back to the plus strand yet -- that is an extension exercise).
    Keep only ORFs whose protein is at least min_aa long (excluding the stop).
    Sort by descending protein length.
    """
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
                            i = j
                            break
                        j += 3
                    else:
                        break
                i += 3
    orfs.sort(key=lambda o: -len(o["protein"]))
    return orfs
    #raise NotImplementedError("TODO")


def sliding_gc(seq, window=50, step=10):
    """Return [(midpoint_index, gc_fraction), ...] for windows along seq.

    The last partial window is dropped. If seq is shorter than window,
    return an empty list."""
    out = []
    for start in range(0, len(seq) - window + 1, step):
        out.append((start + window // 2, gc_content(seq[start:start + window])))
    return out
    #raise NotImplementedError("TODO")
