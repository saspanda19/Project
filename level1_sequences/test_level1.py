"""Tests for Level 1.

By default they run against starter.py (your work). To check the reference:

    LADDER_MOD=solution python -m pytest level1_sequences/test_level1.py -q
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)


def _load(name):
    """Load this level's starter/solution BY PATH, so levels that share module
    names don't collide when pytest collects the whole repo at once."""
    spec = importlib.util.spec_from_file_location(
        f"{os.path.basename(HERE)}_{name}", os.path.join(HERE, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


m = _load(os.environ.get("LADDER_MOD", "starter"))

FASTA = os.path.join(ROOT, "data", "dna_sample.fasta")


def test_read_fasta_shape():
    recs = m.read_fasta(FASTA)
    assert len(recs) == 6
    names = [r[0] for r in recs]
    assert names == ["gene_A", "gene_B", "gene_C", "gene_D", "gene_E", "gene_F"]
    for name, desc, seq in recs:
        assert desc.startswith("synthetic|")
        assert set(seq) <= set("ACGT"), "sequence lines must be joined+uppercased"
        assert "\n" not in seq and " " not in seq


def test_gc_content():
    assert m.gc_content("GGCC") == 1.0
    assert m.gc_content("ATAT") == 0.0
    assert abs(m.gc_content("ACGT") - 0.5) < 1e-12
    assert m.gc_content("") == 0.0
    assert abs(m.gc_content("acgt") - 0.5) < 1e-12, "handle lowercase"


def test_gc_content_matches_generator_targets():
    # each record's description carries the GC the generator aimed for
    for _n, desc, seq in m.read_fasta(FASTA):
        target = float(desc.split("gc_target=")[1].split("|")[0])
        assert abs(m.gc_content(seq) - target) < 0.12


def test_reverse_complement():
    assert m.reverse_complement("ATGC") == "GCAT"
    assert m.reverse_complement("AAAA") == "TTTT"
    assert m.reverse_complement("N") == "N"
    s = "ATGCCGTTAGC"
    assert m.reverse_complement(m.reverse_complement(s)) == s


def test_translate():
    assert m.translate("ATGGCCTAA") == "MA*"
    assert m.translate("ATGGCCTAA", stop_at_stop=True) == "MA"
    assert m.translate("ATGGCCTAAGG") == "MA*", "trailing 2 bases are ignored"
    assert m.translate("") == ""


def test_find_orfs_recovers_the_planted_gene():
    for _n, desc, seq in m.read_fasta(FASTA):
        expected_aa = int(desc.split("orf_aa=")[1])
        orfs = m.find_orfs(seq, min_aa=30)
        assert orfs, "every synthetic gene contains one long ORF"
        top = orfs[0]
        assert top["protein"].startswith("M")
        assert "*" not in top["protein"], "the stop codon is not part of the protein"
        # >= (not ==): an upstream in-frame ATG legitimately extends the ORF,
        # which is exactly why "first ATG" and "longest ORF" are different
        # gene-calling rules. Pick one, document it, stay consistent.
        assert expected_aa <= len(top["protein"]) <= expected_aa + 60
        assert top["strand"] == "+"
        assert seq[top["start"]:top["start"] + 3] == "ATG"
        assert (top["end"] - top["start"]) % 3 == 0
        assert len(top["protein"]) == (top["end"] - top["start"]) // 3 - 1


def test_find_orfs_min_aa_filter():
    _n, _d, seq = m.read_fasta(FASTA)[0]
    assert len(m.find_orfs(seq, min_aa=1)) >= len(m.find_orfs(seq, min_aa=30))
    assert all(len(o["protein"]) >= 40 for o in m.find_orfs(seq, min_aa=40))


def test_find_orfs_searches_both_strands():
    # a minus-strand-only ORF: build it reverse-complemented into the input
    protein_dna = "ATG" + "GCT" * 35 + "TAA"
    seq = "TTTTT" + m.reverse_complement(protein_dna) + "TTTTT"
    orfs = m.find_orfs(seq, min_aa=30)
    assert any(o["strand"] == "-" and len(o["protein"]) == 36 for o in orfs)


def test_sliding_gc():
    assert m.sliding_gc("ACGT", window=50) == []
    pts = m.sliding_gc("GC" * 50, window=20, step=10)
    assert pts and all(abs(v - 1.0) < 1e-12 for _i, v in pts)
    assert [i for i, _v in pts] == sorted(i for i, _v in pts)
    assert all(0 <= i < 100 for i, _v in pts)
