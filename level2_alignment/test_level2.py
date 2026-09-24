"""Tests for Level 2.  LADDER_MOD=solution to check the reference."""
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
from blosum62 import BLOSUM62  # noqa: E402


def test_score_pair():
    assert m.score_pair("W", "W") == 11
    assert m.score_pair("L", "I") == 2
    assert m.score_pair("W", "D") == -4
    assert m.score_pair("A", "?") == 0


def test_score_alignment_referee():
    assert m.score_alignment("ACDEF", "ACDEF") == 4 + 9 + 6 + 5 + 6
    # one gap of length 3: open 11 + 3 * extend 1
    assert m.score_alignment("AAAAA", "AA---") == 4 + 4 - (11 + 3)
    # two separate gaps of length 1 cost more than one gap of length 2
    one = m.score_alignment("AAAA", "A--A")
    two = m.score_alignment("AAAA", "A-A-")
    assert one > two, "affine gaps must reward keeping a gap contiguous"


def test_score_alignment_rejects_gap_gap():
    try:
        m.score_alignment("A-", "A-")
    except ValueError:
        return
    raise AssertionError("gap aligned to gap should raise ValueError")


def test_nw_identical_sequences():
    s = "ACDEFGHIKLMNPQRSTVWY"
    score, a, b = m.needleman_wunsch(s, s)
    assert a == s and b == s
    assert score == sum(BLOSUM62[(c, c)] for c in s)


def test_nw_score_matches_referee():
    pairs = [("ACDEF", "ACDEF"), ("ACDEF", "AEF"), ("WWWW", "W"),
             ("MKTAYIAKQRQISFVK", "MKTAYIAKQRQISFVKSHFSRQ"),
             ("PLLKKMTSRTA", "PKKMTSRTA"), ("A", "Y")]
    for x, y in pairs:
        score, ax, ay = m.needleman_wunsch(x, y)
        assert len(ax) == len(ay)
        assert ax.replace("-", "") == x and ay.replace("-", "") == y, \
            "the traceback must reproduce the input sequences"
        assert score == m.score_alignment(ax, ay), \
            f"DP said {score} but the alignment scores {m.score_alignment(ax, ay)}"


def test_nw_is_global():
    x, y = "AAAACDEFAAAA", "CDEF"
    score, ax, ay = m.needleman_wunsch(x, y)
    assert len(ax) == len(x), "a global alignment spans the whole longer seq"
    assert ay.count("-") == 8


def test_affine_prefers_one_long_gap():
    x = "ACDEFGHIKLMNPQRST"
    y = "ACDEFRST"           # a single contiguous deletion
    score, ax, ay = m.needleman_wunsch(x, y)
    runs = sum(1 for k in range(len(ay))
               if ay[k] == "-" and (k == 0 or ay[k - 1] != "-"))
    assert runs == 1, f"expected one gap run, got {runs}: {ay}"


def test_sw_finds_the_embedded_domain():
    domain = "MKTAYIAKQRQISFVKSHFSRQ"
    x = "PPPPPPPP" + domain + "PPPP"
    y = "WWWWWWWWWWWW" + domain + "WW"
    score, ax, ay, si, sj = m.smith_waterman(x, y)
    assert ax == ay == domain
    assert si == 8 and sj == 12
    assert score == sum(BLOSUM62[(c, c)] for c in domain)


def test_sw_never_scores_below_zero():
    score, ax, ay, si, sj = m.smith_waterman("WWWW", "PPPP")
    assert score >= 0


def test_sw_beats_nw_on_local_similarity():
    domain = "MKTAYIAKQRQISFVKSHFSRQ"
    x = "PPPPPPPP" + domain + "PPPP"
    y = "WWWWWWWWWWWW" + domain + "WW"
    local = m.smith_waterman(x, y)[0]
    glob = m.needleman_wunsch(x, y)[0]
    assert local > glob


def test_percent_identity():
    assert m.percent_identity("ACDE", "ACDE") == 1.0
    assert m.percent_identity("ACDE", "ACDF") == 0.75
    assert m.percent_identity("AC-E", "ACDE") == 0.75, "gap columns count in the denominator"
    assert m.percent_identity("", "") == 0.0


def test_real_data_ordering():
    """close > mid > far identity, on the synthetic family in data/."""
    def read(path):
        recs, name, chunks = {}, None, []
        for line in open(path):
            line = line.strip()
            if line.startswith(">"):
                if name:
                    recs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
        recs[name] = "".join(chunks)
        return recs

    seqs = read(os.path.join(ROOT, "data", "proteins.fasta"))
    ids = []
    for n in ("prot_close", "prot_mid", "prot_far"):
        _s, a, b = m.needleman_wunsch(seqs["prot_parent"], seqs[n])
        ids.append(m.percent_identity(a, b))
    assert ids[0] > ids[1] > ids[2]
    assert ids[0] > 0.8 and ids[2] < 0.7
