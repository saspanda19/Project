"""Tests for Level 3.  LADDER_MOD=solution to check the reference."""
import importlib.util
import os
import sys

import numpy as np

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

CSV = os.path.join(ROOT, "data", "peptide_binders.csv")


def test_load_dataset():
    df = m.load_dataset(CSV)
    assert list(df.columns) == ["peptide_id", "family", "sequence", "plate",
                                "label"]
    assert len(df) == 1200
    assert set(df["label"]) == {0, 1}
    assert df["family"].nunique() == 400


def test_aa_composition_is_a_fraction():
    v = m.aa_composition("AAAA")
    assert len(v) == 20
    assert abs(sum(v) - 1.0) < 1e-9
    assert abs(v[0] - 1.0) < 1e-9
    assert abs(sum(m.aa_composition("ACDEFGHIKLMNPQRSTVWY")) - 1.0) < 1e-9
    # composition must not encode length
    assert m.aa_composition("AC") == m.aa_composition("AACC")


def test_build_vocab_and_kmer_counts():
    seqs = ["RGDWY"] * 2 + ["AAAAAA"] * 6   # RGD appears 2x, AAA appears 24x
    vocab = m.build_vocab(seqs, k=3, min_count=5)
    assert "AAA" in vocab and "RGD" not in vocab, "min_count must be applied"
    assert vocab == sorted(vocab)
    counts = m.kmer_counts("AAAAA", k=3, vocab=vocab)
    assert len(counts) == len(vocab)
    assert counts[vocab.index("AAA")] == 3
    assert m.kmer_counts("WWWWW", k=3, vocab=vocab) == [0] * len(vocab)


def test_biophysical_features():
    f = m.biophysical_features("RGDWWW")
    assert f["length"] == 6
    assert f["has_RGD"] == 1
    assert abs(f["net_charge"] - 0.0) < 1e-9      # R +1, D -1
    assert abs(f["frac_aromatic"] - 0.5) < 1e-9   # WWW of 6
    assert m.biophysical_features("AAA")["has_RGD"] == 0
    assert m.biophysical_features("IIII")["mean_hydropathy"] > 0
    assert m.biophysical_features("RRRR")["mean_hydropathy"] < 0


def test_featurize_shapes():
    df = m.load_dataset(CSV).head(20)
    vocab = m.build_vocab(df["sequence"], k=3, min_count=1)
    X, names = m.featurize(df, vocab)
    assert X.shape == (20, 26 + len(vocab))
    assert len(names) == X.shape[1]
    Xp, namesp = m.featurize(df, vocab, use_plate=True)
    assert Xp.shape[1] == X.shape[1] + 1
    assert np.isfinite(X).all()


def test_grouped_split_has_no_shared_families():
    df = m.load_dataset(CSV)
    tr, te = m.grouped_split(df, test_size=0.25, seed=0)
    assert len(tr) + len(te) == len(df)
    assert not (set(tr["family"]) & set(te["family"])), \
        "this is the whole point of the level: no family in both halves"
    assert 0.15 < len(te) / len(df) < 0.35
    assert te["label"].nunique() == 2


def test_grouped_split_is_deterministic():
    df = m.load_dataset(CSV)
    a1, b1 = m.grouped_split(df, seed=1)
    a2, b2 = m.grouped_split(df, seed=1)
    assert list(b1["peptide_id"]) == list(b2["peptide_id"])
    _a3, b3 = m.grouped_split(df, seed=2)
    assert list(b1["peptide_id"]) != list(b3["peptide_id"])


def test_evaluate():
    y = [0, 0, 1, 1]
    perfect = m.evaluate(y, [0.1, 0.2, 0.8, 0.9])
    assert abs(perfect["auroc"] - 1.0) < 1e-9
    assert abs(perfect["f1"] - 1.0) < 1e-9
    assert perfect["n"] == 4 and perfect["n_pos"] == 2
    inverted = m.evaluate(y, [0.9, 0.8, 0.2, 0.1])
    assert abs(inverted["auroc"] - 0.0) < 1e-9
    chance = m.evaluate(y, [0.5, 0.5, 0.5, 0.5])
    assert abs(chance["auroc"] - 0.5) < 1e-9
    for key in ("auroc", "auprc", "accuracy", "precision", "recall", "f1"):
        assert key in perfect


def test_the_leak_is_real():
    """The headline result: a random split scores higher than a grouped one.

    If this fails, either your split or your featurization is wrong -- or you
    built the k-mer vocabulary on the full dataset instead of on train.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import train_test_split
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    df = m.load_dataset(CSV)

    def run(tr, te):
        vocab = m.build_vocab(tr["sequence"], k=3, min_count=5)
        Xtr, _ = m.featurize(tr, vocab)
        Xte, _ = m.featurize(te, vocab)
        clf = make_pipeline(StandardScaler(),
                            LogisticRegression(max_iter=2000, random_state=0))
        clf.fit(Xtr, tr["label"])
        return m.evaluate(te["label"], clf.predict_proba(Xte)[:, 1])["auroc"]

    rtr, rte = train_test_split(df, test_size=0.25, random_state=0,
                                stratify=df["label"])
    gtr, gte = m.grouped_split(df, test_size=0.25, seed=0)
    leaky, honest = run(rtr, rte), run(gtr, gte)
    assert honest > 0.6, "there is real signal; you should beat chance"
    assert leaky > honest + 0.05, \
        f"expected leakage to inflate the score (leaky={leaky:.3f}, honest={honest:.3f})"
