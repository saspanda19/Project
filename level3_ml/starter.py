"""
LEVEL 3 - Machine learning on peptide sequences, done honestly.

    python -m pytest level3_ml/test_level3.py -q

The dataset (data/peptide_binders.csv) is booby-trapped on purpose:

  1. REAL SIGNAL     an 'RGD' motif plus hydrophobic enrichment -> binder
  2. NEAR-DUPLICATES binders come in families of 3 mutated siblings
  3. BATCH CONFOUND  the 'plate' column predicts the label 80% of the time

A random split leaks family members across train/test, and a model handed
'plate' learns the plate instead of the biology. Both produce a beautiful
number. Your job is to produce an honest one.
"""

AA = "ACDEFGHIKLMNPQRSTVWY"
KD_HYDROPATHY = {  # Kyte-Doolittle
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1, "H": 0.1}


def load_dataset(path):
    """Read the CSV into a pandas DataFrame with columns
    peptide_id, family, sequence, plate, label."""
    raise NotImplementedError("TODO")


def aa_composition(seq):
    """Return a length-20 list: the FRACTION of each amino acid in AA order.

    Fraction, not count - otherwise your model learns peptide length, which in
    this dataset differs slightly between classes. That is a third confound,
    and it is the one people miss."""
    raise NotImplementedError("TODO")


def kmer_counts(seq, k=3, vocab=None):
    """Counts of each k-mer in `vocab` (a list of k-mer strings).

    With k=3 the full vocabulary is 8000 columns for a 1200-row dataset, so
    build `vocab` from the TRAINING SET ONLY and pass it in. Building it from
    the whole dataset is leakage - subtle, common, and fatal."""
    raise NotImplementedError("TODO")


def build_vocab(sequences, k=3, min_count=5):
    """k-mers occurring at least min_count times across `sequences`, sorted."""
    raise NotImplementedError("TODO")


def biophysical_features(seq):
    """Return a dict of interpretable features:
        length, mean_hydropathy, net_charge, frac_aromatic (FWY),
        frac_hydrophobic (AVILMFWY), has_RGD (0/1)
    """
    raise NotImplementedError("TODO")


def featurize(df, vocab, k=3, use_plate=False):
    """Return (X, feature_names) as a numpy array.

    Concatenate: aa_composition (20) + biophysical (6) + kmer_counts(len(vocab)).
    If use_plate, append a one-hot plate column - you will use this to MEASURE
    the confound, not to win.
    """
    raise NotImplementedError("TODO")


def grouped_split(df, group_col="family", test_size=0.25, seed=0):
    """Split so that no group appears in both halves. Return (train_df, test_df).

    Compare against a naive random split -- that comparison IS this level.
    """
    raise NotImplementedError("TODO")


def evaluate(y_true, y_score, threshold=0.5):
    """Return a dict: auroc, auprc, accuracy, precision, recall, f1, n_pos, n.

    Use sklearn.metrics. Report AUPRC alongside AUROC: on imbalanced data
    AUROC flatters a model that AUPRC exposes."""
    raise NotImplementedError("TODO")
