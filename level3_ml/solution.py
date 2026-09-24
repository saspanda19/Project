"""LEVEL 3 reference solution.

    python level3_ml/solution.py

Runs the same model four ways and prints a table that is the whole point of
this level: random split, grouped split, grouped split + plate feature,
grouped split with the motif removed.
"""
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
def _sibling(name):
    """Import this level's own module by path (levels share module names)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        f"{os.path.basename(HERE)}_{name}", os.path.join(HERE, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

_st = _sibling("starter")
AA, KD_HYDROPATHY, CHARGE = _st.AA, _st.KD_HYDROPATHY, _st.CHARGE

AROMATIC = set("FWY")
HYDROPHOBIC = set("AVILMFWY")


def load_dataset(path):
    return pd.read_csv(path)


def aa_composition(seq):
    n = len(seq)
    if n == 0:
        return [0.0] * 20
    c = Counter(seq)
    return [c.get(a, 0) / n for a in AA]


def build_vocab(sequences, k=3, min_count=5):
    c = Counter()
    for s in sequences:
        for i in range(len(s) - k + 1):
            c[s[i:i + k]] += 1
    return sorted(km for km, n in c.items() if n >= min_count)


def kmer_counts(seq, k=3, vocab=None):
    if vocab is None:
        raise ValueError("pass an explicit vocab built from the training set")
    idx = {km: i for i, km in enumerate(vocab)}
    out = [0] * len(vocab)
    for i in range(len(seq) - k + 1):
        j = idx.get(seq[i:i + k])
        if j is not None:
            out[j] += 1
    return out


def biophysical_features(seq):
    n = max(len(seq), 1)
    return {
        "length": float(len(seq)),
        "mean_hydropathy": sum(KD_HYDROPATHY.get(c, 0) for c in seq) / n,
        "net_charge": float(sum(CHARGE.get(c, 0) for c in seq)),
        "frac_aromatic": sum(c in AROMATIC for c in seq) / n,
        "frac_hydrophobic": sum(c in HYDROPHOBIC for c in seq) / n,
        "has_RGD": float("RGD" in seq),
    }


BIO_KEYS = ["length", "mean_hydropathy", "net_charge", "frac_aromatic",
            "frac_hydrophobic", "has_RGD"]


def featurize(df, vocab, k=3, use_plate=False):
    rows, names = [], [f"comp_{a}" for a in AA] + BIO_KEYS + \
        [f"kmer_{v}" for v in vocab]
    if use_plate:
        names.append("plate_is_P1")
    for _, r in df.iterrows():
        s = r["sequence"]
        bio = biophysical_features(s)
        vec = aa_composition(s) + [bio[k_] for k_ in BIO_KEYS] + \
            kmer_counts(s, k=k, vocab=vocab)
        if use_plate:
            vec = vec + [1.0 if r["plate"] == "P1" else 0.0]
        rows.append(vec)
    return np.asarray(rows, dtype=float), names


def grouped_split(df, group_col="family", test_size=0.25, seed=0):
    rng = np.random.default_rng(seed)
    groups = list(df[group_col].unique())
    rng.shuffle(groups)
    n_test = max(1, int(round(len(groups) * test_size)))
    test_groups = set(groups[:n_test])
    mask = df[group_col].isin(test_groups)
    return df[~mask].reset_index(drop=True), df[mask].reset_index(drop=True)


def random_split(df, test_size=0.25, seed=0):
    from sklearn.model_selection import train_test_split
    return train_test_split(df, test_size=test_size, random_state=seed,
                            stratify=df["label"])


def evaluate(y_true, y_score, threshold=0.5):
    from sklearn.metrics import (accuracy_score, average_precision_score,
                                 f1_score, precision_score, recall_score,
                                 roc_auc_score)
    y_true = np.asarray(y_true)
    y_pred = (np.asarray(y_score) >= threshold).astype(int)
    return {
        "auroc": roc_auc_score(y_true, y_score),
        "auprc": average_precision_score(y_true, y_score),
        "accuracy": accuracy_score(y_true, y_pred),
        "precision": precision_score(y_true, y_pred, zero_division=0),
        "recall": recall_score(y_true, y_pred, zero_division=0),
        "f1": f1_score(y_true, y_pred, zero_division=0),
        "n_pos": int(y_true.sum()), "n": int(len(y_true)),
    }


# ------------------------------------------------------------------ demo ---
def _fit_eval(train, test, use_plate=False, strip_motif=False, no_kmers=False,
              seed=0):
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    if strip_motif:
        train = train.copy()
        test = test.copy()
        train["sequence"] = train["sequence"].str.replace("RGD", "GGG")
        test["sequence"] = test["sequence"].str.replace("RGD", "GGG")

    # vocabulary from the TRAINING SET ONLY - building it on all of df is a
    # quiet, very common form of leakage
    vocab = [] if no_kmers else build_vocab(train["sequence"], k=3, min_count=5)
    Xtr, names = featurize(train, vocab, use_plate=use_plate)
    Xte, _ = featurize(test, vocab, use_plate=use_plate)
    clf = make_pipeline(StandardScaler(),
                        LogisticRegression(max_iter=2000, C=1.0,
                                           random_state=seed))
    clf.fit(Xtr, train["label"])
    score = clf.predict_proba(Xte)[:, 1]
    return evaluate(test["label"], score), clf, names, vocab


def main():
    df = load_dataset(os.path.join(ROOT, "data", "peptide_binders.csv"))
    print(f"{len(df)} peptides, {df['label'].sum()} binders, "
          f"{df['family'].nunique()} families\n")

    runs = []
    tr, te = random_split(df)
    runs.append(("random split (LEAKY)", _fit_eval(tr, te)[0]))

    tr, te = grouped_split(df)
    res, clf, names, vocab = _fit_eval(tr, te)
    runs.append(("grouped split (honest)", res))
    runs.append(("grouped + plate feature", _fit_eval(tr, te, use_plate=True)[0]))
    res26, clf26, names26, _v26 = _fit_eval(tr, te, no_kmers=True)
    runs.append(("grouped, 26 features only", res26))
    runs.append(("grouped, 26 feat, no RGD",
                 _fit_eval(tr, te, no_kmers=True, strip_motif=True)[0]))

    hdr = f"{'setting':<28}{'AUROC':>8}{'AUPRC':>8}{'F1':>8}{'recall':>8}"
    print(hdr)
    print("-" * len(hdr))
    for name, r in runs:
        print(f"{name:<28}{r['auroc']:>8.3f}{r['auprc']:>8.3f}"
              f"{r['f1']:>8.3f}{r['recall']:>8.3f}")

    print("\nRead that table, not the first row. The leaky split looks better")
    print("than the honest one because test peptides have siblings in train.")
    print("The plate row shows a confound that would survive peer review.")
    print("Dropping 3000+ junk k-mers barely costs anything - and the last")
    print("row, with the RGD motif erased, is what the model actually knew.\n")

    for title, model, nms in (
            ("3000+ k-mer model - the coefficients are noise", clf, names),
            ("26-feature model - now look at what it leans on", clf26, names26)):
        coefs = model[-1].coef_[0]
        order = np.argsort(np.abs(coefs))[::-1][:8]
        print(f"\n{title}:")
        for i in order:
            print(f"  {nms[i]:<18}{coefs[i]:+.3f}")
    print("\nIf 'length' is at the top, ask why. The generator inserts RGD")
    print("into positives, which makes them 3 residues longer - so length is a")
    print("PROXY for the motif. That is how a model passes your test set and")
    print("fails in the lab. Fixing it is exercise 4 in the README.")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from sklearn.metrics import roc_curve
    except Exception:
        return
    figdir = os.path.join(ROOT, "figures")
    os.makedirs(figdir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(5.2, 5))
    tr_r, te_r = random_split(df)
    for label, (a, b), kw in (("random split (leaky)", (tr_r, te_r), {}),
                              ("grouped split (honest)", grouped_split(df), {})):
        r, c, _n, _v = _fit_eval(a, b)
        X = featurize(b, build_vocab(a["sequence"]))[0]
        s = c.predict_proba(X)[:, 1]
        fpr, tpr, _ = roc_curve(b["label"], s)
        ax.plot(fpr, tpr, label=f"{label} (AUC {r['auroc']:.3f})", lw=1.8)
    ax.plot([0, 1], [0, 1], "k--", lw=0.8)
    ax.set_xlabel("false positive rate")
    ax.set_ylabel("true positive rate")
    ax.set_title("Level 3 - the same model, two splits")
    ax.legend(fontsize=8, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "level3_roc.png"), dpi=150)
    print(f"\nwrote {os.path.join(figdir, 'level3_roc.png')}")


if __name__ == "__main__":
    main()
