# Level 3 — ML on sequences, done honestly

**Time:** 2–4 days · **Prereqs:** Levels 1–2, basic numpy/pandas · **New ideas:** featurization, leakage, grouped splits, confounds, proxy features

## Why this level exists

Getting a model to fit is the easy half. The hard half — the half that decides whether anything you build is worth acting on — is knowing whether the number you just printed is real.

`data/peptide_binders.csv` is booby-trapped on purpose:

| trap | what it is | what it does to you |
|---|---|---|
| near-duplicates | binders come in families of 3 near-identical siblings | a random split puts siblings in train *and* test |
| batch confound | `plate` predicts the label 80% of the time | a model handed `plate` learns the plate |
| proxy feature | the motif insertion makes positives 3 residues longer | the model learns `length` and looks fine |
| label noise | 7% of labels are flipped | your ceiling is not 1.0, and chasing it is overfitting |

All four exist in real screening data. The last two are the ones people miss.

## What you build

Featurization (`aa_composition`, `kmer_counts`, `biophysical_features`), a **grouped** train/test split, and an `evaluate` that reports AUPRC next to AUROC.

```bash
python -m pytest level3_ml/test_level3.py -q
python level3_ml/solution.py
```

## The table you should reproduce

```
setting                        AUROC   AUPRC      F1  recall
random split (LEAKY)           0.841   0.845   0.794   0.770
grouped split (honest)         0.708   0.717   0.667   0.697
grouped + plate feature        0.848   0.857   0.800   0.803
grouped, 26 features only      0.781   0.812   0.708   0.750
grouped, 26 feat, no RGD       0.778   0.813   0.710   0.750
```

Four things worth sitting with:

1. **0.841 → 0.708.** The same model, same data, honest split. The 0.841 is the number that gets into a slide deck.
2. **The plate feature restores 0.848.** A confound can look exactly like success. Nothing in the metrics tells you it's a plate.
3. **26 hand-built features beat 3000+ k-mers** (0.781 vs 0.708). With 900 training rows, a huge sparse feature space mostly fits noise — look at the printed coefficients of both models and compare.
4. **Erasing the RGD motif costs almost nothing** (0.781 → 0.778). The model wasn't using the motif; it was using `length`, which is a *proxy* for the motif. This is the failure mode that survives a clean test set and dies in the lab.

## Checkpoints

1. Your grouped split shares no family between halves — the test asserts it.
2. You build the k-mer vocabulary on the **training set only**. Building it on all of `df` is leakage that costs about 0.02 AUROC here and much more in real data.
3. `aa_composition` returns fractions, so `"AC"` and `"AACC"` give the same vector.
4. You report AUPRC as well as AUROC and can say when they disagree.

## Questions to answer in your notes

- Which of the four traps would you have caught if nobody had told you it was there? Be honest; that's the point.
- The honest AUROC is ~0.71 and 7% of labels are wrong. What *is* the ceiling here, and how would you estimate it?
- You're handed a model with test AUROC 0.93 on a screening dataset. List the five questions you ask before believing it. (Split strategy, duplicate/homology handling, feature provenance, class balance, and what the baseline is.)

## Extensions

1. **Cross-validate properly.** Replace the single grouped split with `GroupKFold` and report mean ± std. One split's number has a standard error you're currently ignoring.
2. **Beat the linear model.** Random forest and gradient boosting, same split. Does the gap between leaky and honest widen? (It usually does — higher-capacity models eat leakage faster.)
3. **Permutation importance** on the honest split, instead of reading raw coefficients — coefficients of correlated, standardized features lie.
4. **Kill the length proxy.** Match the length distributions between classes (subsample, or add `length` as a covariate and regress it out) and re-measure. How much signal is left when the model can't cheat?
5. **Homology-aware splitting.** Cluster sequences by your Level 2 alignment identity at 70%, then split on clusters instead of on the `family` column you were handed. This is what real benchmarks (and reviewers) expect.

## Graduating

You're done when you can explain, without notes, why the first row of that table is the most dangerous number in the repo.
