# Level 1 — Sequence toolkit

**Time:** 2–4 hours · **Prereqs:** basic Python (loops, dicts, slicing) · **New ideas:** frames, strandedness, half-open intervals, edge cases

## Why this level exists

Every computational biology stack sits on top of string handling that looks trivial and isn't. Frames, strands, and off-by-one errors at sequence ends are where beginners silently lose days. You write it once by hand; after that you may use Biopython with a clear conscience.

## What you build

Six functions in `starter.py`: `read_fasta`, `gc_content`, `reverse_complement`, `translate`, `find_orfs`, `sliding_gc`.

## Workflow

```bash
python data/make_data.py                              # once, generates the data
python -m pytest level1_sequences/test_level1.py -q   # red -> green as you go
python level1_sequences/solution.py                   # the demo, once you're done
```

Read `solution.py` **after** your tests pass, not before. Comparing two working implementations teaches; copying one doesn't.

## Checkpoints

1. `read_fasta` handles wrapped sequence lines and a missing description.
2. `gc_content("")` returns something you chose deliberately, not a `ZeroDivisionError`.
3. `find_orfs` finds ORFs on the minus strand. If you only get `+` hits, you forgot to search the reverse complement.
4. Your ORF protein length equals `(end - start) // 3 - 1`. If it's off by one, you're including the stop codon.

## Questions to answer in your notes

- **"First ATG" vs "longest ORF."** The test allows your top ORF to be *longer* than the planted gene, because an upstream in-frame ATG extends it. Which rule does a real gene caller use, and why is the answer "it depends on the organism"?
- Your ORF coordinates are relative to the strand you searched. Write down the formula that maps a minus-strand hit back to plus-strand coordinates. (Then implement it — extension 1.)
- GC content is a single number over a whole sequence. What does the sliding-window plot in `figures/level1_gc.png` show that the single number hides?

## Extensions (do at least one)

1. Map minus-strand ORF coordinates back onto the plus strand and verify by re-extracting the sequence.
2. Add a codon-usage table per sequence and compute a codon adaptation index against `gene_B` as the reference set.
3. Add `--min-aa` and `--fasta` command-line flags with `argparse`, so the module is a real tool, not a script.
4. Swap the synthetic data for a real genome: download an *E. coli* K-12 FASTA from NCBI and see which of your assumptions break at 4.6 Mb (hint: memory, and `N` characters).

## Graduating

You're done when the tests pass, you've written the coordinate-mapping formula down, and the sliding-GC figure makes sense to you.
