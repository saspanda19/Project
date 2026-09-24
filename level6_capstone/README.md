# Level 6 — Capstone

**Time:** 3–6 weeks · **Prereqs:** Levels 1–5 · **New ideas:** none. That's the point — this level is about doing real work with real data and defending it.

There is no starter code here on purpose. Levels 1–5 gave you the components; a capstone is where you choose a question, and the choosing is the skill.

## Pick one track

### Track A — A real binder-design campaign

Replace the Level 5 mocks with real models against a real target (PD-L1, IL-7Rα and TrkA all have published designed binders you can benchmark against).

- **Backbones:** RFdiffusion, target hotspots taken from the structure of a known complex.
- **Sequences:** ProteinMPNN (CPU-viable for short binders; GPU is faster).
- **Refold:** AlphaFold-Multimer, Boltz-2 or Chai-1 (GPU; Colab is enough for a few hundred designs).
- **Filter:** self-consistency RMSD via your Level 4 Kabsch, interface PAE, ipTM, plus developability.
- **Deliverable:** 20–50 ranked designs, a manifest, a funnel figure, and a page on what you'd order and why.

**The honest framing:** you cannot validate these computationally. Say so, then show that your pipeline recovers a *published* binder's metrics when you run it through — that's the closest thing to a positive control you have.

### Track B — A benchmark nobody can cheat

Take a public dataset (SKEMPI ΔΔG on mutations, a peptide-MHC binding set, a solubility set) and build the benchmark you wished existed: homology-aware splits, a stated baseline, a leakage audit, and a leaderboard script.

- Cluster sequences at 30%/50%/70% identity (your Level 2 alignment, or MMseqs2 for scale) and split on clusters.
- Report a **trivial baseline** (composition-only logistic regression) beside every model. Most published gains are smaller than this baseline.
- Audit the dataset for near-duplicates across the split and report how many you found.
- **Deliverable:** a repo anyone can run, a results table with error bars, and a short write-up of what the leakage audit found.

This is the least glamorous track and the one that most often produces something other people actually use.

### Track C — Bridging MD and design

If you already work in molecular dynamics, this is the track where your existing skills compound instead of resetting.

Take a peptide–membrane or peptide–protein system, run short CGMD or all-atom simulations of a handful of designed sequences, and ask whether any **cheap** descriptor predicts a simulation-derived observable (insertion depth, tilt-angle distribution, contact lifetime, PMF at a single point).

- Design 10–20 variants with the Level 5 pipeline; simulate; measure one observable properly, with error bars from independent replicas.
- Build the Level 3 model that predicts the observable from sequence. It will be data-poor — say what *n* you'd need.
- **Deliverable:** a figure of predicted vs simulated, a frank statement of statistical power, and a note on what a sequence-to-simulation surrogate would need to be useful.

The scientific value here is the negative result done well: with *n* = 15 you will probably not have a model, and being able to say precisely why is more impressive than a fitted line through nine points.

## Rubric — what makes any of these good

| | weak | strong |
|---|---|---|
| **question** | "apply X to Y" | a question with a wrong answer you'd recognise |
| **baseline** | none | a trivial baseline reported beside every result |
| **splits** | random | homology-aware, with the audit shown |
| **uncertainty** | a single number | error bars, and a stated source for them |
| **negative results** | omitted | reported, with what you learned |
| **reproducibility** | "run the notebook" | pinned deps, a seed, a manifest, one command |
| **limits** | "future work" | a specific statement of what would falsify your claim |

## Compute

Everything through Level 5 runs on a laptop CPU. For Track A you need a GPU for the folding step: Colab's free tier handles a few hundred designs, and university cluster time is the natural next step. For Track C the simulation dominates — plan replicas before you plan analysis.

## Write-up

Five pages, in this order: question, data, method, what you found, what you don't know. Put the funnel figure and the results table on page one. Whoever reads it — a PI, a hiring manager, a reviewer — is checking whether you know the difference between a result and a plot.

## Graduating

Someone else clones your repo, runs one command, and gets your numbers. Then they ask the hardest question they can think of and you have an answer, including "we can't tell from this."
