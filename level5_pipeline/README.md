# Level 5 — A binder-design campaign (mock models, real pipeline)

**Time:** 1–2 weeks · **Prereqs:** Levels 1–4 · **New ideas:** backend interfaces, funnels, filter ordering, diversity caps, provenance

## Why this level exists

In a real design campaign the models get the credit and the **pipeline** decides the outcome. What you filter on, in what order, at what thresholds, how you keep the shortlist diverse, and whether you can reproduce a run six months later — that's the engineering, and it's what a hiring manager or a PI actually asks you about.

So this level gives you a complete, working campaign with **mock** backends that run in under a second on a CPU. You build and debug the plumbing now; you swap in GPU models later without changing the pipeline.

> A mock model is not a small real model. It tells you your plumbing works. It tells you **nothing** about whether a design binds.

## The four interfaces

| stage | interface | real tool that fills it |
|---|---|---|
| backbone generation | `BackboneGenerator.generate(target, n, seed)` | RFdiffusion, RFdiffusion2, Chroma |
| sequence design | `SequenceDesigner.design(backbone, target, n, temperature, seed)` | ProteinMPNN, LigandMPNN |
| structure prediction | `StructurePredictor.predict(design, backbone, target)` | AlphaFold-Multimer, Boltz, Chai-1 |
| scoring | `Scorer.score(design, backbone, target)` | Rosetta InterfaceAnalyzer, FoldX, a learned scorer |

Swapping one is a constructor argument. `test_a_custom_backend_plugs_in_without_touching_the_pipeline` proves it with a poly-alanine designer; your RFdiffusion wrapper drops in the same way.

```bash
python level5_pipeline/pipeline.py
python level5_pipeline/pipeline.py --n-backbones 500 --seqs-per-backbone 4 --temperature 0.3
python -m pytest level5_pipeline/test_level5.py -q
```

Outputs land in `runs/demo/`: `all_designs.csv`, `shortlist.csv`, `shortlist.fasta`, and a `manifest.json` recording config, backends, funnel and versions. **A run without a manifest is an anecdote.**

## What to look at first

```
FUNNEL
  backbones generated                                     200
  backbones reaching >=2 hotspots                         186
  sequences designed                                      558
  structures predicted + scored                           558
  passed all filters                                      178
  shortlist (top 24, max 2/backbone)                       24

FIRST FAILING FILTER
  plddt          107      i_pae      157      contacts    48
  aggregation     58      net charge  10      iface hydro  0
```

Three habits are built into that report:

1. **Cheap filters run before expensive stages.** The hotspot pre-filter costs microseconds and removes candidates before folding. With real models, folding is minutes per design — filter ordering is most of your compute budget.
2. **Attribute each design's *first* failure.** It tells you which threshold is actually shaping your output. Here `i_pae` is doing the most work and `iface hydrophobicity` is doing none — that last one is a threshold you set and never checked, which is how dead filters survive in real code for years.
3. **Cap designs per backbone in the shortlist.** Rank on a single composite and your top 24 becomes six copies of one backbone and one failure mode. The cap is cheap insurance against a correlated shortlist.

## Checkpoints

1. Same seed, same shortlist — bit for bit. If not, you have an unseeded RNG somewhere; find it.
2. Impossible thresholds give an empty shortlist and a clean message, not a traceback. "Nothing survived" is a result.
3. The funnel never grows after the fan-out stage.
4. Higher `--temperature` gives more diverse sequences and worse scores. Plot that trade-off; it's the single most useful curve in a campaign.

## Questions to answer in your notes

- You have 24 hours of GPU time. Do you generate 10× more backbones, or design 10× more sequences per backbone, or fold everything with a better predictor? What would you measure to decide, and what do you actually know after a mock run? (Careful: the mock can't answer this. Say why.)
- Every threshold in `Filters` is a mock-scale number. Which ones are *structurally* right (the metric and direction) even though the value is fake?
- Your shortlist has 24 designs and you can order 12. What do you drop — the lowest composite, or the least diverse? Defend it.

## Real thresholds, for when you swap the models in

These are the community's rough starting points for a real campaign — treat them as a starting point to re-derive, not as truth:

- **pLDDT of the binder** ≥ 80 (AlphaFold-Multimer, binder chain only)
- **interface PAE** ≤ 10 Å, and worth trusting more than pLDDT for binding
- **pTM / ipTM** ≥ 0.5–0.6 for the complex
- **Rosetta ddG** ≤ −30 REU, and **SASA buried** ≥ 800–1000 Å²
- **self-consistency RMSD** ≤ 2 Å between the designed backbone and the refolded prediction — the single most predictive cheap filter in the literature
- **developability**: no free cysteines, no N-glycosylation sequons you didn't intend, no long hydrophobic or aggregation-prone stretch, pI away from your formulation pH

## Extensions, in the order I'd do them

1. **Add a self-consistency stage.** Refold each designed sequence, superpose against its designed backbone with your **Level 4** `kabsch_rmsd`, and filter on scRMSD. This is levels 4 and 5 clicking together, and it's the most valuable stage in the whole funnel.
2. **Make filtering data-driven.** Instead of hand-set thresholds, hold out designs with known labels and pick thresholds that maximise precision at a fixed shortlist size. Now your funnel has a **Level 3** evaluation behind it.
3. **Add caching.** Key each stage's output on a hash of its inputs and config, so a re-run costs nothing. With real models this is the difference between iterating and waiting.
4. **Add a partial-diffusion loop.** Take the top 10% of designs, re-diffuse around them, re-design, re-score. Track whether the composite improves across rounds, or whether you're just concentrating a shortlist.
5. **Swap in one real backend.** ProteinMPNN is the cheapest to start with (it runs on CPU for short binders). Keep every test green while you do it — that's the exercise.
6. **Parallelise the expensive stage** with a job array, and make the pipeline resumable from partial results. A campaign that can't resume will be re-run from scratch at 3am, and you'll be the one doing it.

## Graduating

You're done when you can hand someone `runs/<name>/manifest.json` and they can reproduce your shortlist exactly, and when you can explain which filter was doing the work.
