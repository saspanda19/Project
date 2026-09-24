# Level 4 — Structure: parsing, geometry, interfaces

**Time:** 2–3 days · **Prereqs:** Levels 1–3, numpy broadcasting, a little linear algebra · **New ideas:** fixed-column formats, contact maps, superposition, chirality

## Why this level exists

Three things bite everyone once:

1. **PDB is column-delimited, not whitespace-delimited.** `line.split()` works on tidy files and silently mangles real ones, where `-123.456-112.345` is two numbers with no space between them. One test in this level is exactly that line.
2. **A contact map without a sequence-separation filter is mostly diagonal.** Neighbouring residues are close because they're bonded; that tells you nothing about the fold.
3. **Kabsch superposition without the determinant sign fix happily superposes a protein onto its mirror image** and reports RMSD ≈ 0. Proteins are chiral. There's a test for that too.

## What you build

`parse_pdb`, `chain_sequence`, `coords`, `distance_matrix`, `contact_map`, `interface_residues`, `radius_of_gyration`, `kabsch_rmsd`.

```bash
python -m pytest level4_structure/test_level4.py -q
python level4_structure/solution.py
```

## About the data

`data/bundle_peptide.pdb` is a **synthetic CA-only model** — an ideal three-helix bundle (chain A) with a 12-residue peptide docked against helix 1 (chain B). Ideal geometry, made-up sequence. The file format, the distance math and the interface logic are the real thing; the molecule isn't.

Point it at a real structure as soon as your tests pass:

```bash
curl -O https://files.rcsb.org/download/1CRN.pdb     # crambin, 46 residues, all-atom
curl -O https://files.rcsb.org/download/1BRS.pdb     # barnase-barstar, a classic interface
python -c "import sys; sys.path.insert(0,'level4_structure'); import solution as s; \
           a=s.parse_pdb('1CRN.pdb'); print(len(a), s.chain_sequence(a,'A'))"
```

Real files will break assumptions your synthetic file let you keep: altloc indicators, insertion codes, multiple models, hydrogens, waters and ligands in `HETATM`, and missing residues that make `resseq` non-contiguous. **Each of those breaking is a lesson, not a bug in your code.** Fix them one at a time.

## Checkpoints

1. Your parser handles the run-together-coordinates line in the test.
2. Your contact map is symmetric and has a zeroed band of width `min_seq_sep` around the diagonal.
3. `cm[10, 14]` is True — in an α-helix, *i* and *i*+4 sit ~6.2 Å apart. If it isn't, your helix geometry or your cutoff is wrong.
4. `kabsch_rmsd(p, mirrored_p) > 0.5`. If it's ~0, you skipped the reflection fix.

## Questions to answer in your notes

- Why is the *long-range* contact count (|i−j| > 12) the number that distinguishes a bundle from one long helix? Look at both panels of `figures/level4_contacts.png`.
- A CA–CA cutoff of 8 Å is a convention. What does it approximate, and why do heavy-atom contact definitions (4.5–5 Å, any atom pair) give different interfaces?
- RMSD is dominated by the worst-fitting residues. When would you prefer TM-score or GDT-TS, and why do they behave differently for a two-domain protein whose hinge moved?

## Extensions

1. **Backbone dihedrals.** Compute φ/ψ from N, CA, C atoms of a real all-atom PDB and draw a Ramachandran plot. (This is why you need a real file — CA-only can't do it.)
2. **Secondary-structure assignment**, the poor man's DSSP: call a helix wherever CA(i)–CA(i+4) < 6.5 Å over a run of ≥4 residues; a sheet from long-range, near-parallel contact ladders. Compare with the PDB's own HELIX/SHEET records.
3. **Buried surface area.** Implement Shrake–Rupley (sphere sampling) and compute ΔSASA on binding for 1BRS. This is the standard way to size an interface.
4. **Interface hotspots.** Rank interface residues by contact count × a hydrophobicity weight, and compare your top 5 to the alanine-scanning hotspots reported for barnase–barstar.
5. **A real superposition task.** Take two structures of the same protein (e.g. an apo/holo pair), align their sequences with your **Level 2** code, superpose the aligned residues with Kabsch, and report per-residue deviation. That's three levels working together — and it's a real analysis, not an exercise.

## Graduating

Tests pass, your parser survives a real PDB entry you downloaded yourself, and you can say what your contact-map cutoff assumes.
