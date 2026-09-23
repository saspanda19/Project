"""
LEVEL 4 - Protein structure: parsing, geometry, interfaces.

    python -m pytest level4_structure/test_level4.py -q

data/bundle_peptide.pdb is a CA-only model: chain A is a three-helix bundle,
chain B a 12-residue peptide docked against helix 1. The coordinates are
idealized, but the file format, the distance math and the interface logic are
exactly what you do on a real structure -- and the README shows you how to
point this code at a real PDB entry in one line.

Coordinates are in angstroms. Everything below is numpy-friendly; you may use
numpy freely.
"""

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def parse_pdb(path):
    """Parse ATOM records into a list of dicts, in file order:

        {"serial": int, "atom": str, "resname": str, "chain": str,
         "resseq": int, "x": float, "y": float, "z": float}

    PDB is a COLUMN-DELIMITED format, not whitespace-delimited. Splitting on
    whitespace works until it doesn't -- coordinates like '-123.456-12.345'
    run together in real files. Slice by column:

        serial  7-11   atom 13-16   resname 18-20   chain 22
        resseq  23-26  x 31-38      y 39-46        z 47-54
    (those are 1-based inclusive ranges from the PDB spec; Python slices are
    0-based half-open, which is the off-by-one you are about to hit.)

    Skip every line that does not start with 'ATOM'.
    """
    raise NotImplementedError("TODO")


def chain_sequence(atoms, chain):
    """One-letter sequence of `chain`, in residue-number order, from CA atoms.
    Unknown residue names become 'X'."""
    raise NotImplementedError("TODO")


def coords(atoms, chain=None, atom_name="CA"):
    """Return an (N, 3) numpy array of coordinates, filtered by chain and
    atom name. chain=None means all chains."""
    raise NotImplementedError("TODO")


def distance_matrix(a, b=None):
    """Pairwise Euclidean distances. If b is None, use a against itself.
    Return shape (len(a), len(b)). Vectorize it -- no Python double loop."""
    raise NotImplementedError("TODO")


def contact_map(atoms, chain, cutoff=8.0, min_seq_sep=3):
    """Boolean (N, N) CA-CA contact map for one chain.

    Residues closer than min_seq_sep in sequence are NOT contacts (they are
    trivially close because they are covalently near each other), so zero out
    the band around the diagonal."""
    raise NotImplementedError("TODO")


def interface_residues(atoms, chain_a="A", chain_b="B", cutoff=8.0):
    """Return (residues_in_a, residues_in_b): sorted lists of residue numbers
    with at least one CA-CA contact across the interface."""
    raise NotImplementedError("TODO")


def radius_of_gyration(xyz):
    """Rg = sqrt(mean(|r_i - r_center|^2)). One number, in angstroms."""
    raise NotImplementedError("TODO")


def kabsch_rmsd(p, q):
    """Optimal-superposition RMSD between two (N, 3) point sets.

    1. center both
    2. H = P^T Q ; U, S, Vt = svd(H)
    3. d = sign(det(V U^T)) -- the reflection fix; without it you will happily
       superpose a structure onto its mirror image
    4. rotate and take the RMSD

    Return a float. Raise ValueError if the shapes disagree.
    """
    raise NotImplementedError("TODO")
