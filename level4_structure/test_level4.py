"""Tests for Level 4.  LADDER_MOD=solution to check the reference."""
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

PDB = os.path.join(ROOT, "data", "bundle_peptide.pdb")


def test_parse_pdb():
    atoms = m.parse_pdb(PDB)
    assert len(atoms) == 94, "82 in chain A + 12 in chain B; TER/REMARK are not ATOM"
    a = atoms[0]
    assert set(a) >= {"serial", "atom", "resname", "chain", "resseq",
                      "x", "y", "z"}
    assert a["atom"] == "CA" and a["chain"] == "A" and a["resseq"] == 1
    assert isinstance(a["x"], float) and isinstance(a["resseq"], int)
    assert {x["chain"] for x in atoms} == {"A", "B"}
    assert all(len(x["resname"]) == 3 for x in atoms)


def test_parse_pdb_is_column_based():
    """Real PDB coordinate columns can run together with no space between
    them. A whitespace split silently mangles those lines; column slicing
    does not."""
    line = ("ATOM    999  CA  TRP B 123    "
            "-123.456-112.345  99.999  1.00 20.00           C")
    tmp = os.path.join(os.path.dirname(PDB), "_tmp_colcheck.pdb")
    with open(tmp, "w") as fh:
        fh.write(line + "\nEND\n")
    try:
        a = m.parse_pdb(tmp)[0]
        assert abs(a["x"] - -123.456) < 1e-6
        assert abs(a["y"] - -112.345) < 1e-6
        assert abs(a["z"] - 99.999) < 1e-6
        assert a["chain"] == "B" and a["resseq"] == 123 and a["resname"] == "TRP"
    finally:
        os.remove(tmp)


def test_chain_sequence():
    atoms = m.parse_pdb(PDB)
    a, b = m.chain_sequence(atoms, "A"), m.chain_sequence(atoms, "B")
    assert len(a) == 82 and len(b) == 12
    assert b == "RGDWYKPTAEIL"
    assert set(a) <= set("ACDEFGHIKLMNPQRSTVWYX")
    assert "GSGSG" in a, "the two helix linkers should survive parsing"


def test_coords():
    atoms = m.parse_pdb(PDB)
    assert m.coords(atoms).shape == (94, 3)
    assert m.coords(atoms, "A").shape == (82, 3)
    assert m.coords(atoms, "B").shape == (12, 3)


def test_distance_matrix():
    a = np.array([[0.0, 0, 0], [3.0, 4, 0]])
    d = m.distance_matrix(a)
    assert d.shape == (2, 2)
    assert abs(d[0, 0]) < 1e-12 and abs(d[0, 1] - 5.0) < 1e-12
    assert np.allclose(d, d.T)
    b = np.array([[0.0, 0, 1]])
    assert m.distance_matrix(a, b).shape == (2, 1)


def test_contact_map():
    atoms = m.parse_pdb(PDB)
    cm = m.contact_map(atoms, "A", cutoff=8.0, min_seq_sep=3)
    n = len(m.chain_sequence(atoms, "A"))
    assert cm.shape == (n, n)
    assert cm.dtype == bool
    assert np.array_equal(cm, cm.T), "a contact map is symmetric"
    i, j = np.indices((n, n))
    assert not cm[np.abs(i - j) < 3].any(), "near-diagonal band must be zeroed"
    assert cm.sum() > 0
    # an ideal alpha helix puts i and i+4 about 6.2 A apart -> in contact
    assert cm[10, 14]
    # a bundle has long-range contacts; a single helix would have none
    assert (cm & (np.abs(i - j) > 12)).sum() > 0


def test_contact_map_cutoff_is_monotone():
    atoms = m.parse_pdb(PDB)
    tight = m.contact_map(atoms, "A", cutoff=6.0).sum()
    loose = m.contact_map(atoms, "A", cutoff=12.0).sum()
    assert loose > tight


def test_interface_residues():
    atoms = m.parse_pdb(PDB)
    ra, rb = m.interface_residues(atoms, "A", "B", cutoff=10.0)
    assert ra and rb
    assert set(ra) <= set(range(1, 83)) and set(rb) <= set(range(1, 13))
    assert ra == sorted(ra) and rb == sorted(rb)
    assert max(ra) < 30, "the peptide docks against helix 1, not the whole bundle"
    tight, _ = m.interface_residues(atoms, "A", "B", cutoff=6.0)
    assert len(tight) <= len(ra)


def test_radius_of_gyration():
    # 8 corners of a cube of side 2 centered at the origin: Rg = sqrt(3)
    pts = np.array([[x, y, z] for x in (-1, 1) for y in (-1, 1)
                    for z in (-1, 1)], dtype=float)
    assert abs(m.radius_of_gyration(pts) - np.sqrt(3)) < 1e-9
    # translation invariant
    assert abs(m.radius_of_gyration(pts + 100) -
               m.radius_of_gyration(pts)) < 1e-9
    atoms = m.parse_pdb(PDB)
    assert m.radius_of_gyration(m.coords(atoms, "A")) > \
        m.radius_of_gyration(m.coords(atoms, "B"))


def test_kabsch_rmsd():
    rng = np.random.default_rng(0)
    p = rng.normal(size=(20, 3))
    assert m.kabsch_rmsd(p, p) < 1e-9
    theta = 0.9
    R = np.array([[np.cos(theta), -np.sin(theta), 0],
                  [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    moved = (R @ p.T).T + np.array([5.0, -2.0, 7.0])
    assert m.kabsch_rmsd(p, moved) < 1e-9, "RMSD must be rotation+translation invariant"
    noisy = p + rng.normal(0, 0.1, p.shape)
    assert 0.0 < m.kabsch_rmsd(p, noisy) < 0.3


def test_kabsch_rejects_mirror_images():
    """Without the reflection fix (the sign of the determinant), a mirrored
    structure superposes onto the original with RMSD ~0. It should not."""
    rng = np.random.default_rng(1)
    p = rng.normal(size=(30, 3))
    mirrored = p * np.array([1.0, 1.0, -1.0])
    assert m.kabsch_rmsd(p, mirrored) > 0.5


def test_kabsch_shape_mismatch():
    try:
        m.kabsch_rmsd(np.zeros((5, 3)), np.zeros((6, 3)))
    except ValueError:
        return
    raise AssertionError("mismatched shapes should raise ValueError")
