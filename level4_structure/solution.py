"""LEVEL 4 reference solution.

    python level4_structure/solution.py

Prints chain sequences, contact statistics and the A/B interface, and writes
figures/level4_contacts.png.
"""
import os
import sys

import numpy as np

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

THREE_TO_ONE = _sibling("starter").THREE_TO_ONE


def parse_pdb(path):
    atoms = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            atoms.append({
                "serial": int(line[6:11]),
                "atom": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21],
                "resseq": int(line[22:26]),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
            })
    return atoms


def chain_sequence(atoms, chain):
    seen = {}
    for a in atoms:
        if a["chain"] == chain and a["atom"] == "CA":
            seen.setdefault(a["resseq"], THREE_TO_ONE.get(a["resname"], "X"))
    return "".join(seen[k] for k in sorted(seen))


def coords(atoms, chain=None, atom_name="CA"):
    sel = [a for a in atoms
           if (chain is None or a["chain"] == chain)
           and (atom_name is None or a["atom"] == atom_name)]
    return np.array([[a["x"], a["y"], a["z"]] for a in sel], dtype=float)


def distance_matrix(a, b=None):
    a = np.asarray(a, dtype=float)
    b = a if b is None else np.asarray(b, dtype=float)
    diff = a[:, None, :] - b[None, :, :]
    return np.sqrt((diff ** 2).sum(-1))


def contact_map(atoms, chain, cutoff=8.0, min_seq_sep=3):
    xyz = coords(atoms, chain)
    d = distance_matrix(xyz)
    n = len(xyz)
    i, j = np.indices((n, n))
    return (d < cutoff) & (np.abs(i - j) >= min_seq_sep)


def interface_residues(atoms, chain_a="A", chain_b="B", cutoff=8.0):
    a_atoms = [x for x in atoms if x["chain"] == chain_a and x["atom"] == "CA"]
    b_atoms = [x for x in atoms if x["chain"] == chain_b and x["atom"] == "CA"]
    d = distance_matrix(coords(atoms, chain_a), coords(atoms, chain_b))
    hit = d < cutoff
    res_a = sorted({a_atoms[i]["resseq"] for i in np.where(hit.any(1))[0]})
    res_b = sorted({b_atoms[j]["resseq"] for j in np.where(hit.any(0))[0]})
    return res_a, res_b


def radius_of_gyration(xyz):
    xyz = np.asarray(xyz, dtype=float)
    center = xyz.mean(0)
    return float(np.sqrt(((xyz - center) ** 2).sum(1).mean()))


def kabsch_rmsd(p, q):
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    if p.shape != q.shape:
        raise ValueError(f"shape mismatch: {p.shape} vs {q.shape}")
    pc = p - p.mean(0)
    qc = q - q.mean(0)
    H = pc.T @ qc
    U, _S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    rotated = (R @ pc.T).T
    return float(np.sqrt(((rotated - qc) ** 2).sum(1).mean()))


# ------------------------------------------------------------------ demo ---
def main():
    pdb = os.path.join(ROOT, "data", "bundle_peptide.pdb")
    atoms = parse_pdb(pdb)
    print(f"{len(atoms)} atoms, chains "
          f"{sorted({a['chain'] for a in atoms})}\n")

    for ch in ("A", "B"):
        seq = chain_sequence(atoms, ch)
        xyz = coords(atoms, ch)
        print(f"chain {ch}: {len(seq)} residues  Rg={radius_of_gyration(xyz):5.2f} A")
        print(f"          {seq}")

    cm = contact_map(atoms, "A", cutoff=8.0, min_seq_sep=3)
    n = cm.shape[0]
    print(f"\nchain A contact map: {int(cm.sum()) // 2} contacts among {n} residues")
    # long-range contacts (|i-j| > 12) are the ones that report on the fold
    i, j = np.indices((n, n))
    lr = int((cm & (np.abs(i - j) > 12)).sum()) // 2
    print(f"  of which long-range (|i-j| > 12): {lr}"
          "  <- these are what say 'bundle' rather than 'single helix'")

    ra, rb = interface_residues(atoms, "A", "B", cutoff=10.0)
    print(f"\ninterface at 10 A: {len(ra)} residues in A, {len(rb)} in B")
    print(f"  A: {ra}")
    print(f"  B: {rb}")

    # Kabsch sanity check: rotate chain B and recover ~0 RMSD
    xyz = coords(atoms, "B")
    theta = 0.7
    R = np.array([[np.cos(theta), -np.sin(theta), 0],
                  [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    moved = (R @ xyz.T).T + np.array([10.0, -5.0, 3.0])
    print(f"\nKabsch RMSD after a rigid rotation+translation: "
          f"{kabsch_rmsd(xyz, moved):.2e} A  (should be ~0)")
    noisy = xyz + np.random.default_rng(0).normal(0, 0.5, xyz.shape)
    print(f"Kabsch RMSD with 0.5 A gaussian noise:            "
          f"{kabsch_rmsd(xyz, noisy):.3f} A")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    figdir = os.path.join(ROOT, "figures")
    os.makedirs(figdir, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))
    axes[0].imshow(cm, cmap="Greys", origin="lower", interpolation="nearest")
    axes[0].set_title("chain A contact map (CA-CA < 8 A, |i-j| >= 3)")
    axes[0].set_xlabel("residue i")
    axes[0].set_ylabel("residue j")
    dab = distance_matrix(coords(atoms, "A"), coords(atoms, "B"))
    im = axes[1].imshow(dab, cmap="viridis_r", origin="lower", aspect="auto")
    axes[1].set_title("A-B distance matrix (A)")
    axes[1].set_xlabel("chain B residue")
    axes[1].set_ylabel("chain A residue")
    fig.colorbar(im, ax=axes[1], label="distance (A)")
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "level4_contacts.png"), dpi=150)
    print(f"\nwrote {os.path.join(figdir, 'level4_contacts.png')}")


if __name__ == "__main__":
    main()
