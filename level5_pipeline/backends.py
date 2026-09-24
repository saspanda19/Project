"""
LEVEL 5 - the backends.

Every stage of a binder-design campaign is one of four jobs:

    BackboneGenerator   target + hotspots      -> candidate backbones
    SequenceDesigner    backbone               -> sequences for that backbone
    StructurePredictor  sequence (+ target)    -> predicted complex + confidence
    Scorer              predicted complex      -> numbers you can filter on

RFdiffusion, ProteinMPNN, AlphaFold/Boltz/Chai and Rosetta each fill exactly
one of those slots. This file defines the four interfaces and provides MOCK
implementations that run in milliseconds on a CPU, so you can build, test and
debug the *pipeline* -- which is where campaigns actually go wrong -- without
a GPU.

The mocks are deliberately crude but not arbitrary: they encode a few real
regularities (hydrophobic/charged complementarity at an interface, a length
sweet spot, aggregation-prone stretches being bad) so that ranking, filtering
and the score-vs-diversity trade-off behave qualitatively like the real thing.

    A mock model is not a small real model. It tells you your plumbing works.
    It tells you NOTHING about whether a design binds.
"""
from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass, field, asdict
from typing import List, Optional

import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"
HYDROPHOBIC = set("AVILMFWY")
CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1}
# a crude "sticky, aggregation-prone" alphabet
AGG = set("VIYFWL")


# --------------------------------------------------------------- records ---
@dataclass
class Backbone:
    """A candidate binder backbone (CA-only here)."""
    id: str
    coords: np.ndarray                 # (L, 3)
    hotspot_contacts: int              # how many target hotspots it reaches
    parent: Optional[str] = None       # for partial-diffusion style lineages

    @property
    def length(self) -> int:
        return len(self.coords)


@dataclass
class Design:
    """One sequence threaded onto one backbone, plus everything learned later."""
    id: str
    backbone_id: str
    sequence: str
    metrics: dict = field(default_factory=dict)

    def row(self) -> dict:
        out = {"design_id": self.id, "backbone_id": self.backbone_id,
               "sequence": self.sequence, "length": len(self.sequence)}
        out.update(self.metrics)
        return out


@dataclass
class Target:
    """The thing you are designing against."""
    name: str
    sequence: str
    coords: np.ndarray                 # (N, 3) CA coordinates
    hotspots: List[int]                # 0-based residue indices


def _seeded_rng(*parts) -> np.random.Generator:
    """Deterministic RNG from any hashable parts -- same inputs, same output,
    which is what makes a pipeline reproducible and a bug findable."""
    h = hashlib.sha256("|".join(map(str, parts)).encode()).digest()
    return np.random.default_rng(int.from_bytes(h[:8], "big"))


# -------------------------------------------------------------- backends ---
class BackboneGenerator:
    """Interface. A real implementation shells out to RFdiffusion."""
    name = "abstract"

    def generate(self, target: Target, n: int, seed: int = 0) -> List[Backbone]:
        raise NotImplementedError


class MockBackboneGenerator(BackboneGenerator):
    """Emits idealized helical backbones placed near the target hotspots.

    'Quality' here = how many hotspots the backbone reaches. Real diffusion
    models produce a similar spread: most samples are mediocre, a few are good.
    """
    name = "mock-diffusion"

    def __init__(self, min_len=14, max_len=34):
        self.min_len, self.max_len = min_len, max_len

    def generate(self, target, n, seed=0):
        out = []
        for i in range(n):
            rng = _seeded_rng("bb", target.name, seed, i)
            L = int(rng.integers(self.min_len, self.max_len + 1))
            anchor = target.coords[rng.choice(target.hotspots)]
            offset = rng.normal(0, 4.0, 3)
            offset *= 9.0 / (np.linalg.norm(offset) + 1e-9)
            axis = rng.normal(size=3)
            axis /= np.linalg.norm(axis)
            perp = np.cross(axis, rng.normal(size=3))
            perp /= np.linalg.norm(perp) + 1e-9
            perp2 = np.cross(axis, perp)
            coords = np.array([
                anchor + offset + axis * 1.5 * k
                + 2.3 * (math.cos(math.radians(100 * k)) * perp
                         + math.sin(math.radians(100 * k)) * perp2)
                for k in range(L)])
            d = np.linalg.norm(
                coords[:, None, :] - target.coords[None, target.hotspots, :],
                axis=-1)
            hot = int((d.min(0) < 11.0).sum())
            out.append(Backbone(id=f"bb{i:04d}", coords=coords,
                                hotspot_contacts=hot))
        return out


class SequenceDesigner:
    """Interface. A real implementation calls ProteinMPNN / LigandMPNN."""
    name = "abstract"

    def design(self, backbone: Backbone, target: Target, n: int,
               temperature: float = 0.1, seed: int = 0) -> List[str]:
        raise NotImplementedError


class MockSequenceDesigner(SequenceDesigner):
    """Position-aware sampling: buried positions (many intra-backbone
    neighbours) get hydrophobics, exposed positions get polars/charges,
    positions facing target hotspots get complementary charge.

    Temperature does what it does in MPNN: low = conservative and repetitive,
    high = diverse and worse. Watch that trade-off in the demo.
    """
    name = "mock-mpnn"
    BURIED = "AVILMFW"
    EXPOSED = "EKRQNSTDH"

    def design(self, backbone, target, n, temperature=0.1, seed=0):
        d_self = np.linalg.norm(
            backbone.coords[:, None, :] - backbone.coords[None, :, :], axis=-1)
        burial = (d_self < 10.0).sum(1) - 1
        d_tgt = np.linalg.norm(
            backbone.coords[:, None, :] - target.coords[None, :, :], axis=-1)
        facing = d_tgt.min(1) < 12.0
        tgt_charge = np.array([
            sum(CHARGE.get(target.sequence[j], 0)
                for j in np.where(d_tgt[i] < 12.0)[0])
            for i in range(backbone.length)])

        seqs = []
        for k in range(n):
            rng = _seeded_rng("seq", backbone.id, seed, k)
            chars = []
            for i in range(backbone.length):
                if facing[i] and tgt_charge[i] > 0:
                    pool = "DDEESTNQ"
                elif facing[i] and tgt_charge[i] < 0:
                    pool = "KKRRSTNQ"
                elif burial[i] >= np.median(burial):
                    pool = self.BURIED
                else:
                    pool = self.EXPOSED
                if rng.random() < temperature:
                    pool = AA            # temperature = chance of going rogue
                chars.append(pool[int(rng.integers(len(pool)))])
            seqs.append("".join(chars))
        return seqs


class StructurePredictor:
    """Interface. A real implementation calls AlphaFold-Multimer / Boltz /
    Chai and returns pLDDT, PAE and the predicted complex."""
    name = "abstract"

    def predict(self, design: Design, backbone: Backbone,
                target: Target) -> dict:
        raise NotImplementedError


class MockStructurePredictor(StructurePredictor):
    """Returns a pLDDT-like confidence and an interface PAE-like number.

    Confidence rewards: designed-for burial actually being hydrophobic, a
    length near 24, and low sequence entropy. It is a smooth function of the
    sequence, which is exactly what makes it a useful *plumbing* test and a
    useless *biology* test.
    """
    name = "mock-folding"

    def predict(self, design, backbone, target):
        s = design.sequence
        rng = _seeded_rng("fold", design.id)
        hyd = sum(c in HYDROPHOBIC for c in s) / len(s)
        # penalize both extremes: all-hydrophobic aggregates, all-polar unfolds
        hyd_term = 1.0 - abs(hyd - 0.42) / 0.42
        len_term = 1.0 - abs(len(s) - 24) / 24
        counts = np.array([s.count(a) for a in AA], dtype=float)
        p = counts[counts > 0] / len(s)
        entropy = float(-(p * np.log(p)).sum()) / math.log(20)
        plddt = 100 * np.clip(
            0.38 + 0.30 * hyd_term + 0.18 * len_term + 0.12 * entropy
            + rng.normal(0, 0.03), 0.05, 0.98)
        contact_term = min(backbone.hotspot_contacts / 6.0, 1.0)
        i_pae = float(np.clip(
            30 - 20 * contact_term - 6 * (plddt / 100) + rng.normal(0, 1.2),
            1.0, 32.0))
        return {"plddt": float(plddt), "i_pae": i_pae}


class Scorer:
    """Interface. A real implementation runs Rosetta InterfaceAnalyzer, or
    any physics/ML score you trust."""
    name = "abstract"

    def score(self, design: Design, backbone: Backbone, target: Target) -> dict:
        raise NotImplementedError


class MockInterfaceScorer(Scorer):
    """Interface contacts, shape-ish complementarity, charge complementarity,
    and two developability red flags (aggregation stretch, extreme pI proxy).
    """
    name = "mock-interface"

    def score(self, design, backbone, target):
        s = design.sequence
        d = np.linalg.norm(
            backbone.coords[:, None, :] - target.coords[None, :, :], axis=-1)
        contacts = int((d < 10.0).sum())
        iface_idx = np.where(d.min(1) < 10.0)[0]
        iface_res = [s[i] for i in iface_idx]
        hydrophobic_iface = (sum(c in HYDROPHOBIC for c in iface_res) /
                             max(len(iface_res), 1))
        q_binder = sum(CHARGE.get(c, 0) for c in iface_res)
        q_target = sum(CHARGE.get(target.sequence[j], 0)
                       for j in np.where(d.min(0) < 10.0)[0])
        charge_comp = -q_binder * q_target / 10.0   # opposite charges -> +
        longest_agg = _longest_run(s, AGG)
        ddg = -(0.16 * contacts + 6.0 * hydrophobic_iface + 1.2 * charge_comp)
        return {
            "n_contacts": contacts,
            "iface_hydrophobic_frac": float(hydrophobic_iface),
            "charge_complementarity": float(charge_comp),
            "longest_agg_run": longest_agg,
            "mock_ddg": float(ddg),
            "net_charge": float(sum(CHARGE.get(c, 0) for c in s)),
        }


def _longest_run(seq, alphabet):
    best = cur = 0
    for c in seq:
        cur = cur + 1 if c in alphabet else 0
        best = max(best, cur)
    return best


# -------------------------------------------------------------- a target ---
def make_demo_target(seed: int = 0) -> Target:
    """A synthetic 'receptor': a helical hairpin with three hotspot residues."""
    rng = np.random.default_rng(seed)
    coords = []
    for arm in range(2):
        for k in range(30):
            ang = math.radians(100 * k)
            coords.append([2.3 * math.cos(ang) + arm * 11.0,
                           2.3 * math.sin(ang),
                           1.5 * k * (1 if arm == 0 else -1) + arm * 44])
    coords = np.array(coords)
    seq = "".join(rng.choice(list("AEKLQRSTVDNIFM")) for _ in range(len(coords)))
    seq = seq[:20] + "EDDE" + seq[24:]          # an acidic patch to complement
    return Target(name="demo_receptor", sequence=seq, coords=coords,
                  hotspots=[20, 21, 22, 23, 35])


__all__ = [
    "Backbone", "Design", "Target", "BackboneGenerator", "SequenceDesigner",
    "StructurePredictor", "Scorer", "MockBackboneGenerator",
    "MockSequenceDesigner", "MockStructurePredictor", "MockInterfaceScorer",
    "make_demo_target", "asdict",
]
