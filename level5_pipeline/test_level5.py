"""Tests for Level 5.

Unlike Levels 1-4 there is no starter.py to fill in: the pipeline is given to
you working. These tests pin the PROPERTIES a campaign must have. Your job in
this level is to extend the pipeline (new stages, new backends, a real model)
and keep every one of them green -- that is what it means for the extension
not to have broken the science.

    python -m pytest level5_pipeline/test_level5.py -q
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

from backends import (Backbone, Design, MockBackboneGenerator,  # noqa: E402
                      MockInterfaceScorer, MockSequenceDesigner,
                      MockStructurePredictor, make_demo_target)
from pipeline import (Campaign, CampaignConfig, Filters,  # noqa: E402
                      build_default_campaign)

SMALL = CampaignConfig(n_backbones=40, seqs_per_backbone=2, top_k=8)


def test_target_is_well_formed():
    t = make_demo_target()
    assert t.coords.shape == (60, 3)
    assert len(t.sequence) == 60
    assert all(0 <= h < 60 for h in t.hotspots)


def test_backbones_are_protein_like():
    t = make_demo_target()
    bbs = MockBackboneGenerator().generate(t, 25, seed=0)
    assert len(bbs) == 25
    assert len({b.id for b in bbs}) == 25
    for b in bbs:
        assert 14 <= b.length <= 34
        step = np.linalg.norm(np.diff(b.coords, axis=0), axis=1)
        assert np.allclose(step, step[0], atol=0.2), \
            "consecutive CA atoms sit at a near-constant ~3.8 A spacing"
        assert 2.0 < step.mean() < 5.0


def test_designed_sequences_are_valid():
    t = make_demo_target()
    b = MockBackboneGenerator().generate(t, 1, seed=1)[0]
    seqs = MockSequenceDesigner().design(b, t, 5, temperature=0.1, seed=0)
    assert len(seqs) == 5
    for s in seqs:
        assert len(s) == b.length
        assert set(s) <= set("ACDEFGHIKLMNPQRSTVWY")


def test_temperature_increases_diversity():
    t = make_demo_target()
    b = MockBackboneGenerator().generate(t, 1, seed=2)[0]

    def diversity(temp):
        seqs = MockSequenceDesigner().design(b, t, 12, temperature=temp, seed=0)
        pairs = [(i, j) for i in range(len(seqs)) for j in range(i + 1, len(seqs))]
        return np.mean([sum(a != c for a, c in zip(seqs[i], seqs[j])) / len(seqs[i])
                        for i, j in pairs])

    assert diversity(0.9) > diversity(0.02), \
        "higher sampling temperature must give more diverse sequences"


def test_scores_are_finite_and_signed_correctly():
    t = make_demo_target()
    b = MockBackboneGenerator().generate(t, 1, seed=3)[0]
    s = MockSequenceDesigner().design(b, t, 1, seed=0)[0]
    d = Design(id="d0", backbone_id=b.id, sequence=s)
    m = {}
    m.update(MockStructurePredictor().predict(d, b, t))
    m.update(MockInterfaceScorer().score(d, b, t))
    assert 0 < m["plddt"] <= 100
    assert 1.0 <= m["i_pae"] <= 32.0
    assert m["n_contacts"] >= 0
    assert m["mock_ddg"] <= 0, "a favourable ddG is negative; keep the sign honest"
    assert all(np.isfinite(v) for v in m.values())


def test_funnel_is_monotonically_non_increasing():
    camp = build_default_campaign(SMALL)
    camp.run(verbose=False)
    counts = [s["n"] for s in camp.funnel]
    design_stage = next(i for i, s in enumerate(camp.funnel)
                        if "sequences designed" in s["stage"])
    # the design stage multiplies candidates; everything after only removes
    after = counts[design_stage:]
    assert after == sorted(after, reverse=True), \
        f"a funnel never grows after fan-out: {after}"


def test_run_is_reproducible():
    a, ra = build_default_campaign(SMALL).run(verbose=False)
    b, rb = build_default_campaign(SMALL).run(verbose=False)
    assert list(a["design_id"]) == list(b["design_id"])
    assert list(a["sequence"]) == list(b["sequence"])
    assert list(ra["design_id"]) == list(rb["design_id"])


def test_different_seed_gives_different_designs():
    cfg2 = CampaignConfig(n_backbones=40, seqs_per_backbone=2, top_k=8, seed=7)
    _a, ra = build_default_campaign(SMALL).run(verbose=False)
    _b, rb = build_default_campaign(cfg2).run(verbose=False)
    assert list(ra["sequence"]) != list(rb["sequence"])


def test_filters_actually_filter():
    camp = build_default_campaign(SMALL)
    df, ranked = camp.run(verbose=False)
    f = camp.cfg.filters
    assert df["passes"].any(), "the default thresholds should let something through"
    assert not df["passes"].all(), "...but not everything"
    kept = df[df["passes"]]
    assert (kept["plddt"] >= f.min_plddt).all()
    assert (kept["i_pae"] <= f.max_i_pae).all()
    assert (kept["longest_agg_run"] <= f.max_agg_run).all()
    assert set(ranked["design_id"]) <= set(kept["design_id"])


def test_impossible_filters_return_empty_not_crash():
    cfg = CampaignConfig(n_backbones=20, seqs_per_backbone=2,
                         filters=Filters(min_plddt=99.9, max_i_pae=0.5))
    df, ranked = build_default_campaign(cfg).run(verbose=False)
    assert len(ranked) == 0, "no survivors is a result, not an exception"
    assert len(df) > 0


def test_diversity_cap_is_respected():
    cfg = CampaignConfig(n_backbones=60, seqs_per_backbone=4, top_k=20,
                         max_per_backbone_in_top=2)
    _df, ranked = build_default_campaign(cfg).run(verbose=False)
    if len(ranked):
        assert ranked["backbone_id"].value_counts().max() <= 2


def test_ranking_prefers_better_designs():
    _df, ranked = build_default_campaign(SMALL).run(verbose=False)
    if len(ranked) > 2:
        assert ranked["composite"].is_monotonic_decreasing


def test_manifest_records_provenance(tmp_path):
    camp = build_default_campaign(SMALL)
    df, ranked = camp.run(verbose=False)
    out = camp.save(df, ranked, str(tmp_path / "run"))
    man = json.load(open(os.path.join(out, "manifest.json")))
    assert man["config"]["seed"] == SMALL.seed
    assert man["backends"]["generator"] == "mock-diffusion"
    assert man["funnel"] and man["n_designs"] == len(df)
    assert os.path.exists(os.path.join(out, "shortlist.fasta"))
    assert os.path.exists(os.path.join(out, "all_designs.csv"))


def test_a_custom_backend_plugs_in_without_touching_the_pipeline():
    """The point of the level: swapping a backend is a constructor argument.

    This fake designer returns poly-alanine. If the campaign still runs, your
    RFdiffusion / ProteinMPNN wrapper will drop in the same way."""
    class PolyAlanineDesigner(MockSequenceDesigner):
        name = "poly-A"

        def design(self, backbone, target, n, temperature=0.1, seed=0):
            return ["A" * backbone.length] * n

    t = make_demo_target()
    camp = Campaign(t, MockBackboneGenerator(), PolyAlanineDesigner(),
                    MockStructurePredictor(), MockInterfaceScorer(),
                    CampaignConfig(n_backbones=15, seqs_per_backbone=1))
    df, _ranked = camp.run(verbose=False)
    assert set(df["sequence"].str.replace("A", "")) == {""}
    assert camp.designer.name == "poly-A"
