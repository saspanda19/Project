"""
LEVEL 5 - the campaign.

    python level5_pipeline/pipeline.py
    python level5_pipeline/pipeline.py --n-backbones 300 --seqs-per-backbone 4

A binder-design campaign is a funnel: many cheap candidates in, few expensive
ones out. The engineering that decides whether it works is not the models --
it is the funnel: what you filter on, in what order, with what thresholds,
and whether you can reproduce the run six months later when a reviewer asks.

This module is deliberately backend-agnostic. Swap MockBackboneGenerator for
an RFdiffusion wrapper and nothing else changes. That is the whole design.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from dataclasses import dataclass, field, asdict
from typing import Optional

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

from backends import (  # noqa: E402
    Design, MockBackboneGenerator, MockInterfaceScorer, MockSequenceDesigner,
    MockStructurePredictor, make_demo_target,
)


# ---------------------------------------------------------------- config ---
@dataclass
class Filters:
    """Every threshold in one place, so a run is described by a config, not by
    edits scattered through the code. Thresholds here are MOCK-SCALE numbers;
    the README gives the real-world ones to use with real models."""
    min_plddt: float = 70.0
    max_i_pae: float = 12.0
    min_contacts: int = 40
    max_agg_run: int = 4
    max_abs_net_charge: float = 8.0
    max_iface_hydrophobic: float = 0.70


@dataclass
class CampaignConfig:
    n_backbones: int = 200
    seqs_per_backbone: int = 3
    temperature: float = 0.15
    seed: int = 0
    top_k: int = 24
    max_per_backbone_in_top: int = 2     # diversity guard
    filters: Filters = field(default_factory=Filters)


# -------------------------------------------------------------- campaign ---
class Campaign:
    def __init__(self, target, generator, designer, predictor, scorer,
                 config: Optional[CampaignConfig] = None):
        self.target = target
        self.generator = generator
        self.designer = designer
        self.predictor = predictor
        self.scorer = scorer
        self.cfg = config or CampaignConfig()
        self.funnel = []          # (stage, n_remaining) - the audit trail

    def _record(self, stage, n):
        self.funnel.append({"stage": stage, "n": int(n)})

    def run(self, verbose=True):
        t0 = time.time()
        cfg = self.cfg

        backbones = self.generator.generate(self.target, cfg.n_backbones,
                                            seed=cfg.seed)
        self._record("backbones generated", len(backbones))

        # cheap geometric pre-filter BEFORE the expensive stages. Ordering the
        # funnel by cost-per-candidate is most of the practical speedup.
        backbones = [b for b in backbones if b.hotspot_contacts >= 2]
        self._record("backbones reaching >=2 hotspots", len(backbones))

        designs, by_backbone = [], {}
        for b in backbones:
            by_backbone[b.id] = b
            for j, seq in enumerate(self.designer.design(
                    b, self.target, cfg.seqs_per_backbone,
                    temperature=cfg.temperature, seed=cfg.seed)):
                designs.append(Design(id=f"{b.id}_s{j}", backbone_id=b.id,
                                      sequence=seq))
        self._record("sequences designed", len(designs))

        for d in designs:
            b = by_backbone[d.backbone_id]
            d.metrics.update(self.predictor.predict(d, b, self.target))
            d.metrics.update(self.scorer.score(d, b, self.target))
        self._record("structures predicted + scored", len(designs))

        df = pd.DataFrame([d.row() for d in designs])
        df = self._apply_filters(df)
        ranked = self._rank(df)
        if verbose:
            self.report(df, ranked, time.time() - t0)
        return df, ranked

    def _apply_filters(self, df):
        f = self.cfg.filters
        checks = [
            ("plddt", df["plddt"] >= f.min_plddt),
            ("i_pae", df["i_pae"] <= f.max_i_pae),
            ("contacts", df["n_contacts"] >= f.min_contacts),
            ("aggregation", df["longest_agg_run"] <= f.max_agg_run),
            ("net charge", df["net_charge"].abs() <= f.max_abs_net_charge),
            ("iface hydrophobicity",
             df["iface_hydrophobic_frac"] <= f.max_iface_hydrophobic),
        ]
        keep = pd.Series(True, index=df.index)
        self.filter_losses = []
        for name, ok in checks:
            # attribute each design's FIRST failure, so the report tells you
            # which threshold is actually shaping your output
            newly_lost = int((keep & ~ok).sum())
            self.filter_losses.append((name, newly_lost))
            keep &= ok
        df = df.copy()
        df["passes"] = keep
        self._record("passed all filters", int(keep.sum()))
        return df

    def _rank(self, df):
        """Composite score, then a diversity cap.

        Ranking on a single metric concentrates your shortlist on one backbone
        and one failure mode. The cap is cheap insurance."""
        ok = df[df["passes"]].copy()
        if ok.empty:
            return ok
        z = lambda s: (s - s.mean()) / (s.std() + 1e-9)  # noqa: E731
        ok["composite"] = (1.0 * z(-ok["mock_ddg"]) + 0.8 * z(ok["plddt"])
                           + 0.8 * z(-ok["i_pae"]) + 0.4 * z(ok["n_contacts"]))
        ok = ok.sort_values("composite", ascending=False)
        picked, seen = [], {}
        for _, row in ok.iterrows():
            b = row["backbone_id"]
            if seen.get(b, 0) >= self.cfg.max_per_backbone_in_top:
                continue
            seen[b] = seen.get(b, 0) + 1
            picked.append(row)
            if len(picked) >= self.cfg.top_k:
                break
        out = pd.DataFrame(picked).reset_index(drop=True)
        self._record(f"shortlist (top {self.cfg.top_k}, "
                     f"max {self.cfg.max_per_backbone_in_top}/backbone)",
                     len(out))
        return out

    # ------------------------------------------------------------ output ---
    def report(self, df, ranked, elapsed):
        print(f"\ntarget: {self.target.name}  "
              f"({len(self.target.sequence)} residues, "
              f"{len(self.target.hotspots)} hotspots)")
        print(f"backends: {self.generator.name} -> {self.designer.name} -> "
              f"{self.predictor.name} -> {self.scorer.name}\n")
        print("FUNNEL")
        for step in self.funnel:
            print(f"  {step['stage']:<52}{step['n']:>7}")
        print("\nFIRST FAILING FILTER (of designs still alive at that point)")
        for name, lost in self.filter_losses:
            print(f"  {name:<52}{lost:>7}")
        if ranked.empty:
            print("\nNothing survived. That is a RESULT, not a crash: loosen a "
                  "threshold\nand watch which one was doing the work.")
            return
        cols = ["design_id", "length", "plddt", "i_pae", "n_contacts",
                "mock_ddg", "composite"]
        print(f"\nTOP {min(8, len(ranked))} DESIGNS")
        print(ranked[cols].head(8).to_string(
            index=False, float_format=lambda v: f"{v:7.2f}"))
        print(f"\n  distinct backbones in the shortlist: "
              f"{ranked['backbone_id'].nunique()} / {len(ranked)}")
        print(f"  wall clock: {elapsed:.2f}s")

    def save(self, df, ranked, outdir):
        os.makedirs(outdir, exist_ok=True)
        df.to_csv(os.path.join(outdir, "all_designs.csv"), index=False)
        ranked.to_csv(os.path.join(outdir, "shortlist.csv"), index=False)
        with open(os.path.join(outdir, "shortlist.fasta"), "w") as fh:
            for _, r in ranked.iterrows():
                fh.write(f">{r['design_id']} composite={r['composite']:.3f} "
                         f"plddt={r['plddt']:.1f}\n{r['sequence']}\n")
        # provenance: the run is not reproducible without it
        manifest = {
            "config": asdict(self.cfg),
            "backends": {"generator": self.generator.name,
                         "designer": self.designer.name,
                         "predictor": self.predictor.name,
                         "scorer": self.scorer.name},
            "target": {"name": self.target.name,
                       "length": len(self.target.sequence),
                       "hotspots": self.target.hotspots},
            "funnel": self.funnel,
            "first_failing_filter": dict(self.filter_losses),
            "n_designs": int(len(df)), "n_shortlist": int(len(ranked)),
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "numpy": np.__version__, "python": sys.version.split()[0],
        }
        with open(os.path.join(outdir, "manifest.json"), "w") as fh:
            json.dump(manifest, fh, indent=2)
        return outdir


def plot_funnel(campaign, df, outpath):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return None
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    stages = [s["stage"] for s in campaign.funnel]
    counts = [s["n"] for s in campaign.funnel]
    axes[0].barh(range(len(stages))[::-1], counts, color="#4C6EF5")
    axes[0].set_yticks(range(len(stages))[::-1])
    axes[0].set_yticklabels([s[:38] for s in stages], fontsize=8)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("candidates remaining (log)")
    axes[0].set_title("the funnel")
    passed = df[df["passes"]]
    axes[1].scatter(df["i_pae"], df["plddt"], s=8, alpha=0.25,
                    color="#adb5bd", label="filtered out")
    axes[1].scatter(passed["i_pae"], passed["plddt"], s=14, alpha=0.85,
                    color="#4C6EF5", label="passed")
    axes[1].axhline(campaign.cfg.filters.min_plddt, ls="--", lw=0.8, c="k")
    axes[1].axvline(campaign.cfg.filters.max_i_pae, ls="--", lw=0.8, c="k")
    axes[1].set_xlabel("interface PAE (lower is better)")
    axes[1].set_ylabel("pLDDT")
    axes[1].set_title("where the thresholds cut")
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    return outpath


def build_default_campaign(cfg=None):
    return Campaign(
        target=make_demo_target(),
        generator=MockBackboneGenerator(),
        designer=MockSequenceDesigner(),
        predictor=MockStructurePredictor(),
        scorer=MockInterfaceScorer(),
        config=cfg or CampaignConfig(),
    )


def main():
    ap = argparse.ArgumentParser(description="mock binder-design campaign")
    ap.add_argument("--n-backbones", type=int, default=200)
    ap.add_argument("--seqs-per-backbone", type=int, default=3)
    ap.add_argument("--temperature", type=float, default=0.15)
    ap.add_argument("--top-k", type=int, default=24)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--outdir", default=os.path.join(ROOT, "runs", "demo"))
    a = ap.parse_args()

    cfg = CampaignConfig(n_backbones=a.n_backbones,
                         seqs_per_backbone=a.seqs_per_backbone,
                         temperature=a.temperature, top_k=a.top_k, seed=a.seed)
    camp = build_default_campaign(cfg)
    df, ranked = camp.run()
    camp.save(df, ranked, a.outdir)
    print(f"\nwrote {a.outdir}/{{all_designs.csv,shortlist.csv,"
          f"shortlist.fasta,manifest.json}}")
    fig = plot_funnel(camp, df, os.path.join(ROOT, "figures",
                                             "level5_funnel.png"))
    if fig:
        print(f"wrote {fig}")


if __name__ == "__main__":
    os.makedirs(os.path.join(ROOT, "figures"), exist_ok=True)
    main()
