# -*- coding: utf-8 -*-
"""Worked example: an MD simulation workflow diagram drawn with figkit
conventions (600 dpi PNG + text-as-path SVG). Boxes carry the standard
stage pipeline -- edit the STEPS list for your own protocol."""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common  # noqa: E402
from config import CONFIG  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch  # noqa: E402

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

STEPS = [
    ("System preparation", "structure preparation\nligand topology and charges (AM1-BCC)", "#eaf2fb"),
    ("Solvation and neutralization", "cubic box >= 1.0 nm - SPC water - charge neutral", "#eaf2fb"),
    ("Energy minimization", "steepest descent - Fmax < 10 kJ mol$^{-1}$ nm$^{-1}$", "#eaf2fb"),
    ("NVT equilibration", "v-rescale - 310 K - 500 ps", "#f5f0e6"),
    ("NPT equilibration", "v-rescale + Parrinello-Rahman - 1 bar - 2 ns", "#f5f0e6"),
    ("Production", "Nose-Hoover + Parrinello-Rahman - 310 K / 1 bar\n2 fs - PME - 100 ns", "#f5f0e6"),
    ("Trajectory analysis", "RMSD / RMSF / Rg / SASA / H-bonds / FEL", "#eaf6ee"),
    ("Binding free energy", "MM-PBSA\nper-residue decomposition", "#eaf6ee"),
]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Draw the example MD workflow figure.")
    common.add_common_args(ap)
    args = ap.parse_args(argv)
    log = common.get_logger("example_md_workflow_figure", args.log_level)
    cfg = CONFIG
    if args.config:
        common.apply_overrides(cfg, args.config)
    common.resolve(cfg)

    out_dir = cfg.out_dir / "figures"
    if args.dry_run:
        log.info("would write %s\\md_workflow_figure.{png,svg}", out_dir)
        return common.EXIT_OK
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6.6, 8.2), dpi=100)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, len(STEPS) * 2.6 + 0.4)
    ax.axis("off")

    n = len(STEPS)
    for i, (title, detail, color) in enumerate(STEPS):
        y = (n - i) * 2.6 - 2.2
        box = FancyBboxPatch((1.2, y), 7.6, 1.9,
                             boxstyle="round,pad=0.12,rounding_size=0.18",
                             fc=color, ec="#3a5a78", lw=1.3)
        ax.add_patch(box)
        ax.text(5.0, y + 1.32, title, ha="center", va="center", fontsize=12.5,
                fontweight="bold", color="#1c3a5a")
        ax.text(5.0, y + 0.55, detail, ha="center", va="center", fontsize=9.0,
                color="#333333")
        if i < n - 1:
            arr = FancyArrowPatch((5.0, y - 0.10), (5.0, y - 0.62),
                                  arrowstyle="-|>", mutation_scale=16, lw=1.4,
                                  color="#3a5a78")
            ax.add_patch(arr)

    for ext, dpi in (("png", 600), ("svg", None)):
        p = out_dir / f"md_workflow_figure.{ext}"
        fig.savefig(p, dpi=dpi, bbox_inches="tight", facecolor="white")
    log.info("saved %s", out_dir)
    return common.EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
