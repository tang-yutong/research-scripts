"""Publication figures from the parsed trajectory data.

One axis per configured pair; variant 1 in blue, variant 2 in red. Figures:
C-alpha/ligand RMSD, chain-aware RMSF (continuous index on the axis, original
residue numbers as ticks, chain boundaries from the data itself), radius of
gyration, SASA, hydrogen-bond counts, a free-energy landscape (2D histogram
of C-alpha RMSD x Rg, dG = -kT ln P), MM-PBSA term bars, and per-residue
decompositions. Also writes the equilibration-window statistics table.

Inputs: out/xvg/*.xvg (gmx_batch_analyses), out/mmpbsa_frames.json +
out/mmpbsa_perresidue.json (mmpbsa_parse). Exit codes: 0 ok | 3 missing input.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from config import CONFIG  # noqa: E402
import common  # noqa: E402
from mmpbsa_parse import mean_std  # noqa: E402

matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["svg.fonttype"] = "path"
matplotlib.rcParams["font.family"] = "serif"
matplotlib.rcParams["font.serif"] = ["Times New Roman"]

C_VAR1 = "#1f77b4"
C_VAR2 = "#d62728"
LW = 0.8
FS = 10
KJ2KCAL = 1 / 4.184


def read_xvg(path: Path) -> np.ndarray:
    """Read a GROMACS .xvg (skipping #/@ comment lines) into a float array."""
    rows = []
    for l in open(path, errors="replace"):
        if l.startswith(("#", "@")):
            continue
        p = l.split()
        if p:
            rows.append([float(x) for x in p])
    return np.array(rows)


def window_stats(d: np.ndarray, window_ns: float, conv: float = 1.0):
    """Mean +- sample std of column 1 over the last `window_ns` of column 0
    (the window end is taken from the data, not assumed)."""
    m = d[:, 0] >= (d[-1, 0] - window_ns)
    return mean_std([v * conv for v in d[m, 1]])


class Figures:
    def __init__(self, cfg):
        self.cfg = cfg
        self.pairs = [(f"pair_{i}", v1, v2) for i, (v1, v2) in enumerate(cfg.md_pairs, 1)]

    def load(self, tag, which):
        return read_xvg(self.cfg.xvg_dir / f"{tag}_{which}.xvg")

    def save(self, fig, name):
        for ext in ("png", "pdf", "svg"):
            fig.savefig(self.cfg.figure_dir / f"{name}.{ext}", dpi=600)
        plt.close(fig)
        print("saved:", name)

    def _style_axis(self, ax, pair, title, xlab, ylab):
        ax.set_title(title, fontsize=FS + 1)
        ax.set_xlabel(xlab, fontsize=FS)
        ax.set_ylabel(ylab, fontsize=FS)
        ax.tick_params(labelsize=FS - 1)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    def rmsf_figure(self):
        """RMSF on a continuous residue index (avoids broken lines where a
        structure renumbers mid-chain or concatenates chains); original
        residue numbers as tick labels, chain boundaries drawn from the
        data with configured segment names."""
        fig, axes = plt.subplots(1, 3, figsize=(16 / 2.54 * 2.1, 5.2))
        for ax, (pair, v1, v2) in zip(axes, self.pairs):
            d1 = self.load(v1, "rmsf")
            d2 = self.load(v2, "rmsf")
            x = np.arange(1, len(d1) + 1)
            ax.plot(x, d1[:, 1] * 10, color=C_VAR1, lw=LW, label="var_1")
            ax.plot(x, d2[:, 1] * 10, color=C_VAR2, lw=LW, label="var_2")
            orig = d1[:, 0]
            bnds = list(np.where(np.diff(orig) != 1)[0])
            starts = [0] + [b + 1 for b in bnds] + [len(d1)]
            labels = self.cfg.segment_labels.get(v1, [])
            for b in bnds:
                ax.axvline(b + 1.5, color="0.45", lw=0.9, ls=(0, (4, 3)),
                           alpha=0.9, zorder=3)
            for i in range(len(starts) - 1):
                if i < len(labels):
                    mid = (starts[i] + starts[i + 1]) / 2 + 0.5
                    ax.text(mid, 0.985, labels[i], fontsize=FS - 2.5, color="0.35",
                            ha="center", va="top",
                            transform=ax.get_xaxis_transform())
            ticks = np.linspace(0, len(d1) - 1, 6).astype(int)
            ax.set_xticks(ticks + 1)
            ax.set_xticklabels([str(int(orig[i])) for i in ticks], fontsize=FS - 2)
            self._style_axis(ax, pair, f"{pair} ({v1} / {v2})",
                             "Residue number", "RMSF (Å)")
        axes[0].legend(fontsize=FS - 1, frameon=False)
        fig.tight_layout()
        self.save(fig, "fig_rmsf")

    def trio_figure(self, which, ylab, conv, fname, xl="Time (ns)"):
        fig, axes = plt.subplots(1, 3, figsize=(16 / 2.54 * 1.9, 5.2))
        for ax, (pair, v1, v2) in zip(axes, self.pairs):
            for tag, color, label in [(v1, C_VAR1, "var_1"), (v2, C_VAR2, "var_2")]:
                d = self.load(tag, which)
                ax.plot(d[:, 0], conv(d[:, 1]), color=color, lw=LW, label=label)
            self._style_axis(ax, pair, f"{pair} ({v1} / {v2})", xl, ylab)
        axes[0].legend(fontsize=FS - 1, frameon=False)
        fig.tight_layout()
        self.save(fig, fname)

    def stats_table(self):
        """Mean +- std over the equilibration window (end taken from data)."""
        out = {}
        eq = self.cfg.equilibration_window_ns
        for _pair, v1, v2 in self.pairs:
            for tag in (v1, v2):
                s = {}
                for which, conv in [("rmsd_ca", 10), ("rmsd_ligand", 10),
                                    ("gyrate", 10), ("sasa", 1), ("hbond", 1)]:
                    s[which] = window_stats(self.load(tag, which), eq, conv)
                out[tag] = s
        with open(self.cfg.table_dir / "equilibration_stats.csv", "w",
                  encoding="utf-8-sig", newline="") as f:
            w = csv.writer(f)
            w.writerow(["system", "variant", f"ca_rmsd_A (last {eq:.0f} ns)",
                        "ligand_rmsd_A", "rg_A", "sasa_nm2", "hbonds"])
            for _pair, v1, v2 in self.pairs:
                for tag, variant in ((v1, "var_1"), (v2, "var_2")):
                    s = out[tag]
                    w.writerow([tag, variant] +
                               [f"{s[k][0]:.2f}+-{s[k][1]:.2f}"
                                for k in ("rmsd_ca", "rmsd_ligand", "gyrate")] +
                               [f"{s['sasa'][0]:.1f}+-{s['sasa'][1]:.1f}",
                                f"{s['hbond'][0]:.1f}+-{s['hbond'][1]:.1f}"])
        (self.cfg.out_dir / "stats_eq.json").write_text(
            json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")

    def fel(self):
        """Free-energy landscape: 2D histogram of C-alpha RMSD x Rg,
        dG = -kT ln P, both variants overlaid as contour sets."""
        kt_kj = 8.314462618e-3 * self.cfg.temperature_k
        fig, axes = plt.subplots(1, 3, figsize=(16 / 2.54 * 1.9, 5.2))
        for ax, (pair, v1, v2) in zip(axes, self.pairs):
            data = {}
            for tag in (v1, v2):
                r1 = self.load(tag, "rmsd_ca")
                r2 = self.load(tag, "gyrate")
                n = min(len(r1), len(r2))
                data[tag] = (r1[:n, 1] * 10, r2[:n, 1] * 10)  # Å
            xlo = min(np.min(v[0]) for v in data.values())
            xhi = max(np.max(v[0]) for v in data.values())
            ylo = min(np.min(v[1]) for v in data.values())
            yhi = max(np.max(v[1]) for v in data.values())
            bx = np.linspace(xlo, xhi, self.cfg.fel_bins)
            by = np.linspace(ylo, yhi, self.cfg.fel_bins)
            levels = np.arange(0, 10.5, 1.0)
            for tag, color in [(v1, C_VAR1), (v2, C_VAR2)]:
                x, y = data[tag]
                h, xe, ye = np.histogram2d(x, y, bins=[bx, by])
                g = -kt_kj * np.log(np.maximum(h.T, 1e-12) / h.max())
                g = np.where(h.T > 0, g, np.nan)
                X, Y = np.meshgrid((xe[:-1] + xe[1:]) / 2, (ye[:-1] + ye[1:]) / 2)
                ax.contourf(X, Y, g, levels=levels, cmap=None, colors=[color],
                            alpha=0.13)
                ax.contour(X, Y, g, levels=levels, colors=color, linewidths=0.7)
            self._style_axis(ax, pair, f"{pair} ({v1} / {v2})",
                             "C-alpha RMSD (Å)", "Rg (Å)")
        from matplotlib.lines import Line2D
        handles = [Line2D([0], [0], color=C_VAR1, lw=1.5, label="var_1"),
                   Line2D([0], [0], color=C_VAR2, lw=1.5, label="var_2")]
        axes[0].legend(handles=handles, fontsize=FS - 1, frameon=False)
        fig.tight_layout()
        self.save(fig, "fig_free_energy_landscape")

    def mmpbsa_figs(self):
        """MM-PBSA term bars and per-residue decomposition bars."""
        frames = json.loads((self.cfg.out_dir / "mmpbsa_frames.json").read_text(
            encoding="utf-8"))
        perres = json.loads((self.cfg.out_dir / "mmpbsa_perresidue.json").read_text(
            encoding="utf-8"))
        comps = ["VDW", "COU", "PB", "SA"]
        comp_labels = [r"$\Delta E_{vdW}$", r"$\Delta E_{ele}$",
                       r"$\Delta G_{polar}$", r"$\Delta G_{nonpolar}$"]
        fig, axes = plt.subplots(1, 3, figsize=(16 / 2.54 * 1.9, 5.2))
        for ax, (pair, v1, v2) in zip(axes, self.pairs):
            wid = 0.38
            for k, (tag, color, label) in enumerate([(v1, C_VAR1, "var_1"),
                                                     (v2, C_VAR2, "var_2")]):
                fr = frames["frames"][tag]
                means, stds = [], []
                for c in comps:
                    m, s = mean_std([f[c] * KJ2KCAL for f in fr])
                    means.append(m)
                    stds.append(s)
                su = frames["summary"][tag]
                tds = su["TdS"][0] * KJ2KCAL
                dh_m, dh_s = mean_std([f["dH"] for f in fr])
                means += [tds, (dh_m + tds)]
                stds += [0.0, dh_s * KJ2KCAL]
                ax.bar(np.arange(len(means)) + (k - 0.5) * wid, means, wid,
                       yerr=stds, capsize=2, color=color, label=label,
                       error_kw=dict(lw=0.7))
            ax.axhline(0, color="black", lw=0.7)
            ax.set_xticks(np.arange(len(comps) + 2))
            ax.set_xticklabels(comp_labels + [r"$-T\Delta S$", r"$\Delta G_{bind}$"],
                               fontsize=FS - 1)
            self._style_axis(ax, pair, f"{pair} ({v1} / {v2})", "",
                             "Energy (kcal/mol)")
        axes[0].legend(fontsize=FS - 1, frameon=False)
        fig.tight_layout()
        self.save(fig, "fig_mmpbsa_terms")

        fig, axes = plt.subplots(1, 3, figsize=(16 / 2.54 * 2.1, 5.2))
        for ax, (pair, v1, v2) in zip(axes, self.pairs):
            top1 = [t for t, _ in perres[v1]["top_kJ"][:8]]
            top2 = [t for t, _ in perres[v2]["top_kJ"][:8]]
            ordered = (top1 + [t for t in top2 if t not in top1])[:12]
            d1 = dict(perres[v1]["top_kJ"])
            d2 = dict(perres[v2]["top_kJ"])
            xs = np.arange(len(ordered))
            wid = 0.38
            ax.bar(xs - wid / 2, [d1.get(t, 0) * KJ2KCAL for t in ordered], wid,
                   color=C_VAR1, label="var_1")
            ax.bar(xs + wid / 2, [d2.get(t, 0) * KJ2KCAL for t in ordered], wid,
                   color=C_VAR2, label="var_2")
            ax.axhline(0, color="black", lw=0.7)
            labels = [t[2:] if t.startswith("P~") else t for t in ordered]
            ax.set_xticks(xs)
            ax.set_xticklabels(labels, rotation=60, fontsize=FS - 2.5)
            self._style_axis(ax, pair, f"{pair} ({v1} / {v2})", "",
                             r"$\Delta G_{residue}$ (kcal/mol)")
        axes[0].legend(fontsize=FS - 1, frameon=False)
        fig.tight_layout()
        self.save(fig, "fig_per_residue")

    def run_all(self):
        self.trio_figure("rmsd_ca", "C-alpha RMSD (Å)", lambda v: v * 10,
                         "fig_ca_rmsd")
        self.trio_figure("rmsd_ligand", "Ligand RMSD (Å)", lambda v: v * 10,
                         "fig_ligand_rmsd")
        self.rmsf_figure()
        self.trio_figure("gyrate", "Rg (Å)", lambda v: v * 10, "fig_rg")
        self.trio_figure("sasa", "SASA (nm$^2$)", lambda v: v, "fig_sasa")
        self.trio_figure("hbond", "Number of H-bonds", lambda v: v, "fig_hbonds")
        self.fel()
        self.mmpbsa_figs()
        self.stats_table()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Trajectory figures from parsed data.")
    common.add_common_args(ap)
    args = ap.parse_args(argv)
    log = common.get_logger("trajectory_plots", args.log_level)
    cfg = CONFIG
    if args.config:
        common.apply_overrides(cfg, args.config)
    common.resolve(cfg)

    pairs = [s for pair in cfg.md_pairs for s in pair]
    xvg_inputs = [cfg.xvg_dir / f"{tag}_{w}.xvg"
                  for tag in pairs for w in ("rmsd_ca", "rmsd_ligand", "rmsf",
                                             "gyrate", "sasa", "hbond")]
    json_inputs = [cfg.out_dir / "mmpbsa_frames.json",
                   cfg.out_dir / "mmpbsa_perresidue.json"]
    missing = [str(p) for p in xvg_inputs + json_inputs if not p.exists()]
    if missing:
        # the MM-PBSA figures need the json inputs; xvg inputs are mandatory
        missing_xvg = [m for m in missing if m.endswith(".xvg")]
        missing_json = [m for m in missing if m.endswith(".json")]
        if missing_xvg or len(missing_json) == len(json_inputs):
            log.error("missing input(s): %s", ", ".join(missing))
            return common.EXIT_INPUT
        log.warning("partial inputs (skipping MM-PBSA figures): %s",
                    ", ".join(missing_json))
    if args.dry_run:
        log.info("would write 9 figures x3 formats to %s + equilibration_stats.csv",
                 cfg.figure_dir)
        return common.EXIT_OK

    cfg.figure_dir.mkdir(parents=True, exist_ok=True)
    cfg.table_dir.mkdir(parents=True, exist_ok=True)
    figs = Figures(cfg)
    figs.run_all()
    common.write_manifest(cfg.out_dir / "work" / "trajectory_plots",
                          Path(__file__).name, sys.argv[1:], cfg,
                          inputs=xvg_inputs + json_inputs)
    return common.EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
