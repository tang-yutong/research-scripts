# -*- coding: utf-8 -*-
"""Shared figure toolkit: fonts, export, and RDKit structure rendering.

The structure helpers are chemically aware: spiro atoms are ring-membership
>= 2 carbons with degree 4; bridgehead atoms are ring-membership >= 2 atoms
with degree 3. Highlight colors can be attached to detected motifs before
rendering.
"""
import io
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CONFIG  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from PIL import Image  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem.Draw import rdMolDraw2D  # noqa: E402

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["svg.fonttype"] = "path"


def save(fig, name: str, out_dir=None):
    """PNG 600 dpi + SVG, text-as-paths (editable in presentation software)."""
    out_dir = out_dir or (CONFIG.out_dir / "figures")
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(os.path.join(out_dir, name + ".png"), dpi=600,
                bbox_inches="tight", facecolor="white")
    fig.savefig(os.path.join(out_dir, name + ".svg"), bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    print("saved", name)


def mol_img(smiles, hl=None, w=440, h=400):
    """Render a molecule; hl = atom indices, 'bridgeheads' or 'spiro' for
    automatic motif detection, with a highlight color per atom."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("bad SMILES: " + smiles)
    if hl == "bridgeheads":
        hl = bridgehead_atoms(mol)
    elif hl == "spiro":
        hl = spiro_atoms(mol)
    colors = {i: (1.00, 0.60, 0.60) for i in (hl or [])}
    d = rdMolDraw2D.MolDraw2DCairo(w, h)
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol, highlightAtoms=list(hl or []),
                                       highlightAtomColors=colors)
    d.FinishDrawing()
    return Image.open(io.BytesIO(d.GetDrawingText()))


def spiro_atoms(mol):
    """Spiro atoms: share two rings, degree 4."""
    ri = mol.GetRingInfo()
    return [a.GetIdx() for a in mol.GetAtoms()
            if a.GetDegree() == 4 and ri.NumAtomRings(a.GetIdx()) >= 2]


def bridgehead_atoms(mol):
    """Bridgehead atoms: ring-membership >= 2 with degree 3."""
    ri = mol.GetRingInfo()
    return [a.GetIdx() for a in mol.GetAtoms()
            if a.GetDegree() == 3 and ri.NumAtomRings(a.GetIdx()) >= 2]


def rxn_row(ax, y, left, arrow, right, h=2.35):
    """Draw one row of a reaction scheme on `ax`: reagents on the left, an
    arrow with a caption, products on the right.

    left/right = [("smi", smiles, rel_width) | ("txt", text, rel_width)];
    horizontal layout: left 0-4.3, arrow 4.5-6.5, right 6.7-10."""
    import numpy as np
    from matplotlib.patches import FancyArrowPatch

    def place(items, x0, x1):
        ws = [it[2] for it in items]
        total = sum(ws)
        span = x1 - x0
        cx = x0
        for it in items:
            w = span * (it[2] / total)
            if it[0] == "smi":
                img = mol_img(it[1], w=int(w * 130), h=int(h * 100))
                ax.imshow(img, extent=(cx + 0.08, cx + w - 0.08, y, y + h))
            else:
                ax.text(cx + w / 2, y + h / 2, it[1], ha="center", va="center",
                        fontsize=9.3, color="#333",
                        bbox=dict(boxstyle="round,pad=0.30", fc="#f4f4f4",
                                  ec="#aaa", lw=0.8))
            cx += w

    place(left, 0.0, 4.3)
    ax.add_patch(FancyArrowPatch((4.45, y + h / 2), (6.55, y + h / 2),
                                 arrowstyle="-|>", mutation_scale=17, lw=1.7,
                                 color="#1c3a5a"))
    ax.text(5.5, y + h * 0.72, arrow, ha="center", va="center", fontsize=8.6,
            color="#1c3a5a")
    place(right, 6.7, 10.0)
