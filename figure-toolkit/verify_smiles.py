# -*- coding: utf-8 -*-
"""Gate before any figure: verify every candidate SMILES.

Each structure is parsed with RDKit, its molecular formula recomputed and
compared against the expected value, and a preview grid is rendered for one
visual pass. A wrong structure in a figure is worse than no figure -- this
catches transcription-class errors before they reach one.

Usage: python verify_smiles.py            # verify the example set
       python verify_smiles.py extra.smi  # verify "name SMILES formula" lines
Exit codes: 0 all verified | 3 any verification failed.
"""
import io
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common  # noqa: E402
from config import CONFIG  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem.Draw import rdMolDraw2D  # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from PIL import Image  # noqa: E402

# example set: name -> (SMILES, expected formula or None to skip the check)
CAND = {
    "benzene": ("c1ccccc1", "C6H6"),
    "ethanol": ("CCO", "C2H6O"),
    "aspirin": ("CC(=O)Oc1ccccc1C(=O)O", "C9H8O4"),
}


def verify(cand: dict):
    """-> dict of verified name -> SMILES; prints one line per compound."""
    ok = {}
    for name, (smi, expect) in cand.items():
        m = Chem.MolFromSmiles(smi)
        if m is None:
            print(f"FAIL  {name}: {smi}")
            continue
        formula = CalcMolFormula(m)
        if expect is not None and formula != expect:
            print(f"FAIL  {name}: formula {formula}, expected {expect}")
            continue
        print(f"OK    {name}: {formula}")
        ok[name] = smi
    return ok


def preview(ok: dict, out_path: Path):
    """Render a preview grid of the verified structures."""
    if not ok:
        return
    cols = min(3, len(ok))
    rows = (len(ok) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(3.4 * cols, 3.1 * rows),
                             dpi=90, squeeze=False)
    for ax, (name, smi) in zip(axes.flat, ok.items()):
        d = rdMolDraw2D.MolDraw2DCairo(380, 340)
        rdMolDraw2D.PrepareAndDrawMolecule(d, Chem.MolFromSmiles(smi))
        d.FinishDrawing()
        ax.imshow(Image.open(io.BytesIO(d.GetDrawingText())))
        ax.set_title(name, fontsize=11)
        ax.axis("off")
    for ax in axes.flat[len(ok):]:
        ax.axis("off")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=130)
    print("preview saved:", out_path)


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="Verify candidate SMILES + formulas.")
    common.add_common_args(ap)
    ap.add_argument("extra", nargs="*", default=[],
                    help='files with "name SMILES [formula]" lines to verify')
    args = ap.parse_args(argv)
    log = common.get_logger("verify_smiles", args.log_level)
    cfg = CONFIG
    if args.config:
        common.apply_overrides(cfg, args.config)
    common.resolve(cfg)

    cand = dict(CAND)
    for f in args.extra:
        p = Path(f)
        if not p.exists():
            common.fail(log, f"file not found: {p}", common.EXIT_INPUT)
        for ln in p.read_text(encoding="utf-8").splitlines():
            parts = ln.split()
            if len(parts) >= 2:
                cand[parts[0]] = (parts[1], parts[2] if len(parts) > 2 else None)
    if args.dry_run:
        log.info("would verify %d compounds; preview -> %s", len(cand),
                 cfg.out_dir / "figures" / "_smiles_preview.png")
        return common.EXIT_OK
    ok = verify(cand)
    preview(ok, CONFIG.out_dir / "figures" / "_smiles_preview.png")
    return common.EXIT_OK if len(ok) == len(cand) else common.EXIT_INPUT


if __name__ == "__main__":
    sys.exit(main())
