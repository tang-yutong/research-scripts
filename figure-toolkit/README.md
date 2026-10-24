# figure-toolkit

Publication figures with chemical structures, built so that every structure
that enters a figure is traceable to an authoritative source and verified
before rendering.

| Script | What it does |
|---|---|
| `figkit.py` | Shared module: font/export conventions (PNG 600 dpi + SVG with text-as-paths), RDKit structure rendering with chemically aware motif highlights (spiro atoms = ring membership >= 2 with degree 4; bridgehead atoms = ring membership >= 2 with degree 3), and a proportional-width reaction-scheme layout engine. |
| `fetch_smiles_pubchem.py` | Canonical SMILES + molecular formula from the **PubChem REST API** by compound name -- structures come from the authority, never hand-typed. |
| `verify_smiles.py` | Gate before any figure: every SMILES parsed by RDKit, molecular formula recomputed and compared against expectation, preview grid for one visual pass. Accepts extra "name SMILES [formula]" files. |
| `example_md_workflow_figure.py` | Worked example: a standard MD simulation workflow diagram (preparation -> solvation -> minimization -> NVT/NPT -> production -> analysis -> binding free energy). |

## Why the fetch-then-verify pattern

A wrong structure in a figure is worse than no figure. The pattern --
canonical SMILES from PubChem, formula re-computed with RDKit and compared
against an expected value, preview grid eyeballed once -- catches
transcription-class errors before they reach a figure.
