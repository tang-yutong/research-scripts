# md-analysis

Batch trajectory analysis for protein-ligand molecular dynamics in GROMACS,
plus MM-PBSA post-processing and publication figures. Built for a paired
study design: several systems, each simulated in two variants (`var_1` /
`var_2`), grouped into pairs (`pair_1`, ...) defined in `config.py`.

| Script | What it does |
|---|---|
| `gmx_batch_analyses.py` | Runs the GROMACS analysis suite over every system: protein C-alpha RMSD, ligand RMSD, RMSF, radius of gyration, SASA, protein-ligand hydrogen bonds. Probes index-group numbers via `make_ndx` and feeds selections by number; writes a human-readable command log. |
| `mmpbsa_parse.py` | Parses per-frame MM-PBSA output (energy terms) and per-residue decomposition into tidy CSV/JSON; cross-checks frame-averaged dH against the summary row; derives system composition (atoms, water, ions, chains, box) from the `.gro` file instead of hard-coding it. |
| `trajectory_plots.py` | All figures from the parsed data: paired variant panels, chain-aware RMSF, free-energy landscape (RMSD x Rg histogram), MM-PBSA term bars, per-residue bars; plus the equilibration-window statistics table. |
| `summary_table.py` | One CSV row per system: structural stats + recomputed dG_bind + per-pair dG difference + ADMET highlights. |

## Input layout (see `config.py`)

```
data/md/<system>_out/{fit.xtc, md_0_1.tpr}        # PBC-removed, fitted trajectory
data/md/<system>_int/complex_newbox.gro           # for system composition
data/mmpbsa/<pair>/<variant>/_pid~MMPBSA.dat      # + _pid~res_MMPBSA.dat
data/admet.json                                   # optional ADMET highlights
```

`gmx` must be on PATH. On Windows the GROMACS console may emit GBK; the
runner talks to subprocesses in bytes and decodes defensively.

## Notes

- MM-PBSA sign convention: the summary row's `-TdS=` value already IS -T*dS,
  so dG = dH + (-T*dS) is a **sum**. Getting this backwards flips every dG
  sign (learned the hard way).
- RMSF is plotted on a continuous residue index (some structures renumber
  mid-chain or concatenate chains); original numbering is used for tick
  labels and chain boundaries are detected from the data itself.
- Standard deviations are sample std (ddof=1) throughout.
