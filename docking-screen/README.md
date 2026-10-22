# docking-screen

Pipeline for a large-ligand-library docking screen: prepare a SMILES library
into 3D PDBQT, dock it against one receptor with AutoDock Vina using N
parallel workers, validate the setup by redocking the co-crystal reference
ligand, fingerprint the top poses geometrically, rank from the raw logs, and
triage by pose quality + ADMET predictions.

## Pipeline order

1. `prep_ligands.ps1` -- ligand list (`ligands.smi`: `SMILES<TAB>L0001`) to
   3D PDBQT via OpenBabel; resumable; failures to `failures.csv`.
2. `dock_parallel.ps1` -- N parallel Vina workers, each taking ligands where
   `index % N == worker`; per-ligand timeout kills and quarantines hung
   jobs; existing outputs skipped. Receptor, box centre/size,
   exhaustiveness and timeout are parameters (box values in the script are
   examples -- set them to your binding site).
3. `dock_rescue.ps1` -- serial second pass over still-missing ligands with
   all cores per ligand.
4. `redock_rmsd.py` -- validation inside PyMOL: atom-name-matched RMSD
   between the co-crystal ligand and the redocked pose (runs in PyMOL, not
   from the shell).
5. `pose_fingerprint.py` -- for the top fraction of ranked ligands: centroid
   offset from the reference ligand, pocket occupancy, contacted residues,
   contacts to anchor residues (auto-derived from the reference ligand's
   contacts, not hand-picked), H-bond-like counts.
6. `rank_from_logs.py` -- best-mode affinity from every Vina log, bad logs
   skipped and reported -> ranked CSV.
7. `merge_screen_results.py` -- ranking x pose fingerprint x ADMET export
   with explicit pass/fail gates.

## Notes

- ADMET platforms return exports **without ligand ids** -- results must be
  matched back by row order, and the id spelling differs between files
  (`L0017` vs `lig_0017`). `merge_screen_results.py` encodes both traps.
- No DILI cutoff by design: predictors are systematically harsh on large
  lipophilic molecules; DILI is reported and judged with context.
- Ligand element types the scoring engine does not parameterize make Vina
  hang mid-run; check your library's element coverage and exclude such
  compounds up front.
- PowerShell 5.1+; `obabel` and `vina` on PATH; for conda OpenBabel builds
  set `BABEL_DATADIR` or `--gen3d` silently yields all-zero coordinates.
- Exit codes: 0 ok | 3 missing input | 4 tool failures (PowerShell scripts
  and Python entry points share this contract).
