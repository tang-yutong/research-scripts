# multiwfn-batch

Batch wavefunction analysis of Gaussian `.fchk` files with
[Multiwfn](https://multiwfn.mu-lab.net): Mayer bond orders, AIM critical
points, and IRI/IGMH weak-interaction analysis -- the analysis layer of
reaction-mechanism DFT studies.

| Script | What it does |
|---|---|
| `mwfn_driver.py` | The core: an interactive-session driver for a program with no batch interface (`Popen` + incremental stdin + log-marker polling + rename-before-overwrite). Full IRI sessions (grid, scatter, cubes, covalent screening) and IGMH sessions (three fragments, per-pair dG integration for pairs (1,2) and (1,3)). Batch orchestration over configured species lives in `main`. |
| `batch_mayer_aim.py` | Batch Mayer bond orders + minimal AIM topology per species; resumable; success judged by output files, not return codes. |
| `parse_bond_orders.py` | Bond-order matrix parser as a library + CLI. Handles open-shell (alpha/beta/total sections) and closed-shell (single matrix) files. |
| `match_bcp.py` | Geometric matching of (3,-1) critical points to atom pairs (nearest to midpoint, with a ratio warning) from `CPprop.txt` + fchk coordinates. |
| `extract_dg_pairs.py` | Per-atom-pair dG extraction from IGMH `atmdg` files; the pair registry is input data, so tags come from the registry, not the code. |
| `auto_fragments.py` | Derives IGMH fragment ranges from fchk connectivity (graph growth from seed atoms, complement as third fragment) with disjointness checks. |

## Hard-won notes (Multiwfn 3.8 dev, Windows)

- Normal exit codes are **not reliable** (some builds return 24, some 0 on
  success) -- judge success by the output file existing and being non-empty.
- The interactive-session marker strings polled by `mwfn_driver.py` depend
  on the version's prompt wording (e.g. `"exported to func2.cub"`,
  `"Screen out covalent bond regions"`, `"Post-processing menu"`); verify
  them when upgrading. The IRI post-processing menu has no dedicated header,
  so option 9's menu text doubles as the completion marker.
- ~200-atom systems: IRI grid mode 4 at 0.20 Bohr (finer grids run out of
  memory); IGMH uses the high-quality preset (mode 3).
- When running with a working directory, Multiwfn looks for `settings.ini`
  next to the executable; copy it into the work directory or accept the
  default 4 threads.
- The two IGMH fragment-pair calls overwrite each other's `atmdg.txt` /
  `IBSIW.txt` / `atmdg.pdb` -- the driver copies and renames between them.
