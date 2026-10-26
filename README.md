# research-scripts

Automation pipelines from my computational chemistry work: batch MD
trajectory analysis, a large-scale docking-screen workflow, batch
wavefunction analysis, publication figures, and a small job-vacancy radar.
Scripts are written with AI assistance and maintained by me.

## Modules

| Folder | What it does |
|---|---|
| `md-analysis/` | Batch GROMACS structural analyses (RMSD/RMSF/Rg/SASA/H-bonds) over configured systems, MM-PBSA output parsing with consistency checks, paired-variant figures, summary tables |
| `docking-screen/` | SMILES-to-PDBQT preparation, N-worker parallel AutoDock Vina with timeout isolation, redocking validation, geometric pose fingerprinting, log-based ranking, multi-criteria triage |
| `multiwfn-batch/` | Interactive-session driver for IRI/IGMH weak-interaction analysis, batch Mayer bond orders and AIM critical points, matrix/BCP parsers, fragment-range derivation |
| `figure-toolkit/` | RDKit structure rendering with motif highlights, PubChem-sourced and formula-verified SMILES, reaction-scheme and workflow-diagram layout |
| `job-radar/` | Daily PhD-vacancy digest over three official boards, with testable pure-function filters |

Every folder has its own README with usage, expected inputs and notes.

## Design principles

- **Idempotent and resumable** -- existing outputs are skipped, so any run
  can be interrupted and re-entered; nothing destructive happens by default.
- **Inputs are read-only** -- all outputs go to `out/`; expected data lives
  under `data/` and is never bundled.
- **External tools are isolated** -- gmx / vina / obabel / Multiwfn are
  invoked via subprocess with explicit timeouts and non-zero exits mapped to
  a shared exit-code contract (0 ok / 2 args / 3 input / 4 tool).
- **Loud failures** -- missing inputs and empty results stop with a readable
  message, not a traceback or a silent skip.
- **Run provenance** -- pipeline entries write `manifest.json` (timestamp,
  argv, config snapshot, input sizes, tool versions, git commit).

## Three things worth reading

1. **`multiwfn-batch/mwfn_driver.py`** -- driving an interactive
   Windows-only program that has no batch interface: incremental stdin
   choreography, log-marker polling instead of fixed sleeps, and
   rename-between-overwrites when two calls would clobber the same files.
2. **`docking-screen/dock_parallel.ps1`** -- N-worker sharding by index
   modulo with per-ligand timeout kill and quarantine, so one hung docking
   run cannot stall a multi-thousand-compound screen.
3. **`docking-screen/rank_from_logs.py` + `md-analysis/mmpbsa_parse.py`** --
   fault-tolerant parsing of dirty tool output at scale: skip-and-report for
   broken logs, paired-value regexes, and a sign convention (`-TdS` already
   IS -T*dS) encoded where it cannot be forgotten.

## Quick start

```
pip install pytest          # the only dev dependency
make test                   # or: python -m pytest -q
python md-analysis/gmx_batch_analyses.py --help
python docking-screen/rank_from_logs.py --help
python multiwfn-batch/mwfn_driver.py --dry-run batch
```

Configure paths, batch lists and thresholds in `config.py` (all values there
are an example layout), or override fields per run with `--config my.json`.

## Known limitations

- Windows-first: GBK-tolerant subprocess decoding, PowerShell 5.1+ scripts,
  CJK-safe file naming.
- Multiwfn menu-marker strings are version-pinned (see the module README).
- Scripts share `config.py`/`common.py` via a relative `sys.path` insert;
  there is no installer on purpose -- this is a working toolbox, not a
  package.

## License

MIT -- see [LICENSE](LICENSE).
