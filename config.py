"""Central configuration for every pipeline in this repository.

All paths resolve relative to this file (the repository root). Batch lists
and thresholds below describe an EXAMPLE project layout -- edit them to
match your own data, or override individual fields at runtime with
``--config my_values.json`` (see common.apply_overrides).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
OUT = ROOT / "out"


@dataclass
class Config:
    # -- paths ---------------------------------------------------------------
    data_dir: Path = DATA
    out_dir: Path = OUT
    work_dir: Path = OUT / "work"
    md_dir: Path = DATA / "md"             # <system>_out/ trajectories, <system>_int/ gro
    mmpbsa_dir: Path = DATA / "mmpbsa"     # <pair>/<variant>/_pid~MMPBSA.dat
    xvg_dir: Path = OUT / "xvg"
    figure_dir: Path = OUT / "figures"
    table_dir: Path = OUT / "tables"
    docking_dir: Path = DATA / "docking"   # receptor, ligand list, log/, poses/
    dft_dir: Path = DATA / "dft"           # fchk wavefunctions for Multiwfn
    registry: Path = DATA / "multiwfn" / "pair_registry.json"
    # -- external tools (expected on PATH; full paths also work) -------------
    gmx: str = "gmx"
    vina: str = "vina"
    obabel: str = "obabel"
    multiwfn: str = "Multiwfn"
    gmx_timeout_s: int = 3600             # per gmx invocation
    # -- md-analysis batches (example layout: 3 pairs x 2 variants) ----------
    md_pairs: list = field(default_factory=lambda: [
        ("sys_03", "sys_01"),   # pair_1 = (var_1 system, var_2 system)
        ("sys_06", "sys_02"),   # pair_2
        ("sys_05", "sys_04"),   # pair_3
    ])
    segment_labels: dict = field(default_factory=lambda: {
        # chains/subunits per system, in topological order (RMSF tick labels)
        "sys_03": ["chain_A", "extra_domain"],
        "sys_05": ["chain_A", "chain_B", "chain_C", "chain_D", "chain_E"],
    })
    ligand_resname: str = "MOL"
    water_resnames: tuple = ("SOL", "WAT", "HOH")
    ion_resnames: tuple = ("NA", "CL", "K", "MG")
    # -- docking-screen --------------------------------------------------------
    receptor: str = "receptor.pdbqt"
    ligand_list: str = "ligands.smi"       # lines: "SMILES<TAB>L0001"
    ligand_map: str = "ligand_map.csv"     # columns: num (L0001), name, ...
    ref_ligand: str = "ref_ligand.pdb"     # co-crystal pose of the reference ligand
    box_center: tuple = (30.0, 30.0, 30.0)  # example -- set to your binding site
    box_size: tuple = (24.0, 24.0, 24.0)    # example
    exhaustiveness: int = 8
    num_modes: int = 9
    workers: int = 8
    vina_timeout_s: int = 1200
    rescue_timeout_s: int = 900
    top_fraction: float = 0.10             # triage the top fraction of ranked ligands
    anchor_count: int = 8                  # pocket residues treated as anchors
    pocket_pad_a: float = 5.0              # pocket radius padding, A
    contact_d2: float = 16.0               # heavy-atom contact cutoff squared (4.0 A)
    hbond_d2: float = 12.25                # H-bond-like cutoff squared (3.5 A)
    pocket_d2: float = 25.0                # pocket-membership cutoff squared (5.0 A)
    max_centroid_off_a: float = 4.0        # pose gate: centroid offset, A
    bbb_min: float = 0.5                   # ADMET gates (predictor probabilities)
    herg_max: float = 0.7
    ames_max: float = 0.5
    # -- multiwfn-batch ---------------------------------------------------------
    run1_species: dict = field(default_factory=lambda: {
        # example: run label -> wavefunction file + IGMH fragment ranges
        "run1_sp_01": {"fchk": "run1_sp_01.fchk", "frags": ["60-85", "1-22"]},
        "run1_sp_03": {"fchk": "run1_sp_03.fchk", "frags": ["30-42,50-55", "56-59"]},
        "run1_sp_05": {"fchk": "run1_sp_05.fchk", "frags": ["60-85", "1-22"]},
    })
    run1_iri: list = field(default_factory=lambda: ["run1_sp_02", "run1_sp_06"])
    run1_test: str = "run1_sp_04"
    run1_test_frags: tuple = ("60-85", "1-22")
    fragment_seeds: dict = field(default_factory=lambda: {
        # example: run label -> wavefunction file + seed atom indices (1-based)
        # for auto_fragments.py: frag_1/frag_2 grow from their seeds,
        # frag_3 = the complement; "metal" marks a check atom inside frag_3
        "run1_sp_01": {"fchk": "run1_sp_01.fchk", "frag_1": [12, 15],
                       "frag_2": [4, 7], "metal": 9},
        "run1_sp_04": {"fchk": "run1_sp_04.fchk", "frag_1": [8, 11],
                       "frag_2": [3, 5], "metal": 6},
    })
    run2_species: list = field(default_factory=lambda: [
        f"run2_sp_{i:02d}" for i in range(1, 12)])
    run2_igmh: list = field(default_factory=lambda: [
        "run2_sp_02", "run2_sp_09", "run2_sp_04", "run2_sp_11", "run2_sp_05"])
    iri_grid_bohr: float = 0.20            # IRI grid spacing; finer costs RAM
    mwfn_timeout_s: int = 3600
    bcp_midpoint_ratio_warn: float = 0.35  # diagnostic warning ratio
    # -- trajectory plots ---------------------------------------------------------
    temperature_k: float = 310.0
    equilibration_window_ns: float = 50.0
    fel_bins: int = 40
    rmsf_smooth: int = 5
    rmsf_diff_tol_a: float = 0.15
    rmsf_min_run: int = 5
    # -- job-radar -------------------------------------------------------------
    radar_state_cap: int = 4000            # oldest seen-keys dropped beyond this


CONFIG = Config()
