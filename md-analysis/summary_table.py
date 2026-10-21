"""One summary CSV per configured system: structural stability, binding
free energy, the per-pair dG difference, and ADMET highlights.

Inputs: out/stats_eq.json (trajectory_plots), out/mmpbsa_frames.json
(mmpbsa_parse), and data/admet.json -- a user-provided ADMET export keyed
"<system>_<variant>" with fields Molecular Weight / Lipinski Rule / hERG.
dG_bind is recomputed here from per-frame dH + (-TdS), matching the
mmpbsa_parse convention. Exit codes: 0 ok | 3 missing input.
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CONFIG  # noqa: E402
import common  # noqa: E402

KJ2KCAL = 1 / 4.184


def binding_dg(frames_entry) -> float:
    """dG (kJ/mol) = mean(per-frame dH) + (-TdS); see mmpbsa_parse."""
    dh = statistics.mean(f["dH"] for f in frames_entry["frames"])
    return dh + frames_entry["summary"]["TdS"][0]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Aggregate per-system summary CSV.")
    common.add_common_args(ap)
    args = ap.parse_args(argv)
    log = common.get_logger("summary_table", args.log_level)
    cfg = CONFIG
    if args.config:
        common.apply_overrides(cfg, args.config)
    common.resolve(cfg)

    inputs = [cfg.out_dir / "stats_eq.json", cfg.out_dir / "mmpbsa_frames.json",
              cfg.data_dir / "admet.json"]
    common.require_inputs(log, *inputs)
    if args.dry_run:
        log.info("would write %s", cfg.table_dir / "summary_all.csv")
        return common.EXIT_OK

    stats = json.loads(inputs[0].read_text(encoding="utf-8"))
    frames = json.loads(inputs[1].read_text(encoding="utf-8"))
    admet = json.loads(inputs[2].read_text(encoding="utf-8"))["params"]

    cfg.table_dir.mkdir(parents=True, exist_ok=True)
    dg = {}
    variants = {}
    for i, (v1, v2) in enumerate(cfg.md_pairs, 1):
        for variant, sysname in (("var_1", v1), ("var_2", v2)):
            dg[sysname] = binding_dg({"frames": frames["frames"][sysname],
                                      "summary": frames["summary"][sysname]}) * KJ2KCAL
            variants[sysname] = variant
    ddg = {f"pair_{i}": dg[v2] - dg[v1] for i, (v1, v2) in enumerate(cfg.md_pairs, 1)}
    pair_of = {s: f"pair_{i}" for i, (v1, v2) in enumerate(cfg.md_pairs, 1)
               for s in (v1, v2)}

    with open(cfg.table_dir / "summary_all.csv", "w", encoding="utf-8-sig",
              newline="") as f:
        w = csv.writer(f)
        w.writerow(["system", "variant", "ca_rmsd_A", "ligand_rmsd_A", "rg_A",
                    "sasa_nm2", "hbonds", "dG_bind_kcal", "ddG_pair_kcal",
                    "mw_da", "lipinski_violations", "herg_block_prob"])
        for i, (v1, v2) in enumerate(cfg.md_pairs, 1):
            for sysname, variant in ((v1, "var_1"), (v2, "var_2")):
                s = stats[sysname]
                a = admet.get(f"{sysname}_{variant}", {})
                dd = f"{ddg[pair_of[sysname]]:+.2f}" if variant == "var_2" else "-"
                w.writerow([sysname, variant,
                            f"{s['rmsd_ca'][0]:.2f}+-{s['rmsd_ca'][1]:.2f}",
                            f"{s['rmsd_ligand'][0]:.2f}+-{s['rmsd_ligand'][1]:.2f}",
                            f"{s['gyrate'][0]:.2f}+-{s['gyrate'][1]:.2f}",
                            f"{s['sasa'][0]:.1f}+-{s['sasa'][1]:.1f}",
                            f"{s['hbond'][0]:.1f}+-{s['hbond'][1]:.1f}",
                            f"{dg[sysname]:.2f}", dd,
                            a.get("Molecular Weight", ""),
                            a.get("Lipinski Rule", ""), a.get("hERG", "")])
    common.write_manifest(cfg.out_dir / "work" / "summary_table",
                          Path(__file__).name, sys.argv[1:], cfg, inputs=inputs)
    log.info("wrote %s", cfg.table_dir / "summary_all.csv")
    for k, v in ddg.items():
        log.info("ddG %s: %+.2f kcal/mol", k, v)
    return common.EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
