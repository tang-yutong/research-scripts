"""Rank AutoDock Vina results straight from their log files.

Every Vina log carries a scoring table; the mode-1 row is that ligand's
best-mode binding affinity. This walks a folder of per-ligand logs
(``log/lig_####.log``), extracts that number, joins ligand metadata from a
mapping CSV, sorts ascending, and writes a ranked CSV.

Malformed / empty / truncated logs are skipped and reported, never fatal --
in a screen of thousands of ligands a few dozen bad logs are routine.
Positive affinities are garbage poses of failed dockings: they are kept
as-is on purpose and only trimmed later at the histogram/plotting stage.

Exit codes: 0 ok | 2 bad args | 3 missing input / nothing ranked.
"""
from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CONFIG  # noqa: E402
import common  # noqa: E402

# first data row after the scoring-table header ("mode | affinity | ...")
ROW = re.compile(r"^\s*1\s+(-?\d+\.?\d*)")


def parse_log(text: str):
    """Best-mode affinity from one Vina log; None if no scoring table."""
    in_table = False
    for ln in text.splitlines():
        if "dist from best mode" in ln:
            in_table = True
            continue
        if in_table:
            m = ROW.match(ln)
            if m:
                return float(m.group(1))
    return None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Build a ranked CSV from AutoDock Vina log files.")
    common.add_common_args(ap)
    ap.add_argument("--logs-dir", type=Path, default=None,
                    help="folder of lig_####.log files (default: <docking_dir>/log)")
    ap.add_argument("--out", type=Path, default=None,
                    help="output CSV (default: <docking_dir>/final_ranking.csv)")
    args = ap.parse_args(argv)
    log = common.get_logger("rank_from_logs", args.log_level)
    cfg = CONFIG
    if args.config:
        common.apply_overrides(cfg, args.config)
    common.resolve(cfg)
    logs_dir = args.logs_dir or (cfg.docking_dir / "log")
    out_csv = args.out or (cfg.docking_dir / "final_ranking.csv")
    ligand_map = cfg.docking_dir / cfg.ligand_map
    common.require_inputs(log, ligand_map)
    if not logs_dir.is_dir():
        common.fail(log, f"logs directory not found: {logs_dir}", common.EXIT_INPUT)

    meta = {}
    with open(ligand_map, encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f):
            meta[r["num"]] = r                          # keyed "L0017"

    rows, bad = [], []
    for p in sorted(logs_dir.glob("*.log")):
        lig = p.stem                                    # "lig_0017"
        aff = parse_log(p.read_text(errors="ignore"))
        if aff is None:
            bad.append(lig)
            continue
        # log filename key "lig_0017" <-> mapping key "L0017"
        m = meta.get("L" + lig.split("_", 1)[1], {})
        rows.append({"ligand": lig, "affinity": aff, "name": m.get("name", ""),
                     "indication": m.get("indication", ""),
                     "target": m.get("target", "")})

    if not rows:
        common.fail(log, f"no parseable Vina logs in {logs_dir}", common.EXIT_INPUT)
    rows.sort(key=lambda r: r["affinity"])
    for i, r in enumerate(rows, 1):
        r["rank"] = i

    cols = ["ligand", "affinity", "name", "indication", "target", "rank"]
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    log.info("ranked %d ligands -> %s | skipped %d bad logs",
             len(rows), out_csv, len(bad))
    if bad:
        log.warning("skipped: %s", ", ".join(bad[:20]) +
                    (" ..." if len(bad) > 20 else ""))
    common.write_manifest(cfg.docking_dir, Path(__file__).name,
                          sys.argv[1:], cfg, inputs=[ligand_map])
    return common.EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
