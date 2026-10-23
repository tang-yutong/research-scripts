"""Batch Mayer bond orders + minimal AIM topology for configured species.

For every species in ``run2_species`` (config.py) this drives two headless
Multiwfn runs in a per-species work directory:

  mayer input stream  9 / 1 / y           -> bndmat.txt   (bond-order matrix)
  topo  input stream  2 / 2 / 3 / 7 / -1  -> CPprop.txt   (critical points)

Resumable: existing outputs are skipped. Success is judged by the output
file existing and being non-empty -- NOT by the return code (some Multiwfn
builds exit with non-zero codes on success; see README). Failures are
collected and reported; the run exits 4 if anything failed.

Exit codes: 0 ok | 2 bad args | 3 missing input | 4 tool failure.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CONFIG  # noqa: E402
import common  # noqa: E402

MAYER_STREAM = "9\n1\ny\n"          # last 'y' confirms the bndmat.txt export
TOPO_STREAM = "2\n2\n3\n7\n-1\n"


def run_one(log, cfg, name: str, streams) -> bool:
    """Produce bond-order + critical-point outputs for one species."""
    wd = cfg.work_dir / "mwfn_batch" / name
    wd.mkdir(parents=True, exist_ok=True)
    fchk = cfg.dft_dir / f"{name}.fchk"
    if not fchk.exists():
        log.error("[%s] missing wavefunction file: %s", name, fchk)
        return False
    targets = [("bond_orders.txt", "bndmat.txt", MAYER_STREAM, "mayer"),
               ("critical_points.txt", "CPprop.txt", TOPO_STREAM, "topo")]
    ok = True
    for out_name, tmp, stream, tag in targets:
        dst = cfg.out_dir / "multiwfn" / f"{name}_{out_name}"
        if dst.exists():
            log.info("[%s] %s: outputs exist, skipping", name, tag)
            continue
        t0 = time.time()
        log_f = wd / f"{tag}.log"
        with open(log_f, "w", encoding="utf-8", errors="ignore") as lf:
            p = subprocess.run([cfg.multiwfn, str(fchk)], input=stream,
                               cwd=wd, stdout=lf, stderr=subprocess.STDOUT,
                               text=True, timeout=cfg.mwfn_timeout_s)
        src = wd / tmp
        if src.exists() and src.stat().st_size > 0:
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(src, dst)
            log.info("[%s] %s ok (%.0fs, rc=%s)", name, tag,
                     time.time() - t0, p.returncode)
        else:
            log.error("[%s] %s FAILED rc=%s tmp_exists=%s log=%s",
                      name, tag, p.returncode, src.exists(), log_f)
            ok = False
    return ok


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Batch Mayer bond orders + AIM topology.")
    common.add_common_args(ap)
    ap.add_argument("--species", nargs="*", default=None,
                    help="subset of run2_species (default: all configured)")
    args = ap.parse_args(argv)
    log = common.get_logger("batch_mayer_aim", args.log_level)
    cfg = CONFIG
    if args.config:
        common.apply_overrides(cfg, args.config)
    common.resolve(cfg)

    species = args.species or cfg.run2_species
    if args.dry_run:
        for name in species:
            log.info("would run mayer+topo for %s (%s) -> %s",
                     name, cfg.dft_dir / f"{name}.fchk",
                     cfg.out_dir / "multiwfn")
        return common.EXIT_OK

    failed = [n for n in species if not run_one(log, cfg, n, (MAYER_STREAM, TOPO_STREAM))]
    common.write_manifest(cfg.out_dir / "work" / "batch_mayer_aim",
                          Path(__file__).name, sys.argv[1:], cfg,
                          inputs=[cfg.dft_dir / f"{n}.fchk" for n in species])
    if failed:
        log.error("failed species: %s", failed)
        return common.EXIT_TOOL
    log.info("all species done")
    return common.EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
