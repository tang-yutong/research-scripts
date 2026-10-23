"""Popen-based driver for Multiwfn's interactive interface.

Multiwfn has no batch interface for the weak-interaction analyses (IRI /
IGMH), so this driver choreographs one interactive session per task:
incremental stdin writes, log-file marker polling (no fixed sleeps), and
rename-before-overwrite copies between the two fragment-pair calls that
would otherwise clobber each other's output.

Implemented sessions:
  IRI   grid (manual spacing) -> scatter -> normal cubes -> covalent
        screening -> screened cubes
  IGMH  three fragments -> high-quality grid -> cubes -> per-pair dG
        integration for pairs (1,2) and (1,3)

Marker strings depend on the Multiwfn version's prompt wording; verify them
when upgrading (see README). Exit codes: 0 ok | 2 bad args | 3 missing
input | 4 task failure. ``--dry-run`` prints the planned stdin sequences.
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


class MultiwfnDriver:
    """Drives one Multiwfn process with incremental stdin and log polling."""

    def __init__(self, exe: str, fchk: Path, workdir: Path):
        self.workdir = workdir
        workdir.mkdir(parents=True, exist_ok=True)
        self.log_path = workdir / "driver_log.txt"
        self.log_f = open(self.log_path, "w", encoding="utf-8", errors="replace")
        self.proc = subprocess.Popen(
            [exe, str(fchk)],
            stdin=subprocess.PIPE,
            stdout=self.log_f,
            stderr=subprocess.STDOUT,
            cwd=workdir,
            text=True,
        )

    def send(self, s: str) -> None:
        self.proc.stdin.write(s + "\n")
        self.proc.stdin.flush()

    def wait_marker(self, marker: str, count: int = 1, timeout: int = 3600,
                    tag: str = "") -> None:
        """Poll the log until `marker` appears >= count times; fail loudly if
        the process dies or the timeout elapses."""
        start = time.time()
        while time.time() - start < timeout:
            self.log_f.flush()
            try:
                n = self.log_path.read_text(encoding="utf-8",
                                            errors="replace").count(marker)
            except (IOError, OSError):
                n = 0
            if n >= count:
                return
            if self.proc.poll() is not None:
                raise RuntimeError(
                    f"Multiwfn exited (rc={self.proc.returncode}) while waiting "
                    f"for '{marker}' x{count} [{tag}]")
            time.sleep(3)
        raise TimeoutError(
            f"timeout ({timeout}s) waiting for '{marker}' x{count} [{tag}]")

    def finish(self) -> None:
        """Exit the menus and reap the process; safe if it already died."""
        try:
            if self.proc.poll() is None:
                self.send("0")     # exit post-processing menu
                self.send("0")     # exit weak-interaction menu
                self.send("q")     # quit
                self.proc.stdin.close()
                self.proc.wait(timeout=30)
            else:
                self.proc.stdin.close()
        except Exception:
            self.proc.kill()
        finally:
            self.log_f.close()


def _copy_if_ok(log, src: Path, dst: Path, what: str) -> None:
    if src.exists() and src.stat().st_size > 0:
        shutil.copy2(src, dst)
        log.info("  saved %s (%d B)", dst.name, dst.stat().st_size)
    else:
        raise RuntimeError(f"{what} missing or empty after export: {src}")


def run_iri(cfg, log, fchk: Path, dest_name: str, dry: bool = False) -> bool:
    """Full IRI session; returns True on success."""
    task = f"iri_{dest_name}"
    wd = cfg.work_dir / "mwfn" / task
    dest = cfg.out_dir / "multiwfn" / "IRI" / dest_name
    if dry:
        log.info("[IRI] %s: would drive %s -> %s (grid %.2f Bohr)",
                 dest_name, fchk.name, dest, cfg.iri_grid_bohr)
        return True
    log.info("[IRI] %s: starting", dest_name)
    dest.mkdir(parents=True, exist_ok=True)
    drv = MultiwfnDriver(cfg.multiwfn, fchk, wd)
    try:
        drv.send("20")     # weak interaction analysis
        drv.send("4")      # IRI
        drv.send("4")      # manual grid spacing, whole system
        drv.send(str(cfg.iri_grid_bohr))
        # the IRI post-processing menu has no dedicated header; the menu text
        # of option 9 doubles as the completion marker (count 1 = first menu)
        drv.wait_marker("Screen out covalent bond regions", 1,
                        timeout=cfg.mwfn_timeout_s, tag="grid")
        log.info("[IRI] %s: grid done", dest_name)

        drv.send("2")      # scatter output
        drv.wait_marker("Column 1/2/3", 1, timeout=300, tag="scatter")
        drv.send("3")      # normal cube export
        drv.wait_marker("exported to func2.cub", 1, timeout=300, tag="cubes")
        for fn in ("func1.cub", "func2.cub"):
            _copy_if_ok(log, wd / fn, dest / fn, fn)
        log.info("[IRI] %s: normal cubes saved", dest_name)

        drv.send("9")      # screen out covalent regions; the menu is
        # re-printed afterwards -- wait for its SECOND appearance instead of
        # sleeping a fixed interval (race-free)
        drv.wait_marker("Screen out covalent bond regions", 2,
                        timeout=cfg.mwfn_timeout_s, tag="screen")
        drv.send("3")      # screened cube export (func2.cub 2nd time)
        drv.wait_marker("exported to func2.cub", 2, timeout=300,
                        tag="cubes_screened")
        for fn, dst_fn in (("func1.cub", "func1s.cub"), ("func2.cub", "func2s.cub")):
            _copy_if_ok(log, wd / fn, dest / dst_fn, fn)
        if (wd / "output.txt").exists():
            shutil.copy2(wd / "output.txt", dest / "output.txt")
        drv.finish()
        log.info("[IRI] %s: DONE", dest_name)
        return True
    except Exception as e:
        log.error("[IRI] %s: FAILED - %s", dest_name, e)
        try:
            drv.finish()
        except Exception:
            pass
        return False


def run_igmh(cfg, log, fchk: Path, dest_name: str, frag1: str, frag2: str,
             dry: bool = False) -> bool:
    """Full IGMH session with fragment pairs (1,2) and (1,3)."""
    task = f"igmh_{dest_name}"
    wd = cfg.work_dir / "mwfn" / task
    dest = cfg.out_dir / "multiwfn" / "IGMH" / dest_name
    if dry:
        log.info("[IGMH] %s: would drive %s -> %s (frags %s / %s / c)",
                 dest_name, fchk.name, dest, frag1, frag2)
        return True
    log.info("[IGMH] %s: starting", dest_name)
    dest.mkdir(parents=True, exist_ok=True)
    drv = MultiwfnDriver(cfg.multiwfn, fchk, wd)
    try:
        drv.send("20")     # weak interaction analysis
        drv.send("11")     # IGMH
        drv.send("3")      # 3 fragments
        drv.send(frag1)
        drv.send(frag2)
        drv.send("c")      # fragment 3 = complement
        drv.send("3")      # high-quality grid
        drv.wait_marker("Post-processing menu", 1, timeout=cfg.mwfn_timeout_s,
                        tag="grid")
        log.info("[IGMH] %s: grid done", dest_name)

        drv.send("3")      # cube export
        drv.wait_marker("exported to sl2r.cub", 1, timeout=300, tag="cubes")
        for fn in ("dg.cub", "dg_inter.cub", "dg_intra.cub", "sl2r.cub"):
            _copy_if_ok(log, wd / fn, dest / fn, fn)
        log.info("[IGMH] %s: cubes saved", dest_name)

        # pair (1,2): integrate, answer y for the pdb, copy IMMEDIATELY --
        # the second call overwrites atmdg.txt / IBSIW.txt / atmdg.pdb
        for pair, want in ((("1", "2"), 1), (("1", "3"), 2)):
            drv.send("6")
            drv.send(",".join(pair))
            drv.send("2")  # integration grid quality
            drv.wait_marker("outputted to atmdg.txt", want, timeout=2400,
                            tag=f"pair{pair[0]}{pair[1]}")
            drv.send("y")
            drv.wait_marker("Done! atmdg.pdb", want, timeout=60,
                            tag=f"pdb{pair[0]}{pair[1]}")
            tag = f"atmdg_{pair[0]}{pair[1]}"
            for src_fn, dst_fn in (("atmdg.txt", f"{tag}.txt"),
                                   ("IBSIW.txt", f"IBSIW_{pair[0]}{pair[1]}.txt"),
                                   ("atmdg.pdb", f"{tag}.pdb")):
                _copy_if_ok(log, wd / src_fn, dest / dst_fn, src_fn)
            log.info("[IGMH] %s: pair %s,%s done", dest_name, pair[0], pair[1])
        drv.finish()
        log.info("[IGMH] %s: DONE", dest_name)
        return True
    except Exception as e:
        log.error("[IGMH] %s: FAILED - %s", dest_name, e)
        try:
            drv.finish()
        except Exception:
            pass
        return False


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Drive Multiwfn IRI/IGMH sessions.")
    common.add_common_args(ap)
    ap.add_argument("mode", choices=["batch", "run-iri", "run-igmh", "test"],
                    help="batch: all configured species; run-iri/run-igmh: "
                         "single task; test: configured smoke-test species")
    ap.add_argument("fchk", nargs="?", help="wavefunction file (run-iri/run-igmh)")
    ap.add_argument("dest", nargs="?", help="output folder name (run-iri/run-igmh)")
    ap.add_argument("frag1", nargs="?", help="fragment 1 range, e.g. 174-199 (run-igmh)")
    ap.add_argument("frag2", nargs="?", help="fragment 2 range (run-igmh)")
    args = ap.parse_args(argv)
    log = common.get_logger("mwfn_driver", args.log_level)
    cfg = CONFIG
    if args.config:
        common.apply_overrides(cfg, args.config)
    common.resolve(cfg)

    tasks = []  # (kind, fchk, dest, frag1, frag2)
    if args.mode == "batch":
        for name, spec in cfg.run1_species.items():
            tasks.append(("igmh", cfg.dft_dir / spec["fchk"], name,
                          spec["frags"][0], spec["frags"][1]))
        for name in cfg.run1_iri:
            tasks.append(("iri", cfg.dft_dir / f"{name}.fchk", name, None, None))
    elif args.mode == "test":
        tasks.append(("igmh", cfg.dft_dir / f"{cfg.run1_test}.fchk", cfg.run1_test,
                      cfg.run1_test_frags[0], cfg.run1_test_frags[1]))
    elif args.mode == "run-iri":
        if not (args.fchk and args.dest):
            ap.error("run-iri needs <fchk> <dest>")
        tasks.append(("iri", Path(args.fchk), args.dest, None, None))
    elif args.mode == "run-igmh":
        if not (args.fchk and args.dest and args.frag1 and args.frag2):
            ap.error("run-igmh needs <fchk> <dest> <frag1> <frag2>")
        tasks.append(("igmh", Path(args.fchk), args.dest, args.frag1, args.frag2))

    common.require_inputs(log, *[t[1] for t in tasks])
    results = {}
    for kind, fchk, dest, f1, f2 in tasks:
        if kind == "iri":
            results[f"IRI_{dest}"] = run_iri(cfg, log, fchk, dest, dry=args.dry_run)
        else:
            results[f"IGMH_{dest}"] = run_igmh(cfg, log, fchk, dest, f1, f2,
                                               dry=args.dry_run)
    log.info("===== SUMMARY =====")
    for k, v in results.items():
        log.info("  %s: %s", k, "OK" if v else "FAILED")
    common.write_manifest(cfg.out_dir / "work" / "mwfn_driver",
                          Path(__file__).name, sys.argv[1:], cfg,
                          inputs=[t[1] for t in tasks])
    return common.EXIT_OK if all(results.values()) else common.EXIT_TOOL


if __name__ == "__main__":
    sys.exit(main())
