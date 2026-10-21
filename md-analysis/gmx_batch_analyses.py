"""Batch GROMACS trajectory analyses over all configured systems.

For every system directory ``data/md/<system>_out/`` this runs the standard
structural suite -- protein C-alpha RMSD, ligand RMSD, RMSF, radius of
gyration, SASA, and protein-ligand hydrogen bonds -- writing one .xvg per
analysis plus a human-readable command log. Index groups are probed with
``gmx make_ndx`` and selected by number (unambiguous); the ligand group is
matched by the configured residue name.

Inputs are read-only; runs are resumable (existing outputs are kept). A
manifest with timestamps, tool versions and input sizes is written to
``out/work/gmx_batch_analyses/manifest.json``.

Exit codes: 0 ok | 2 bad args | 3 missing input | 4 tool/analysis failure.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CONFIG  # noqa: E402
import common  # noqa: E402


class _Placeholder(dict):
    """Dry-run stand-in for the probed group table."""

    def __missing__(self, key):
        return "?"


def _decode(b: bytes) -> str:
    """Windows consoles may emit GBK; try UTF-8 first, fall back defensively."""
    for enc in ("utf-8", "gbk"):
        try:
            return b.decode(enc)
        except UnicodeDecodeError:
            continue
    return b.decode("utf-8", errors="replace")


def probe_groups(gmx: str, tpr: Path, workdir: Path) -> dict:
    """List default index groups via make_ndx; returns {group name: number}."""
    workdir.mkdir(parents=True, exist_ok=True)
    p = subprocess.run([gmx, "make_ndx", "-f", str(tpr), "-o", str(workdir / "probe.ndx")],
                       input=b"q\n", capture_output=True, timeout=300)
    if p.returncode != 0:
        raise RuntimeError(f"make_ndx failed (rc={p.returncode}) for {tpr}")
    out = _decode(p.stdout) + _decode(p.stderr)
    groups: dict = {}
    for m in re.finditer(r"^\s*(\d+)\s+(\S+)\s*:", out, re.M):
        groups.setdefault(m.group(2), m.group(1))
    return groups


def build_jobs(tag: str, tpr: Path, fit: Path, g: dict, xvg_dir: Path,
               ligand_group: str):
    """The six standard analyses for one system.

    Returns (output_name, gmx args, index-group selection) triples; group
    numbers come from the probed table so selections are unambiguous.
    """
    def x(name):
        return str(xvg_dir / f"{tag}_{name}.xvg")

    return [
        ("rmsd_ca", ["rms", "-s", str(tpr), "-f", str(fit), "-o", x("rmsd_ca"), "-tu", "ns"],
         [g["C-alpha"], g["C-alpha"]]),
        ("rmsd_ligand", ["rms", "-s", str(tpr), "-f", str(fit), "-o", x("rmsd_ligand"), "-tu", "ns"],
         [g["Protein"], g[ligand_group]]),
        ("rmsf", ["rmsf", "-s", str(tpr), "-f", str(fit), "-o", x("rmsf"), "-res"],
         [g["Protein"]]),
        ("gyrate", ["gyrate", "-s", str(tpr), "-f", str(fit), "-o", x("gyrate")],
         [g["non-Water"]]),
        ("sasa", ["sasa", "-s", str(tpr), "-f", str(fit), "-o", x("sasa"), "-surface"],
         [g["Protein"]]),
        ("hbond", ["hbond", "-s", str(tpr), "-f", str(fit), "-num", x("hbond")],
         [g["Protein"], g[ligand_group]]),
    ]


def process_system(log, cfg, tag, tpr, fit, log_lines) -> bool:
    """Run all six analyses for one system; True iff every job produced output."""
    g = probe_groups(cfg.gmx, tpr, cfg.work_dir / tag)
    need = ["Protein", "C-alpha", cfg.ligand_resname, "non-Water"]
    missing = [n for n in need if n not in g]
    if missing:
        log.error("[%s] missing index groups %s (available: %s)",
                  tag, missing, sorted(g)[:20])
        return False
    log.info("[%s] groups: Protein=%s C-alpha=%s %s=%s non-Water=%s",
             tag, g["Protein"], g["C-alpha"], cfg.ligand_resname,
             g[cfg.ligand_resname], g["non-Water"])
    names = {v: k for k, v in g.items()}
    ok_all = True
    for out_name, jargs, sel in build_jobs(tag, tpr, fit, g, cfg.xvg_dir,
                                           cfg.ligand_resname):
        cmd = [cfg.gmx] + jargs
        stdin_text = "\n".join(sel) + "\n"
        p = subprocess.run(cmd, input=stdin_text.encode(), capture_output=True,
                           timeout=cfg.gmx_timeout_s)
        out_path = Path(jargs[jargs.index("-o") + 1])
        ok = p.returncode == 0 and out_path.exists() and out_path.stat().st_size > 0
        log.info("    %-14s %s", out_name, "OK" if ok else f"FAIL rc={p.returncode}")
        tail = (_decode(p.stdout) + _decode(p.stderr))[-300:].replace("\n", " ")
        readable = "gmx " + " ".join(jargs[0:1]) + f" ... -o {out_name}"
        sel_txt = " | ".join(f"{names.get(x, x)}({x})" for x in sel)
        log_lines.append(f"### {tag} ###")
        log_lines.append(f"{readable}    # selection: {sel_txt}")
        log_lines.append(f"    status: {'OK' if ok else 'FAIL'}")
        if not ok:
            log_lines.append(f"    tail: {tail}")
            ok_all = False
    return ok_all


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Batch GROMACS trajectory analyses.")
    common.add_common_args(ap)
    args = ap.parse_args(argv)
    log = common.get_logger("gmx_batch_analyses", args.log_level)
    cfg = CONFIG
    if args.config:
        common.apply_overrides(cfg, args.config)
    common.resolve(cfg)

    systems = [s for pair in cfg.md_pairs for s in pair]
    layout = {s: (cfg.md_dir / f"{s}_out" / "md_0_1.tpr",
                  cfg.md_dir / f"{s}_out" / "fit.xtc") for s in systems}
    common.require_inputs(log, *[p for s in systems for p in layout[s]])

    cfg.xvg_dir.mkdir(parents=True, exist_ok=True)
    cfg.work_dir.mkdir(parents=True, exist_ok=True)

    if args.dry_run:
        for s in systems:
            tpr, fit = layout[s]
            for out_name, jargs, sel in build_jobs(s, tpr, fit, _Placeholder(),
                                                   cfg.xvg_dir, cfg.ligand_resname):
                log.info("would run: gmx %s   # groups %s", " ".join(jargs), sel)
        return common.EXIT_OK

    log_lines = [
        "GROMACS batch analysis command log",
        f"generated: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        "group selections were probed with make_ndx and fed by number",
        "=" * 80, "",
    ]
    failures = []
    for s in systems:
        tpr, fit = layout[s]
        if not process_system(log, cfg, s, tpr, fit, log_lines):
            failures.append(s)
    (cfg.out_dir / "gmx_command_log.txt").write_text("\n".join(log_lines),
                                                     encoding="utf-8")
    common.write_manifest(cfg.work_dir / "gmx_batch_analyses",
                          Path(__file__).name, sys.argv[1:], cfg,
                          inputs=[p for s in systems for p in layout[s]])
    if failures:
        log.error("failed systems/jobs: %s", failures)
        return common.EXIT_TOOL
    log.info("all systems done; command log: %s",
             cfg.out_dir / "gmx_command_log.txt")
    return common.EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
