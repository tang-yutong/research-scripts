"""Extract per-atom-pair interaction energies (delta-G pairs) from IGMH
atmdg files produced by mwfn_driver.py.

For every configured IGMH species this reads atmdg_12.txt, pulls the pairs
listed in the registry JSON (whatever tags the registry defines -- they are
data, not code), and lists the top-10 strongest pairs. Registry layout:
{"<species>": {"pairs": {"<tag>": [atom_i, atom_j]}}}.

Exit codes: 0 ok | 2 bad args | 3 missing input.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CONFIG  # noqa: E402
import common  # noqa: E402

PAIR_RE = re.compile(
    r"^\s+(\d+)\s+(\d+)\s+:\s+([\d.]+)\s+\(\s*([\d.]+)\s*%\s*\)")


def parse_atmdg(path: Path):
    """[(atom_i, atom_j, dG, percent)] rows from one atmdg file."""
    out = []
    for ln in open(path, encoding="utf-8", errors="ignore"):
        m = PAIR_RE.match(ln)
        if m:
            out.append((int(m.group(1)), int(m.group(2)),
                        float(m.group(3)), float(m.group(4))))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Extract IGMH delta-G pairs.")
    common.add_common_args(ap)
    ap.add_argument("--species", nargs="*", default=None,
                    help="subset of run2_igmh (default: all configured)")
    args = ap.parse_args(argv)
    log = common.get_logger("extract_dg_pairs", args.log_level)
    cfg = CONFIG
    if args.config:
        common.apply_overrides(cfg, args.config)
    common.resolve(cfg)

    common.require_inputs(log, cfg.registry)
    reg = json.loads(cfg.registry.read_text(encoding="utf-8"))
    species = args.species or cfg.run2_igmh
    igmh_dir = cfg.out_dir / "multiwfn" / "IGMH"
    if args.dry_run:
        for sp in species:
            log.info("would read %s", igmh_dir / sp / "atmdg_12.txt")
        return common.EXIT_OK

    lines = ["# IGMH delta-G pairs (full precision, auto-extracted)",
             "# atmdg_12.txt = interactions between fragments 1 and 2", ""]
    found = 0
    for sp in species:
        f12 = igmh_dir / sp / "atmdg_12.txt"
        if not f12.exists():
            lines.append(f"## {sp}  [atmdg_12.txt missing]")
            lines.append("")
            log.warning("%s: %s not found", sp, f12)
            continue
        pairs = parse_atmdg(f12)
        lines.append(f"## {sp}  ({len(pairs)} pairs)")
        want = []
        for tag, pair in reg.get(sp, {}).get("pairs", {}).items():
            want.append((tag, pair[0], pair[1]))
        for tag, a, b in want:
            hit = next((p for p in pairs if set(p[:2]) == {a, b}), None)
            if hit:
                lines.append(f"  {tag} #{a}-#{b}: dG = {hit[2]} ({hit[3]}%)")
            else:
                lines.append(f"  {tag} #{a}-#{b}: NOT FOUND (about 0)")
        top = sorted(pairs, key=lambda p: -p[2])[:10]
        lines.append("  Top-10:")
        for a, b, v, pct in top:
            lines.append(f"    #{a}-#{b}: {v} ({pct}%)")
        lines.append("")
        found += 1
    out = cfg.out_dir / "dg_pairs.md"
    out.write_text("\n".join(lines), encoding="utf-8")
    log.info("saved %s (%d/%d species with data)", out, found, len(species))
    common.write_manifest(cfg.out_dir / "work" / "extract_dg_pairs",
                          Path(__file__).name, sys.argv[1:], cfg,
                          inputs=[cfg.registry])
    return common.EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
