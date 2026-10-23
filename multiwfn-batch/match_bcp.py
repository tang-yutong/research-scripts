"""Match (3,-1) bond critical points in a Multiwfn CPprop.txt to atom pairs
by geometry (nearest to the pair midpoint), reading coordinates from the
matching .fchk file.

Usage:
  python match_bcp.py <CPprop.txt> <wavefunction.fchk> [pairs like 10:12 ...]

A BCP whose distance to the midpoint exceeds bcp_midpoint_ratio_warn x the
A-B distance (config.py, example 0.35) is flagged -- the pair may not have
a BCP at all. Exit codes: 0 ok | 2 bad args | 3 missing input.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CONFIG  # noqa: E402

BOHR_PER_A = 1.8897259886  # physical constant
HEADER = re.compile(r"-+\s+CP\s+(\d+),\s+Type \((3,[+-]\d)\)\s+-+\s*\n")


def read_fchk_coords(path: Path):
    """[(x, y, z), ...] in Bohr from 'Current cartesian coordinates'."""
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    for i, ln in enumerate(lines):
        if ln.startswith("Current cartesian coordinates"):
            n = int(re.search(r"N=\s*(\d+)", ln).group(1))
            vals, j = [], i + 1
            while len(vals) < n:
                vals += [float(x) for x in lines[j].split()]
                j += 1
            return [(vals[3 * k], vals[3 * k + 1], vals[3 * k + 2])
                    for k in range(n // 3)]
    raise ValueError("no 'Current cartesian coordinates' block found")


def read_cpprop(path: Path):
    """[CP dicts] with type, position (Bohr), rho, laplacian, energy density."""
    txt = path.read_text(encoding="utf-8", errors="ignore")
    marks = list(HEADER.finditer(txt))
    cps = []
    for k, m in enumerate(marks):
        seg_end = marks[k + 1].start() if k + 1 < len(marks) else len(txt)
        seg = txt[m.end():seg_end]
        pm = re.search(
            r"Position \(Bohr\):\s*([-\d.+E]+)\s+([-\d.+E]+)\s+([-\d.+E]+)", seg)
        rho = re.search(r"Density of all electrons:\s*([-\d.+E]+)", seg)
        lap = re.search(r"Laplacian of electron density:\s*([-\d.+E]+)", seg)
        h = re.search(r"Energy density E\(r\) or H\(r\):\s*([-\d.+E]+)", seg)
        if pm:
            cps.append({
                "idx": int(m.group(1)), "type": m.group(2),
                "pos": tuple(float(pm.group(i)) for i in (1, 2, 3)),
                "rho": float(rho.group(1)) if rho else None,
                "lap": float(lap.group(1)) if lap else None,
                "h": float(h.group(1)) if h else None,
            })
    return cps


def dist(p, q):
    return sum((a - b) ** 2 for a, b in zip(p, q)) ** 0.5


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Match BCPs to atom pairs by geometry.")
    ap.add_argument("cpprop", type=Path, help="CPprop.txt from Multiwfn")
    ap.add_argument("fchk", type=Path, help="wavefunction file with coordinates")
    ap.add_argument("pairs", nargs="*", default=None,
                    help="atom pairs as I:J (1-based)")
    args = ap.parse_args(argv)
    for p in (args.cpprop, args.fchk):
        if not p.exists():
            print(f"error: file not found: {p}", file=sys.stderr)
            return 3
    try:
        pairs = [tuple(int(x) for x in a.split(":"))
                 for a in (args.pairs or ["1:2"])]
    except ValueError:
        print(f"error: bad pair spec (expected I:J integers): {args.pairs}",
              file=sys.stderr)
        return 2

    coords = read_fchk_coords(args.fchk)
    cps = read_cpprop(args.cpprop)
    bcps = [c for c in cps if c["type"] == "3,-1"]
    print(f"CPprop parsed: {len(cps)} CPs total, {len(bcps)} BCPs")
    for a, b in pairs:
        if not (1 <= a <= len(coords) and 1 <= b <= len(coords)):
            print(f"--- atoms {a}-{b}: index out of range "
                  f"(1..{len(coords)}) ---")
            continue
        pa, pb = coords[a - 1], coords[b - 1]
        mid = tuple((x + y) / 2 for x, y in zip(pa, pb))
        dab = dist(pa, pb)
        print(f"--- atoms {a}-{b} (d = {dab:.3f} Bohr = {dab / BOHR_PER_A:.3f} A) ---")
        if not bcps:
            print("  no (3,-1) critical points in file")
            continue
        near = sorted(bcps, key=lambda c: dist(c["pos"], mid))[:2]
        best = near[0]
        dmid = dist(best["pos"], mid)
        print(f"  nearest BCP #{best['idx']}  dist_to_midpoint = {dmid:.3f} "
              f"Bohr ({dmid / BOHR_PER_A:.3f} A)")
        print(f"    rho = {best['rho']}")
        print(f"    lap = {best['lap']}")
        if best["h"] is not None:
            print(f"    H(r) = {best['h']}")
        if dmid / dab > CONFIG.bcp_midpoint_ratio_warn:
            print("    warning: BCP far from the midpoint -- it may not belong "
                  "to this pair (the bond may have no BCP)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
