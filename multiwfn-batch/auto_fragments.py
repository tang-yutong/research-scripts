"""Derive IGMH fragment ranges from wavefunction-file connectivity.

Instead of hand-counting atom ranges, build a bond graph from the .fchk
coordinates (distance cutoff = 1.25 x sum of covalent radii), grow each of
the two fragments from its seed atoms, and assign the complement as the
third fragment. Validates that the two grown fragments are disjoint and
that the configured check atom (``metal``) lands in the complement.

Seed definitions live in ``fragment_seeds`` (config.py); all values there
are an example layout. Exit codes: 0 ok | 3 missing input | 4 check failed.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CONFIG  # noqa: E402
import common  # noqa: E402

# generic element table -- extend for your chemistry
Z2SYM = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 14: "Si", 15: "P",
         16: "S", 17: "Cl", 35: "Br", 30: "Zn"}
RCOV = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Si": 1.11,
        "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "Zn": 1.22}
BOHR_TO_A = 0.529177  # physical constant


def read_fchk(path: Path):
    """[(symbol, x, y, z) in Angstrom] from atomic numbers + coordinates."""
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()

    def arr(key):
        for i, ln in enumerate(lines):
            if ln.startswith(key):
                n = int(re.search(r"N=\s*(\d+)", ln).group(1))
                vals, j = [], i + 1
                while len(vals) < n:
                    vals += [float(x) for x in lines[j].split()]
                    j += 1
                return vals
        raise ValueError(f"missing array '{key}' in {path}")

    zs = [int(round(z)) for z in arr("Atomic numbers")]
    xyz = arr("Current cartesian coordinates")
    return [(Z2SYM.get(z, "X"), xyz[3 * i] * BOHR_TO_A,
             xyz[3 * i + 1] * BOHR_TO_A, xyz[3 * i + 2] * BOHR_TO_A)
            for i, z in enumerate(zs)]


def bonds(atoms):
    """Adjacency sets by the covalent-radius distance criterion."""
    n = len(atoms)
    adj = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = sum((atoms[i][k + 1] - atoms[j][k + 1]) ** 2 for k in range(3)) ** 0.5
            if d < 1.25 * (RCOV[atoms[i][0]] + RCOV[atoms[j][0]]):
                adj[i].add(j)
                adj[j].add(i)
    return adj


def component(start: int, adj):
    """Connected component containing `start` (0-based), via DFS."""
    seen, stack = set(), [start]
    while stack:
        u = stack.pop()
        if u in seen:
            continue
        seen.add(u)
        stack.extend(adj[u] - seen)
    return seen


def ranges(idxs) -> str:
    """Compress a 0-based index set into Multiwfn's '1,5-8,12' format."""
    s = sorted(i + 1 for i in idxs)
    if not s:
        return ""
    out, run = [], [s[0], s[0]]
    for v in s[1:]:
        if v == run[1] + 1:
            run[1] = v
        else:
            out.append(run)
            run = [v, v]
    out.append(run)
    return ",".join(str(a) if a == b else f"{a}-{b}" for a, b in out)


def fragment_ranges(fchk: Path, seeds: dict):
    """-> (frag_1, frag_2, frag_3, checks) Multiwfn range strings."""
    atoms = read_fchk(fchk)
    adj = bonds(atoms)
    f1 = set()
    for s in seeds["frag_1"]:
        f1 |= component(s - 1, adj)
    f2 = set()
    for s in seeds["frag_2"]:
        f2 |= component(s - 1, adj)
    overlap = f1 & f2
    rest = set(range(len(atoms))) - f1 - f2
    checks = {
        "atoms": len(atoms),
        "overlap_atoms": len(overlap),
        "metal_in_complement": (seeds["metal"] - 1) in rest,
        "frag1_seeds_in_frag1": all(s - 1 in f1 for s in seeds["frag_1"]),
        "frag2_seeds_in_frag2": all(s - 1 in f2 for s in seeds["frag_2"]),
    }
    return ranges(f1), ranges(f2), ranges(rest), checks


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Derive IGMH fragment ranges.")
    common.add_common_args(ap)
    ap.add_argument("--species", nargs="*", default=None,
                    help="subset of fragment_seeds keys (default: all)")
    args = ap.parse_args(argv)
    log = common.get_logger("auto_fragments", args.log_level)
    cfg = CONFIG
    if args.config:
        common.apply_overrides(cfg, args.config)
    common.resolve(cfg)

    names = args.species or list(cfg.fragment_seeds)
    common.require_inputs(log, *[cfg.dft_dir / cfg.fragment_seeds[n]["fchk"]
                                 for n in names])
    bad = []
    for name in names:
        seeds = cfg.fragment_seeds[name]
        f1, f2, f3, checks = fragment_ranges(
            cfg.dft_dir / seeds["fchk"], seeds)
        print(f"{name} ({seeds['fchk']}):")
        print(f"  frag_1: {f1}")
        print(f"  frag_2: {f2}")
        print(f"  frag_3 (complement, or 'c' at the prompt): {f3}")
        print(f"  checks: atoms={checks['atoms']}, "
              f"fragments disjoint={'yes' if checks['overlap_atoms'] == 0 else 'NO'}, "
              f"metal in frag_3={'yes' if checks['metal_in_complement'] else 'NO'}")
        if checks["overlap_atoms"] or not checks["metal_in_complement"]:
            bad.append(name)
    if bad:
        log.error("fragment checks failed for: %s (distance cutoff or seeds "
                  "need adjusting)", bad)
        return 4
    return 0


if __name__ == "__main__":
    sys.exit(main())
