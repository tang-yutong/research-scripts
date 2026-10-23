"""Parse Multiwfn Mayer bond-order matrix files and extract atom-pair values.

The matrix file can contain several sections (open-shell systems have
alpha/beta/total; closed-shell has a single matrix); section names are taken
from the "Bond order matrix for ... electrons" headers.

Usage:
  python parse_bond_orders.py <bndmat.txt> [pairs like 10:12 ...]

Exit codes: 0 ok | 2 bad arguments | 3 unreadable input.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


def parse(path: Path) -> dict:
    """-> {section name: {(atom_i, atom_j): bond order}}."""
    sections: dict = {}
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    i = 0
    while i < len(lines):
        m = re.search(r"Bond order matrix for (.+?) electrons", lines[i])
        if m:
            name = m.group(1).strip()
            i += 1
            mat: dict = {}
            # blocks: a header line of column indices, then rows
            # "idx v1 v2 v3 v4 v5"
            while i < len(lines) and "Bond order matrix" not in lines[i]:
                header = re.findall(r"^\s*((?:\d+\s+)+)\s*$", lines[i])
                if header:
                    cols = [int(x) for x in header[0].split()]
                    i += 1
                    while i < len(lines):
                        row = re.match(
                            r"^\s*(\d+)\s+((?:-?\d+\.\d+\s+){1,6})\s*$",
                            lines[i])
                        if not row:
                            break
                        r_ = int(row.group(1))
                        vals = [float(x) for x in row.group(2).split()]
                        for c, v in zip(cols, vals):
                            mat[(r_, c)] = v
                        i += 1
                    continue
                i += 1
            sections[name] = mat
            continue
        i += 1
    return sections


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Extract atom-pair Mayer bond orders.")
    ap.add_argument("matrix", type=Path, help="bndmat.txt from Multiwfn")
    ap.add_argument("pairs", nargs="*", default=None,
                    help="atom pairs as I:J (default demo pair 10:12)")
    args = ap.parse_args(argv)

    if not args.matrix.exists():
        print(f"error: file not found: {args.matrix}", file=sys.stderr)
        return 3
    pairs = []
    try:
        for a in (args.pairs or ["10:12"]):
            x, y = a.split(":")
            pairs.append((int(x), int(y)))
    except ValueError:
        print(f"error: bad pair spec (expected I:J integers): {args.pairs}",
              file=sys.stderr)
        return 2

    secs = parse(args.matrix)
    print("file:", args.matrix.name)
    print("sections found:", list(secs.keys()))
    if not secs:
        print("error: no bond-order matrix sections parsed", file=sys.stderr)
        return 3
    for (x, y) in pairs:
        print(f"--- atoms {x} - {y} ---")
        for name, mat in secs.items():
            v = mat.get((x, y), mat.get((y, x)))
            print(f"  {name:12s}: {v if v is not None else 'N/A'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
