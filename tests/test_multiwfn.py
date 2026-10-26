"""Tests for multiwfn-batch: bond-order matrix sections, BCP matching,
fragment growth and range compression. Value rows carry a trailing space,
matching real Multiwfn output."""
from pathlib import Path

from parse_bond_orders import parse
from match_bcp import read_cpprop, read_fchk_coords, dist
from auto_fragments import component, ranges

BONDMAT_OPEN_SHELL = (
    "Mayer bond order indices\n"
    "\n"
    " Bond order matrix for alpha electrons\n"
    "\n"
    "    1    2    3 \n"
    " 1  0.000000  0.900000  0.100000 \n"
    " 2  0.000000  0.000000  0.800000 \n"
    " 3  0.000000  0.000000  0.000000 \n"
    "\n"
    " Bond order matrix for beta electrons\n"
    "\n"
    "    1    2    3 \n"
    " 1  0.000000  0.100000  0.050000 \n"
    " 2  0.000000  0.000000  0.700000 \n"
    " 3  0.000000  0.000000  0.000000 \n"
    "\n"
    " Bond order matrix for total electrons\n"
    "\n"
    "    1    2    3 \n"
    " 1  0.000000  1.000000  0.150000 \n"
    " 2  0.000000  0.000000  1.500000 \n"
    " 3  0.000000  0.000000  0.000000 \n"
)

BONDMAT_CLOSED_SHELL = (
    " Bond order matrix for total electrons\n"
    "\n"
    "    1    2 \n"
    " 1  0.000000  1.050000 \n"
    " 2  0.000000  0.000000 \n"
)

CPPROP = """\
------ CP 1, Type (3,+1) ------
Position (Bohr):          5.0000      5.0000      5.0000
Density of all electrons: 0.0100
Laplacian of electron density: 0.0500

------ CP 2, Type (3,-1) ------
Position (Bohr):          0.5000      0.5000      0.5000
Density of all electrons: 0.2500
Laplacian of electron density: -0.8000
Energy density E(r) or H(r): -0.1500

------ CP 3, Type (3,-1) ------
Position (Bohr):          9.0000      9.0000      9.0000
Density of all electrons: 0.0010
Laplacian of electron density: 0.0020
"""

FCHK = """\
Generic fchk stub
Current cartesian coordinates               N=           9
  0.00000000    0.00000000    0.00000000
  1.00000000    0.00000000    0.00000000
  9.00000000    9.00000000    9.00000000
"""


def _tmp(text, suffix):
    import tempfile
    f = tempfile.NamedTemporaryFile("w", suffix=suffix, delete=False)
    f.write(text)
    f.close()
    return Path(f.name)


def test_bond_order_sections():
    p = _tmp(BONDMAT_OPEN_SHELL, ".txt")
    try:
        secs = parse(p)
    finally:
        p.unlink()
    assert set(secs) == {"alpha", "beta", "total"}
    assert secs["total"][(1, 2)] == 1.0
    assert secs["alpha"][(2, 3)] == 0.8
    assert secs["beta"][(1, 2)] == 0.1
    assert secs["total"][(2, 3)] == 1.5           # full-width rows: both triangles

    p = _tmp(BONDMAT_CLOSED_SHELL, ".txt")
    try:
        secs = parse(p)
    finally:
        p.unlink()
    assert list(secs) == ["total"]
    assert secs["total"][(1, 2)] == 1.05


def test_bcp_matching():
    pcp = _tmp(CPPROP, ".txt")
    fchk = _tmp(FCHK, ".fchk")
    try:
        cps = read_cpprop(pcp)
        coords = read_fchk_coords(fchk)
    finally:
        pcp.unlink()
        fchk.unlink()
    assert len(cps) == 3
    bcps = [c for c in cps if c["type"] == "3,-1"]
    assert len(bcps) == 2
    assert coords[0] == (0.0, 0.0, 0.0)
    # atoms 1-2 midpoint = (0.5, 0, 0); nearest BCP is CP 2 at (0.5, 0.5, 0.5)
    mid = tuple((a + b) / 2 for a, b in zip(coords[0], coords[1]))
    best = min(bcps, key=lambda c: dist(c["pos"], mid))
    assert best["idx"] == 2
    assert best["rho"] == 0.25 and best["h"] == -0.15


def test_fragment_growth_and_ranges():
    # atoms 0-2-3-4 connected, 1 isolated: component growth + range compression
    adj = [set() for _ in range(5)]
    for i, j in ((0, 2), (2, 3), (3, 4)):
        adj[i].add(j)
        adj[j].add(i)
    comp = component(4, adj)
    assert comp == {0, 2, 3, 4}
    assert 1 not in comp
    assert ranges({0, 1, 2, 5, 6, 9}) == "1-3,6-7,10"
    assert ranges(set()) == ""                  # empty set is a legal no-op
