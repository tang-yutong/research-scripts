"""Tests for md-analysis: Vina log parsing, MM-PBSA parsing, xvg stats."""
from mmpbsa_parse import parse_mmpbsa, parse_res, mean_std, parse_gro
from rank_from_logs import parse_log
from trajectory_plots import read_xvg, window_stats

# ---- synthetic fixtures ----------------------------------------------------

VINA_LOG = """\
AutoDock Vina v1.2.5
Scoring function : vina
mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1       -11.87          0          0
   2        -8.31      1.918      3.204
"""
VINA_LOG_BAD = "empty run, no scoring table here"

MMPBSA_DAT = """\
some header
_pid~frame_001 dH= -30.000(-29.000) | MM= -35.000(-34.000) PB= 6.000 SA= -3.500 | COU= -20.000(-19.500) VDW= -15.000
_pid~frame_002 dH= -28.000(-27.200) | MM= -33.000(-32.000) PB= 5.000 SA= -3.000 | COU= -19.000(-18.500) VDW= -14.000
dH= -29.000(-28.100) | MM= -34.000(-33.000) PB= 5.500 SA= -3.250 | COU= -19.500(-19.000) VDW= -14.500
-TdS= 4.000(4.100)
dG= -25.000(-24.000) -5.975(-5.738)
"""

RES_DAT = """\
#Frame  P~12ALA(with DH)  P~45GLY(with DH)  L~1LIG(with DH)
mean   -2.500(0.100)      1.250(0.050)     -0.750(0.010)
"""

GRO = """\
system title
           14
    1PRO      N      1    0.000   0.000   0.000
    1PRO     CA      2    0.100   0.000   0.000
    2PRO      N      3    0.200   0.000   0.000
    2PRO     CA      4    0.300   0.000   0.000
    1PRO      N      5    0.400   0.000   0.000
    1PRO     CA      6    0.500   0.000   0.000
    1 MOL     C1     7    0.600   0.000   0.000
    1SOL     OW      8    0.700   0.000   0.000
    1SOL    HW1      9    0.800   0.000   0.000
    1SOL    HW2     10    0.900   0.000   0.000
    2SOL     OW     11    1.000   0.000   0.000
    2SOL    HW1     12    1.100   0.000   0.000
    2SOL    HW2     13    1.200   0.000   0.000
    1NA      NA     14    1.300   0.000   0.000
   4.00000   4.00000   4.00000
"""


class _Cfg:
    ligand_resname = "MOL"
    water_resnames = ("SOL", "WAT", "HOH")
    ion_resnames = ("NA", "CL")


# ---- tests -----------------------------------------------------------------

def test_vina_log_parse():
    assert parse_log(VINA_LOG) == -11.87
    assert parse_log(VINA_LOG_BAD) is None


def test_mmpbsa_frames_and_summary():
    import tempfile, pathlib
    with tempfile.NamedTemporaryFile("w", suffix=".dat", delete=False) as f:
        f.write(MMPBSA_DAT)
        p = pathlib.Path(f.name)
    try:
        frames, summary = parse_mmpbsa(p)
    finally:
        p.unlink()
    assert len(frames) == 2
    assert frames[0]["dH"] == -30.0 and frames[1]["VDW"] == -14.0
    assert summary["dH"][0] == -29.0
    assert summary["TdS"][0] == 4.0
    # sign convention: dG = dH + (-TdS) is a SUM -- the fixture is consistent
    assert abs(summary["dH"][0] + summary["TdS"][0] - summary["dG_kJ"][0]) < 1e-6
    assert abs(summary["dG_kcal"][0] + 5.975) < 1e-6


def test_mmpbsa_per_residue_and_stats():
    import tempfile, pathlib
    with tempfile.NamedTemporaryFile("w", suffix=".dat", delete=False) as f:
        f.write(RES_DAT)
        p = pathlib.Path(f.name)
    try:
        res = parse_res(p)
    finally:
        p.unlink()
    assert res == [("P~12ALA", -2.5), ("P~45GLY", 1.25), ("L~1LIG", -0.75)]
    # sample std (ddof=1), 0.0 for single values
    assert abs(mean_std([1.0, 2.0, 3.0])[1] - 1.0) < 1e-12
    assert mean_std([5.0])[1] == 0.0


def test_gro_composition():
    import tempfile, pathlib
    with tempfile.NamedTemporaryFile("w", suffix=".gro", delete=False) as f:
        f.write(GRO)
        p = pathlib.Path(f.name)
    try:
        info = parse_gro(p, _Cfg())
    finally:
        p.unlink()
    assert info["atoms"] == 14
    assert info["water_atoms"] == 6          # two 3-site waters
    assert info["ion_atoms"] == 1
    assert info["ligand_atoms"] == 1
    assert info["protein_atoms"] == 6
    assert info["chains"] == 2               # residue numbering resets once
    assert info["box"] == [4.0, 4.0, 4.0]


def test_xvg_read_and_window_stats():
    import tempfile, pathlib
    text = "# legend\n@ title\n" + "".join(
        f"{i * 0.5:6.2f} {10 + i:6.3f}\n" for i in range(41))  # t = 0..20 ns
    with tempfile.NamedTemporaryFile("w", suffix=".xvg", delete=False) as f:
        f.write(text)
        p = pathlib.Path(f.name)
    try:
        d = read_xvg(p)
    finally:
        p.unlink()
    assert d.shape == (41, 2)
    # last 10 ns window = t >= 10 -> 21 points, values 30..50
    m, s = window_stats(d, 10.0)
    assert abs(m - 40.0) < 1e-9
    assert abs(s - 6.2048368229954285) < 1e-9          # sample std (ddof=1): sqrt(38.5)
