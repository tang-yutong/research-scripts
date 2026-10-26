"""Tests for docking-screen: id spelling mapping and triage gates."""
from merge_screen_results import lig_key, to_float, triage


class _Cfg:
    max_centroid_off_a = 4.0
    bbb_min = 0.5
    herg_max = 0.7
    ames_max = 0.5


def test_ligand_id_mapping():
    # the two id spellings used across files: SMILES lists use L####,
    # pose/log files use lig_####
    assert lig_key("L0017") == "lig_0017"
    assert lig_key("L0234") == "lig_0234"
    assert to_float("1.5") == 1.5
    assert to_float("not-a-number") is None
    assert to_float(None) is None


def test_triage_gates():
    cfg = _Cfg()
    pose_good = {"centroid_off": "2.5", "anchor_hits": "3/8", "hbonds_like": "4"}
    pose_far = {"centroid_off": "7.0", "anchor_hits": "3/8"}
    pose_noanchor = {"centroid_off": "1.0", "anchor_hits": "0/8"}
    admet_pass = {"BBB": "0.9", "hERG": "0.2", "Ames": "0.1"}

    row = triage({"ligand": "lig_0001"}, pose_good, admet_pass, cfg)
    assert row["pass_pose"] and row["pass_admet"] and row["final"] == "KEEP"

    row = triage({"ligand": "lig_0002"}, pose_far, admet_pass, cfg)
    assert not row["pass_pose"] and row["final"] == ""

    # zero anchor contacts disqualify the pose even at low centroid offset
    row = triage({"ligand": "lig_0003"}, pose_noanchor, admet_pass, cfg)
    assert not row["pass_pose"]

    # ADMET gate failures block the KEEP even with a good pose
    row = triage({"ligand": "lig_0004"}, pose_good,
                 {"BBB": "0.2", "hERG": "0.2", "Ames": "0.1"}, cfg)
    assert row["pass_pose"] and not row["pass_admet"] and row["final"] == ""

    # malformed ADMET values must not crash the triage (counted as failing)
    row = triage({"ligand": "lig_0005"}, pose_good,
                 {"BBB": "n/a", "hERG": "", "Ames": "0.1"}, cfg)
    assert not row["pass_admet"] and row["final"] == ""
