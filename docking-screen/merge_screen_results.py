"""Merge docking ranking, pose fingerprints and ADMET predictions into one
triage table with transparent pass/fail flags.

Inputs (all CSV / TSV):
  ranking : per-ligand best affinity + rank   (rank_from_logs.py output)
  poses   : per-ligand pose fingerprint       (pose_fingerprint.py output)
  smiles  : the submitted ADMET input list, one line "SMILES<TAB>L####";
            row ORDER must match the ADMET export
  admet   : raw ADMET predictor export, one row per submitted SMILES, in
            the SAME order, WITHOUT ligand ids

Two traps this script encodes: ADMET platforms return results without your
ligand ids, so rows must be matched back by position; and the pose table
keys on "lig_0017" while the SMILES list keys on "L0017" -- the two
spellings are converted explicitly (see lig_key).

Triage thresholds come from config.py. No DILI cutoff by design: predictors
are systematically harsh on large lipophilic molecules, so DILI is reported
as a column and judged with context instead of a threshold.

Exit codes: 0 ok | 2 bad args | 3 missing/mismatched input.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CONFIG  # noqa: E402
import common  # noqa: E402

# map the flags below to your predictor's actual column names
ADMET_COLS = {
    "BBB": "BBB", "hERG": "hERG", "Ames": "Ames", "DILI": "DILI",
    "Carcinogenicity": "Carcinogenicity", "CYP3A4inh": "CYP3A4-inh",
    "TPSA": "TPSA", "logP": "logP",
}
FIELDNAMES = ["ligand", "name", "rank", "affinity", "centroid_off", "anchor_hits",
              "hbonds", "BBB", "hERG", "Ames", "DILI", "Carcinogenicity",
              "CYP3A4inh", "TPSA", "logP", "pass_pose", "pass_admet", "final"]


def lig_key(l_num: str) -> str:
    """'L0017' -> 'lig_0017' (the two id spellings used across files)."""
    return "lig_" + l_num[1:]


def to_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def triage(row, pose, adm, cfg):
    """Attach pass_pose / pass_admet / final flags to one ranked row."""
    off = to_float(pose.get("centroid_off"))
    anchors = pose.get("anchor_hits", "0/8")
    bbb = to_float(adm.get(ADMET_COLS["BBB"]))
    herg = to_float(adm.get(ADMET_COLS["hERG"]))
    ames = to_float(adm.get(ADMET_COLS["Ames"]))
    pass_pose = (off is not None and off <= cfg.max_centroid_off_a
                 and anchors != "0/8")
    pass_admet = (bbb is not None and bbb >= cfg.bbb_min
                  and herg is not None and herg <= cfg.herg_max
                  and ames is not None and ames <= cfg.ames_max)
    row.update({
        "centroid_off": off if off is not None else "",
        "anchor_hits": anchors, "hbonds": pose.get("hbonds_like", ""),
        "BBB": bbb if bbb is not None else "", "hERG": herg if herg is not None else "",
        "Ames": ames if ames is not None else "",
        "DILI": adm.get(ADMET_COLS["DILI"], ""),
        "Carcinogenicity": adm.get(ADMET_COLS["Carcinogenicity"], ""),
        "CYP3A4inh": adm.get(ADMET_COLS["CYP3A4inh"], ""),
        "TPSA": adm.get(ADMET_COLS["TPSA"], ""),
        "logP": adm.get(ADMET_COLS["logP"], ""),
        "pass_pose": pass_pose, "pass_admet": pass_admet,
        "final": "KEEP" if (pass_pose and pass_admet) else "",
    })
    return row


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Merge ranking+poses+ADMET triage.")
    common.add_common_args(ap)
    ap.add_argument("--top", type=int, default=None,
                    help="triage only the top N ranked ligands "
                         "(default: top_fraction from config)")
    args = ap.parse_args(argv)
    log = common.get_logger("merge_screen_results", args.log_level)
    cfg = CONFIG
    if args.config:
        common.apply_overrides(cfg, args.config)
    common.resolve(cfg)

    ranking_csv = cfg.docking_dir / "final_ranking.csv"
    poses_csv = cfg.docking_dir / "interaction_report.csv"
    smiles = cfg.docking_dir / "admet_input.smi"
    admet_csv = cfg.docking_dir / "admet_export.csv"
    out_csv = cfg.docking_dir / "final_candidates.csv"
    common.require_inputs(log, ranking_csv, poses_csv, smiles, admet_csv)

    # 1) ADMET rows matched back to ligand ids by SMILES-list order
    with open(smiles, encoding="utf-8-sig") as f:
        order = [ln.split("\t")[1].strip() for ln in f if ln.strip()]
    with open(admet_csv, encoding="utf-8-sig", newline="") as f:
        admet_rows = list(csv.DictReader(f))
    if len(order) != len(admet_rows):
        common.fail(log, f"order mismatch: {len(order)} smiles vs "
                         f"{len(admet_rows)} admet rows", common.EXIT_INPUT)
    admet = {lig_key(l): row for l, row in zip(order, admet_rows)}

    # 2) pose fingerprints keyed by ligand id
    with open(poses_csv, encoding="utf-8-sig", newline="") as f:
        poses = {r["ligand"]: r for r in csv.DictReader(f)}

    # 3) walk the ranking, attach everything, apply gates
    with open(ranking_csv, encoding="utf-8-sig", newline="") as f:
        ranked = list(csv.DictReader(f))
    if not ranked:
        common.fail(log, f"empty ranking: {ranking_csv}", common.EXIT_INPUT)
    n_top = args.top if args.top is not None else \
        int(round(cfg.top_fraction * len(ranked)))
    ranked = ranked[:n_top]

    out = []
    for r in ranked:
        lig = r["ligand"]
        row = {"ligand": lig, "name": r.get("name", ""), "rank": r["rank"],
               "affinity": r["affinity"]}
        out.append(triage(row, poses.get(lig, {}), admet.get(lig, {}), cfg))

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES)
        w.writeheader()
        w.writerows(out)
    keeps = sum(1 for r in out if r["final"] == "KEEP")
    log.info("%d rows -> %s | pass_pose %d | pass_admet %d | KEEP %d",
             len(out), out_csv,
             sum(1 for r in out if r["pass_pose"]),
             sum(1 for r in out if r["pass_admet"]), keeps)
    common.write_manifest(cfg.docking_dir, Path(__file__).name,
                          sys.argv[1:], cfg,
                          inputs=[ranking_csv, poses_csv, smiles, admet_csv])
    return common.EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
