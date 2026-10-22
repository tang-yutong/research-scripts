"""Pose-quality fingerprints for the top-ranked docking poses.

For each pose of the top fraction of the ranked ligands (see
``top_fraction`` in config.py) this computes, relative to the co-crystal
reference ligand (``ref_ligand.pdb``):

  centroid_off    distance between ligand and reference centroids, A
  in_pocket_pct   share of ligand heavy atoms inside the pocket sphere
  pocket_contacts receptor residues contacted (any heavy atom within 4.0 A)
  anchor_hits     contacts to the N pocket residues that contact the
                  reference ligand most (auto-derived, not hand-picked)
  hbonds_like     ligand-receptor polar pairs within 3.5 A

The pocket is defined geometrically from the reference ligand; no
residue is hard-coded. Output: interaction_report.csv with a fixed column
layout (header-only if nothing qualified, never a crash on empty input).

Exit codes: 0 ok | 2 bad args | 3 missing input | 4 (unused here).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CONFIG  # noqa: E402
import common  # noqa: E402

FIELDNAMES = ["ligand", "affinity", "rank", "name", "centroid_off",
              "in_pocket_pct", "pocket_contacts", "anchor_hits", "hbonds_like"]


def parse_atoms(path: Path, model1_only: bool = True):
    """[(chain, resi, resn, x, y, z, atomtype), ...] from PDBQT/PDB ATOM rows."""
    atoms, seen_model = [], False
    for line in open(path, errors="replace"):
        if line.startswith("MODEL"):
            if seen_model and model1_only:
                break
            seen_model = True
            continue
        if line.startswith(("ATOM", "HETATM")):
            atype = line[77:79].strip()
            atoms.append((line[21],                    # chain
                          line[22:26].strip(),         # resi
                          line[17:20].strip(),         # resn
                          float(line[30:38]), float(line[38:46]), float(line[46:54]),
                          atype))
    return atoms


def heavy(a):
    return not a[6].startswith("H")


def d2(p, q):
    return (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2


def fingerprint(lig_atoms, ref_center, r_pocket2, pocket_atoms, anchors, rec_pol, cfg):
    """(centroid_off, in_pocket_frac, n_contacts, n_anchor_hits, n_hbonds_like)."""
    lh = [a for a in lig_atoms if heavy(a)]
    lxyz = [(a[3], a[4], a[5]) for a in lh]
    lcen = tuple(sum(c[i] for c in lxyz) / len(lxyz) for i in range(3))
    off = math.sqrt(d2(lcen, ref_center))
    in_pok = sum(1 for c in lxyz if d2(c, ref_center) <= r_pocket2) / len(lxyz)

    cont_res, anchor_hit = set(), set()
    for a in pocket_atoms:
        ax = (a[3], a[4], a[5])
        for c in lxyz:
            if d2(ax, c) <= cfg.contact_d2:
                cont_res.add((a[0], a[1], a[2]))
                if (a[0], a[1], a[2]) in anchors:
                    anchor_hit.add((a[0], a[1], a[2]))
                break
    hb = sum(1 for a in lh if a[6][0] in "NOS"
             and any(d2((a[3], a[4], a[5]), r) <= cfg.hbond_d2 for r in rec_pol))
    return off, in_pok, len(cont_res), len(anchor_hit), hb


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Fingerprint top docking poses.")
    common.add_common_args(ap)
    args = ap.parse_args(argv)
    log = common.get_logger("pose_fingerprint", args.log_level)
    cfg = CONFIG
    if args.config:
        common.apply_overrides(cfg, args.config)
    common.resolve(cfg)

    ranking = cfg.docking_dir / "final_ranking.csv"
    receptor = cfg.docking_dir / cfg.receptor
    ref_lig = cfg.docking_dir / cfg.ref_ligand
    common.require_inputs(log, ranking, receptor, ref_lig)
    if args.dry_run:
        log.info("would fingerprint top %.0f%% poses -> %s",
                 cfg.top_fraction * 100, cfg.docking_dir / "interaction_report.csv")
        return common.EXIT_OK

    rec = parse_atoms(receptor)
    rec_heavy = [a for a in rec if heavy(a)]
    ref = [a for a in parse_atoms(ref_lig) if heavy(a)]
    if not ref:
        common.fail(log, f"reference ligand has no heavy atoms: {ref_lig}",
                    common.EXIT_INPUT)
    ref_xyz = [(a[3], a[4], a[5]) for a in ref]
    cen = tuple(sum(c[i] for c in ref_xyz) / len(ref_xyz) for i in range(3))
    r_pocket2 = max(d2(c, cen) for c in ref_xyz) + cfg.pocket_pad_a ** 2

    # pocket residues: any receptor heavy atom within pocket_d2 of the reference
    pocket_res: dict = {}
    for a in rec_heavy:
        if min(d2((a[3], a[4], a[5]), c) for c in ref_xyz) <= cfg.pocket_d2:
            pocket_res.setdefault((a[0], a[1], a[2]), []).append(a)
    pocket_atoms = [a for v in pocket_res.values() for a in v]
    # anchors: the residues contacting the reference ligand most
    anchor_rank = sorted(pocket_res.items(),
                         key=lambda kv: -sum(1 for a in kv[1]
                                             if min(d2((a[3], a[4], a[5]), c)
                                                    for c in ref_xyz) <= cfg.contact_d2))
    anchors = {k for k, _v in anchor_rank[:cfg.anchor_count]}
    rec_pol = [(a[3], a[4], a[5]) for a in rec if a[6][0] in "NO"]
    log.info("pocket residues: %d, anchors: %s", len(pocket_res),
             sorted(k[2] + k[1] for k in anchors))

    with open(ranking, encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    n_top = int(round(cfg.top_fraction * len(rows)))
    out = []
    for r in rows[:n_top]:
        lig = r["ligand"]
        pose = cfg.docking_dir / "poses" / f"{lig}.pdbqt"
        if not pose.exists():
            continue
        off, inp, nc, na, hb = fingerprint(parse_atoms(pose), cen, r_pocket2,
                                           pocket_atoms, anchors, rec_pol, cfg)
        out.append({"ligand": lig, "affinity": r["affinity"], "rank": r["rank"],
                    "name": r["name"], "centroid_off": round(off, 2),
                    "in_pocket_pct": round(inp * 100, 0), "pocket_contacts": nc,
                    "anchor_hits": f"{na}/{cfg.anchor_count}", "hbonds_like": hb})
    if not out:
        common.fail(log, "no poses fingerprinted (expected under "
                         f"{cfg.docking_dir / 'poses'}); run docking first",
                    common.EXIT_INPUT)

    with open(cfg.docking_dir / "interaction_report.csv", "w", newline="",
              encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES)
        w.writeheader()
        w.writerows(out)
    well = sum(1 for o in out if o["centroid_off"] <= cfg.max_centroid_off_a
               and o["anchor_hits"] != f"0/{cfg.anchor_count}")
    log.info("wrote %d rows (%d well-positioned) -> interaction_report.csv",
             len(out), well)
    common.write_manifest(cfg.docking_dir, Path(__file__).name,
                          sys.argv[1:], cfg, inputs=[ranking, receptor, ref_lig])
    return common.EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
