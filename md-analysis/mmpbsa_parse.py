"""Parse MM-PBSA output files into tidy tables and JSON.

For every configured system this reads the per-frame energy file
(``_pid~MMPBSA.dat``: enthalpy terms per frame plus the summary block) and
the per-residue decomposition file, cross-checks the frame-averaged dH
against the summary row, and writes:

  out/tables/mmpbsa_terms.csv / _kJ.csv   per-term mean +- std, per system
  out/mmpbsa_perresidue.json              top residues by |dG|
  out/mmpbsa_frames.json                  per-frame values (plotting input)
  out/tables/system_info.csv              composition read from the .gro file

Sign convention (the trap this parser encodes): the summary row's ``-TdS=``
value already IS -T*dS, so dG = dH + (-TdS) is a SUM -- writing it as a
difference flips the sign of every binding free energy.

System composition (total atoms, water, ions, ligand atoms, chain count) is
parsed from ``data/md/<system>_int/complex_newbox.gro`` at run time rather
than hard-coded. Exit codes: 0 ok | 2 bad args | 3 missing input.
"""
from __future__ import annotations

import argparse
import json
import re
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CONFIG  # noqa: E402
import common  # noqa: E402

NUM = r"(-?\d+\.\d+)"
PAIR = re.compile(NUM + r"\(\s*" + NUM + r"\)")  # value(with-DH) pairs
KJ2KCAL = 1 / 4.184                              # physical conversion


def parse_mmpbsa(path: Path):
    """Parse one _pid~MMPBSA.dat -> (per-frame dicts, summary dict)."""
    frames, summary = [], {}
    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith("_pid~"):
                parts = line.split("|")
                # parts: [dH&MM..., PB SA, COU VDW, PB triple, SA triple]
                p0 = PAIR.search(parts[0])           # dH(plain)(DH)
                p1 = PAIR.search(parts[1])           # MM
                nums1 = re.findall(NUM, parts[1])    # MM pair + PB + SA
                p2 = PAIR.search(parts[2])           # COU
                vdw = float(re.findall(NUM, parts[2])[-1])
                frames.append({
                    "time": parts[0].split()[0],
                    "dH": float(p0.group(1)), "dH_DH": float(p0.group(2)),
                    "MM": float(p1.group(1)), "MM_DH": float(p1.group(2)),
                    "PB": float(nums1[2]), "SA": float(nums1[3]),
                    "COU": float(p2.group(1)), "COU_DH": float(p2.group(2)),
                    "VDW": vdw,
                })
            elif line.strip().startswith("dH="):
                parts = line.split("|")
                p0 = PAIR.search(parts[0])
                p1 = PAIR.search(parts[1])
                nums1 = re.findall(NUM, parts[1])
                p2 = PAIR.search(parts[2])
                summary["dH"] = (float(p0.group(1)), float(p0.group(2)))
                summary["MM"] = (float(p1.group(1)), float(p1.group(2)))
                summary["PB"] = float(nums1[2])
                summary["SA"] = float(nums1[3])
                summary["COU"] = (float(p2.group(1)), float(p2.group(2)))
                summary["VDW"] = float(re.findall(NUM, parts[2])[-1])
            elif line.strip().startswith("-TdS="):
                p = PAIR.search(line)
                summary["TdS"] = (float(p.group(1)), float(p.group(2)))
            elif line.strip().startswith("dG="):
                vals = re.findall(NUM, line)
                # kJ plain, kJ DH, kcal plain, kcal DH
                summary["dG_kJ"] = (float(vals[0]), float(vals[1]))
                summary["dG_kcal"] = (float(vals[2]), float(vals[3]))
    return frames, summary


def parse_res(path: Path):
    """Per-residue decomposition -> [(tag, mean_kJ), ...]; tags look like
    'P~12ALA(with DH)' / 'L~...' in the header row."""
    lines = path.read_text(errors="replace").splitlines()
    header = None
    mean_row = None
    for l in lines:
        if l.startswith("#Frame"):
            header = re.findall(r"([PL]~\d+\S*?)\(with DH\)", l)
        elif l.startswith("mean"):
            mean_row = l
    if header is None or mean_row is None:
        return []
    toks = re.findall(PAIR, mean_row)
    if len(toks) != len(header):
        print(f"  warning: {path.name} column mismatch {len(toks)} vs {len(header)}")
    return [(tag, float(v)) for tag, (v, _dh) in zip(header, toks)]


def mean_std(vals):
    """Mean and sample standard deviation (ddof=1); 0.0 for single values."""
    m = statistics.mean(vals)
    s = statistics.stdev(vals) if len(vals) > 1 else 0.0
    return m, s


def parse_gro(path: Path, cfg):
    """System composition straight from the coordinate file: total atoms,
    water/ion/ligand atom counts, protein chain count (residue numbering
    resets at each chain), and the box vectors."""
    water = set(cfg.water_resnames)
    ions = set(cfg.ion_resnames)
    counts = {"water_atoms": 0, "ion_atoms": 0, "ligand_atoms": 0, "protein_atoms": 0}
    ion_by_type: dict = {}
    chains, prev = 0, None
    with open(path, errors="replace") as f:
        f.readline()                                   # title
        natoms = int(f.readline().strip())
        for _ in range(natoms):
            ln = f.readline()
            resnr = int(ln[0:5])
            resname = ln[5:10].strip()
            if resname in water:
                counts["water_atoms"] += 1
            elif resname in ions:
                counts["ion_atoms"] += 1
                ion_by_type[resname] = ion_by_type.get(resname, 0) + 1
            elif resname == cfg.ligand_resname:
                counts["ligand_atoms"] += 1
            else:
                counts["protein_atoms"] += 1
                if prev is not None and resnr < prev:
                    chains += 1
                prev = resnr
        box = f.readline().split()
    return {"atoms": natoms, **counts,
            "chains": chains + (1 if counts["protein_atoms"] else 0),
            "ions_detail": " + ".join(f"{v} {k}" for k, v in ion_by_type.items()) or "0",
            "box": [float(x) for x in box]}


def iter_systems(cfg):
    """Yield (system, pair_name, variant, dat_path, res_path, gro_path) for
    every configured system: data/mmpbsa/<pair>/<variant>/_pid~MMPBSA.dat."""
    for i, (v1, v2) in enumerate(cfg.md_pairs, 1):
        pair = f"pair_{i}"
        for variant, sysname in (("var_1", v1), ("var_2", v2)):
            base = cfg.mmpbsa_dir / pair / variant
            yield (sysname, pair, variant,
                   base / "_pid~MMPBSA.dat",
                   base / "_pid~res_MMPBSA.dat",
                   cfg.md_dir / f"{sysname}_int" / "complex_newbox.gro")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Parse MM-PBSA outputs to tables/JSON.")
    common.add_common_args(ap)
    args = ap.parse_args(argv)
    log = common.get_logger("mmpbsa_parse", args.log_level)
    cfg = CONFIG
    if args.config:
        common.apply_overrides(cfg, args.config)
    common.resolve(cfg)

    jobs = list(iter_systems(cfg))
    common.require_inputs(log, *[j[3] for j in jobs], *[j[4] for j in jobs])
    if args.dry_run:
        for sysname, pair, variant, dat, _r, _g in jobs:
            log.info("would parse %s (%s/%s): %s", sysname, pair, variant, dat)
        return common.EXIT_OK

    cfg.table_dir.mkdir(parents=True, exist_ok=True)
    all_frames, all_summary, all_res = {}, {}, {}
    for sysname, pair, variant, dat, rdat, _gro in jobs:
        frames, summary = parse_mmpbsa(dat)
        res = parse_res(rdat)
        all_frames[sysname] = frames
        all_summary[sysname] = summary
        all_res[sysname] = res
        # cross-check: frame-averaged dH vs the summary row
        m = statistics.mean(f["dH"] for f in frames)
        log.info("%s: %d frames, dH mean %.3f vs summary %.3f, dG %.3f kJ (%.3f kcal)",
                 sysname, len(frames), m, summary["dH"][0],
                 summary["dG_kJ"][0], summary["dG_kcal"][0])

    # ---- energy terms: mean +- std (sample std, ddof=1) ----
    rows = []
    for sysname, pair, variant, *_ in jobs:
        fr = all_frames[sysname]
        su = all_summary[sysname]

        def ms(key):
            return mean_std([f[key] for f in fr])

        vdw_m, vdw_s = ms("VDW")
        cou_m, cou_s = ms("COU")
        pb_m, pb_s = ms("PB")
        sa_m, sa_s = ms("SA")
        dh_m, dh_s = ms("dH")
        tds = su["TdS"][0]      # single ensemble value, already equals -T*dS
        dg = dh_m + tds         # dG = dH + (-T*dS): a SUM (see module docstring)
        # -T*dS contributes no per-frame variance, so the dG spread is the dH spread
        rows.append({
            "system": sysname, "pair": pair, "variant": variant,
            "dE_vdw": (vdw_m, vdw_s), "dE_ele": (cou_m, cou_s),
            "dG_polar": (pb_m, pb_s), "dG_nonpolar": (sa_m, sa_s),
            "dH": (dh_m, dh_s), "minusTdS": (tds, 0.0),
            "dG_bind": (dg, dh_s),
            "dG_summary_kJ": su["dG_kJ"][0], "dG_summary_kcal": su["dG_kcal"][0],
        })

    n_frames = len(all_frames[jobs[0][0]])
    terms = ["dE_vdw", "dE_ele", "dG_polar", "dG_nonpolar", "dH", "minusTdS", "dG_bind"]
    for suffix, scale in (("kJ", 1.0), ("kcal", KJ2KCAL)):
        with open(cfg.table_dir / f"mmpbsa_terms_{suffix}.csv", "w",
                  encoding="utf-8-sig", newline="") as f:
            f.write(f"system,variant," + ",".join(terms) +
                    f" (mean+-std, {suffix}/mol, {n_frames} frames)\n")
            for r in rows:
                cells = [r["system"], r["variant"]]
                cells += [f"{m * scale:.2f}+-{s * scale:.2f}" for m, s in
                          (r[k] for k in terms)]
                f.write(",".join(cells) + "\n")

    # ---- per-pair dG difference (variant 2 minus variant 1), kcal/mol ----
    print("\nddG (var_2 - var_1, kcal/mol):")
    ddgs = {}
    for i, (v1, v2) in enumerate(cfg.md_pairs, 1):
        r1 = [r for r in rows if r["system"] == v1][0]
        r2 = [r for r in rows if r["system"] == v2][0]
        dd = (r2["dG_bind"][0] - r1["dG_bind"][0]) * KJ2KCAL
        ddgs[f"pair_{i}"] = dd
        print(f"  pair_{i} ({v2} - {v1}): {dd:+.2f}")

    # ---- per-residue top contributors ----
    res_out = {}
    for sysname, pair, variant, *_ in jobs:
        res = all_res[sysname]
        prot = [(t, v) for t, v in res if t.startswith("P~")]
        lig = [(t, v) for t, v in res if t.startswith("L~")]
        prot_sorted = sorted(prot, key=lambda x: abs(x[1]), reverse=True)
        res_out[sysname] = {
            "top_kJ": [(t, round(v, 2)) for t, v in prot_sorted[:15]],
            "ligand_row": [(t, round(v, 2)) for t, v in lig],
        }
        top = ", ".join(f"{t[2:]}:{v * KJ2KCAL:.2f}" for t, v in prot_sorted[:8])
        log.info("%s top residues (kcal): %s", sysname, top)
    (cfg.out_dir / "mmpbsa_perresidue.json").write_text(
        json.dumps(res_out, ensure_ascii=False, indent=1), encoding="utf-8")

    # ---- system info from the .gro files (cached; big files) ----
    cache = cfg.out_dir / "sysinfo_cache.json"
    info = json.loads(cache.read_text(encoding="utf-8")) if cache.exists() else {}
    for sysname, pair, variant, _d, _r, gro in jobs:
        if sysname in info:
            continue
        common.require_inputs(log, gro)
        info[sysname] = parse_gro(gro, cfg)
    cache.write_text(json.dumps(info, ensure_ascii=False, indent=1), encoding="utf-8")
    with open(cfg.table_dir / "system_info.csv", "w", encoding="utf-8-sig",
              newline="") as f:
        f.write("system,variant,total_atoms,water_molecules,counter_ions,"
                "protein_chains,ligand_atoms,box_type,box_length_nm\n")
        for sysname, pair, variant, *_ in jobs:
            c = info[sysname]
            b = c["box"]
            cubic = abs(b[0] - b[1]) < 1e-3 and abs(b[0] - b[2]) < 1e-3
            f.write(f"{sysname},{variant},{c['atoms']},"
                    f"{c['water_atoms'] // 3},{c['ions_detail']},{c['chains']},"
                    f"{c['ligand_atoms']},{'cubic' if cubic else 'rectangular'},"
                    f"{b[0]:.2f}\n")

    (cfg.out_dir / "mmpbsa_frames.json").write_text(
        json.dumps({"frames": all_frames,
                    "summary": all_summary, "ddG_kcal": ddgs},
                   ensure_ascii=False, indent=1), encoding="utf-8")
    common.write_manifest(cfg.out_dir / "work" / "mmpbsa_parse",
                          Path(__file__).name, sys.argv[1:], cfg,
                          inputs=[j[3] for j in jobs] + [j[4] for j in jobs])
    log.info("done -> %s", cfg.table_dir)
    return common.EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
