# -*- coding: utf-8 -*-
"""Redocking validation: atom-name-matched RMSD between the co-crystal
reference ligand and the best redocked pose.

Runs INSIDE PyMOL (File > Run script..., or `pymol -c -r this_file.py`),
with the co-crystal ligand loaded as object "ref" and the redocked pose
(MODEL 1 only, heavy atoms) loaded as "redock_0001".

Guards: ligands with zero matched atoms abort with a clear message instead
of dividing by zero; duplicated atom names keep their first occurrence and
the collision count is reported.

Rule of thumb: a redock RMSD below RMSD_GOOD_A means the docking setup
reproduces the crystal binding mode, so screen scores are interpretable.
"""
import math

from pymol import cmd  # noqa: F401  (PyMOL builtin)

RMSD_GOOD_A = 2.0  # acceptance threshold for a redocking setup


def collect(model):
    """[(name, [x, y, z]), ...] for the first occurrence of each atom name."""
    seen, out, dups = set(), [], 0
    for a in cmd.get_model(model).atom:
        if a.name in seen:
            dups += 1
            continue
        seen.add(a.name)
        out.append((a.name, list(a.coord)))
    return out, dups


def main():
    ref, ref_dups = collect("ref")
    lig, lig_dups = collect("redock_0001 and not elem H")
    if ref_dups or lig_dups:
        print(f"note: duplicate atom names ignored (ref: {ref_dups}, "
              f"ligand: {lig_dups})")
    ref_by_name = dict(ref)
    pairs = [(n, c, ref_by_name[n]) for n, c in lig if n in ref_by_name]
    unmatched = len(lig) - len(pairs)
    if not pairs:
        print("ERROR: no atoms matched by name -- check that both objects "
              "use the same ligand atom naming.")
        return
    total = sum(sum((c[i] - r[i]) ** 2 for i in range(3)) for _n, c, r in pairs)
    rmsd = math.sqrt(total / len(pairs))
    print(f"matched heavy atoms: {len(pairs)} (unmatched ligand atoms: {unmatched})")
    print("atom names:", [p[0] for p in pairs[:6]], "...")
    print(f"RMSD = {rmsd:.3f} A  ({'PASS' if rmsd < RMSD_GOOD_A else 'FAIL'}, "
          f"threshold {RMSD_GOOD_A} A)")


main()
