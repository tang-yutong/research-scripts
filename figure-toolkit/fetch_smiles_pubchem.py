# -*- coding: utf-8 -*-
"""Fetch canonical SMILES + molecular formula from the PubChem REST API.

Structures that enter figures are fetched from the authority by compound
name -- never hand-typed -- so the figure and the database agree. Read-only;
exit codes: 0 ok (some failures tolerated, reported) | 3 all failed.
"""
import json
import sys
import urllib.parse
import urllib.request

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))
import common  # noqa: E402

# example lookups -- replace with the compounds your figures need
NAMES = ["benzene", "ethanol", "aspirin"]


def fetch_one(name: str):
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
           + urllib.parse.quote(name) + "/property/CanonicalSMILES,"
           "MolecularFormula/JSON")
    d = json.loads(urllib.request.urlopen(url, timeout=20).read())
    p = d["PropertyTable"]["Properties"][0]
    return p["CanonicalSMILES"], p["MolecularFormula"]


def main() -> int:
    log = common.get_logger("fetch_smiles_pubchem")
    failed = 0
    for n in NAMES:
        try:
            smi, formula = fetch_one(n)
            print(f"{n}: {smi}  ({formula})")
        except Exception as e:  # network errors reported, never fatal
            failed += 1
            print(f"{n}: FAIL {e}")
    return common.EXIT_OK if failed < len(NAMES) else common.EXIT_INPUT


if __name__ == "__main__":
    sys.exit(main())
