"""Download true experimental structures (PDB + CIF) from RCSB."""

import os
import sys
import requests

sys.path.insert(0, os.path.dirname(__file__))
from utils import load_config, file_exists_skip

PDB_DIR = "data/pdb_cif"
os.makedirs(PDB_DIR, exist_ok=True)

RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"


def download_file(url: str, out_path: str, label: str) -> bool:
    if file_exists_skip(out_path, label):
        return True

    print(f"[fetch] {label} from {url}")
    resp = requests.get(url, timeout=60)
    if resp.status_code != 200:
        print(f"[error] {label}: HTTP {resp.status_code}")
        return False

    with open(out_path, "w") as f:
        f.write(resp.text)
    print(f"[ok] saved {out_path}")
    return True


def main():
    config = load_config()
    proteins = config["proteins"]

    for name, cfg in proteins.items():
        if cfg.get("skip_by_default"):
            print(f"[skip] {name} (skip_by_default)")
            continue

        for pdb_key in ("apo_pdb", "holo_pdb"):
            pdb_id = cfg[pdb_key]

            pdb_path = os.path.join(PDB_DIR, f"{pdb_id}_true.pdb")
            download_file(
                RCSB_PDB_URL.format(pdb_id=pdb_id),
                pdb_path,
                f"{pdb_id}.pdb",
            )

            cif_path = os.path.join(PDB_DIR, f"{pdb_id}_true.cif")
            download_file(
                RCSB_CIF_URL.format(pdb_id=pdb_id),
                cif_path,
                f"{pdb_id}.cif",
            )

    print("\nStructure fetching complete.")


if __name__ == "__main__":
    main()
