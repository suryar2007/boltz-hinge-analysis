"""Generate Boltz input YAMLs (apo and holo) for each protein."""

import os
import sys
import yaml
import requests

sys.path.insert(0, os.path.dirname(__file__))
from utils import load_config

FASTA_DIR = "data/fasta"
APO_YAML_DIR = "data/yaml/apo"
HOLO_YAML_DIR = "data/yaml/holo"

os.makedirs(APO_YAML_DIR, exist_ok=True)
os.makedirs(HOLO_YAML_DIR, exist_ok=True)


def read_canonical_sequence(name: str) -> str | None:
    """Read the canonical (standardized) sequence for a protein."""
    path = os.path.join(FASTA_DIR, f"{name}_canonical.fasta")
    if not os.path.exists(path):
        print(f"[error] {name}: canonical FASTA not found at {path}")
        return None

    seq_lines = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                seq_lines.append(line)
    return "".join(seq_lines)


def extract_ligand_chain_sequence(pdb_id: str, chain: str) -> str | None:
    """Extract a chain's sequence from the RCSB FASTA for use as a ligand chain."""
    fasta_path = os.path.join(FASTA_DIR, f"{pdb_id}.fasta")
    if not os.path.exists(fasta_path):
        print(f"[warn] {pdb_id}: FASTA not found, attempting fetch")
        url = f"https://www.rcsb.org/fasta/entry/{pdb_id}/display"
        resp = requests.get(url, timeout=30)
        if resp.status_code != 200:
            print(f"[error] {pdb_id}: HTTP {resp.status_code}")
            return None
        with open(fasta_path, "w") as f:
            f.write(resp.text)

    sequences = {}
    current_header = None
    current_seq = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header is not None:
                    sequences[current_header] = "".join(current_seq)
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            sequences[current_header] = "".join(current_seq)

    for header, seq in sequences.items():
        header_upper = header.upper()
        pdb_upper = pdb_id.upper()
        if pdb_upper in header_upper:
            if (f"CHAIN {chain.upper()}" in header_upper
                    or f"CHAINS {chain.upper()}" in header_upper
                    or f"{pdb_upper}_{chain.upper()}" in header_upper):
                return seq

    # Fallback: return any chain that is NOT chain A
    for header, seq in sequences.items():
        if pdb_id.upper() in header.upper():
            if f"CHAIN A" not in header.upper() and f"CHAINS A" not in header.upper():
                return seq

    return None


def make_apo_yaml(name: str, sequence: str) -> str:
    """Create apo YAML (protein only)."""
    doc = {
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": sequence, "msa": "empty"}}
        ],
    }
    path = os.path.join(APO_YAML_DIR, f"{name}_apo.yaml")
    with open(path, "w") as f:
        yaml.dump(doc, f, default_flow_style=False, sort_keys=False)
    return path


def make_holo_yaml_ccd(name: str, sequence: str, ccd_code: str) -> str:
    """Create holo YAML with a CCD ligand."""
    doc = {
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": sequence, "msa": "empty"}},
            {"ligand": {"id": "B", "ccd": ccd_code}},
        ],
    }
    path = os.path.join(HOLO_YAML_DIR, f"{name}_holo.yaml")
    with open(path, "w") as f:
        yaml.dump(doc, f, default_flow_style=False, sort_keys=False)
    return path


def make_holo_yaml_chain(name: str, protein_seq: str, ligand_seq: str) -> str:
    """Create holo YAML with a separate protein/nucleotide chain as ligand."""
    doc = {
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": protein_seq, "msa": "empty"}},
            {"protein": {"id": "B", "sequence": ligand_seq, "msa": "empty"}},
        ],
    }
    path = os.path.join(HOLO_YAML_DIR, f"{name}_holo.yaml")
    with open(path, "w") as f:
        yaml.dump(doc, f, default_flow_style=False, sort_keys=False)
    return path


def main():
    config = load_config()
    proteins = config["proteins"]

    for name, cfg in proteins.items():
        if cfg.get("skip_by_default"):
            print(f"[skip] {name} (skip_by_default)")
            continue

        sequence = read_canonical_sequence(name)
        if not sequence:
            continue

        # Apo YAML
        apo_path = make_apo_yaml(name, sequence)
        print(f"[ok] {name} apo: {apo_path}")

        # Holo YAML
        ligand_ccd = cfg.get("ligand_ccd")
        ligand_chain = cfg.get("ligand_chain")

        if ligand_ccd:
            holo_path = make_holo_yaml_ccd(name, sequence, ligand_ccd)
            print(f"[ok] {name} holo (CCD={ligand_ccd}): {holo_path}")
        elif ligand_chain:
            holo_pdb = cfg["holo_pdb"]
            ligand_seq = extract_ligand_chain_sequence(holo_pdb, ligand_chain)
            if ligand_seq:
                holo_path = make_holo_yaml_chain(name, sequence, ligand_seq)
                print(f"[ok] {name} holo (chain {ligand_chain}, {len(ligand_seq)} residues): {holo_path}")
            else:
                print(f"[error] {name}: could not extract ligand chain {ligand_chain} from {holo_pdb}")
        else:
            # No ligand info — just create apo-style for holo too
            holo_path = make_apo_yaml(name, sequence)
            print(f"[warn] {name}: no ligand info, holo YAML is same as apo")

    print("\nYAML generation complete.")


if __name__ == "__main__":
    main()
