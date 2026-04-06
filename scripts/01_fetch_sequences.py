"""Fetch FASTA sequences from RCSB for all apo/holo PDB entries."""

import os
import sys
import requests

sys.path.insert(0, os.path.dirname(__file__))
from utils import load_config, file_exists_skip

FASTA_DIR = "data/fasta"
RCSB_FASTA_URL = "https://www.rcsb.org/fasta/entry/{pdb_id}/display"

os.makedirs(FASTA_DIR, exist_ok=True)


def fetch_fasta(pdb_id: str, chain: str) -> str | None:
    """Fetch FASTA from RCSB and extract the specified chain's sequence."""
    out_path = os.path.join(FASTA_DIR, f"{pdb_id}.fasta")
    if file_exists_skip(out_path, f"{pdb_id}.fasta"):
        return out_path

    url = RCSB_FASTA_URL.format(pdb_id=pdb_id)
    print(f"[fetch] {pdb_id} from {url}")
    resp = requests.get(url, timeout=30)
    if resp.status_code != 200:
        print(f"[error] {pdb_id}: HTTP {resp.status_code}")
        return None

    with open(out_path, "w") as f:
        f.write(resp.text)
    print(f"[ok] saved {out_path}")
    return out_path


def extract_chain_sequence(fasta_path: str, pdb_id: str, chain: str,
                           min_protein_length: int = 50) -> str | None:
    """Parse a multi-chain FASTA and return the sequence for the given chain.

    Applies min_protein_length filter to avoid returning short DNA/RNA
    sequences when a protein chain is expected.
    """
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

    pdb_upper = pdb_id.upper()

    # Priority 1: explicit chain label match with length check
    for header, seq in sequences.items():
        header_upper = header.upper()
        if pdb_upper in header_upper:
            if (f"CHAIN {chain.upper()}" in header_upper
                    or f"CHAINS {chain.upper()}" in header_upper
                    or f"{pdb_upper}_{chain.upper()}" in header_upper):
                if len(seq) >= min_protein_length:
                    return seq

    # Priority 2: single sequence
    if len(sequences) == 1:
        seq = list(sequences.values())[0]
        if len(seq) >= min_protein_length:
            print(f"[warn] {pdb_id}: chain {chain} not explicitly found, using only available sequence")
            return seq

    # Priority 3: longest sequence matching PDB ID (protein chains > DNA/RNA)
    candidates = [
        (header, seq) for header, seq in sequences.items()
        if pdb_upper in header.upper() and len(seq) >= min_protein_length
    ]
    if candidates:
        best_header, best_seq = max(candidates, key=lambda x: len(x[1]))
        print(f"[warn] {pdb_id}: chain {chain} ambiguous, using longest match: "
              f"{best_header[:60]} ({len(best_seq)} residues)")
        return best_seq

    print(f"[error] {pdb_id}: could not find chain {chain} with length >= {min_protein_length}")
    return None


def main():
    config = load_config()
    proteins = config["proteins"]

    for name, cfg in proteins.items():
        if cfg.get("skip_by_default"):
            print(f"[skip] {name} (skip_by_default)")
            continue

        chain = cfg.get("chain", "A")

        for pdb_key in ("apo_pdb", "holo_pdb"):
            pdb_id = cfg[pdb_key]
            fasta_path = fetch_fasta(pdb_id, chain)
            if fasta_path is None:
                continue
            seq = extract_chain_sequence(fasta_path, pdb_id, chain)
            if seq:
                print(f"  {name}/{pdb_key} chain {chain}: {len(seq)} residues")
            else:
                print(f"  [error] {name}/{pdb_key}: no sequence extracted")

    print("\nSequence fetching complete.")


if __name__ == "__main__":
    main()
