"""Check apo/holo sequence identity and standardize on apo sequence."""

import os
import sys
import csv

from Bio import pairwise2

sys.path.insert(0, os.path.dirname(__file__))
from utils import load_config

FASTA_DIR = "data/fasta"
REPORT_CSV = "results/sequence_match_report.csv"

os.makedirs("results", exist_ok=True)


def extract_chain_sequence(fasta_path: str, pdb_id: str, chain: str) -> str | None:
    """Parse a multi-chain FASTA and return the sequence for the given chain."""
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
            if f"CHAIN {chain.upper()}" in header_upper or f"CHAINS {chain.upper()}" in header_upper:
                return seq
            if f"{pdb_upper}_{chain.upper()}" in header_upper:
                return seq

    if len(sequences) == 1:
        return list(sequences.values())[0]

    for header, seq in sequences.items():
        if pdb_id.upper() in header.upper():
            return seq

    return None


def check_and_standardize(apo_seq: str, holo_seq: str, name: str) -> tuple[str, float, str]:
    """Align apo/holo sequences and return (canonical_seq, identity, strategy)."""
    if apo_seq == holo_seq:
        return apo_seq, 1.0, "identical"

    alignments = pairwise2.align.globalms(apo_seq, holo_seq, 2, -1, -5, -0.5)
    best = alignments[0]
    matched = sum(
        a == b for a, b in zip(best.seqA, best.seqB)
        if a != "-" and b != "-"
    )
    identity = matched / max(len(apo_seq), len(holo_seq))

    if identity >= 0.95:
        strategy = "high_identity_use_apo"
    else:
        strategy = "low_identity_review"
        print(f"[warn] {name}: low sequence identity {identity:.1%} — review manually")
        print(f"       apo: {len(apo_seq)} residues, holo: {len(holo_seq)} residues")

    return apo_seq, identity, strategy


def main():
    config = load_config()
    proteins = config["proteins"]

    report_rows = []

    for name, cfg in proteins.items():
        if cfg.get("skip_by_default"):
            print(f"[skip] {name} (skip_by_default)")
            continue

        chain = cfg.get("chain", "A")
        apo_pdb = cfg["apo_pdb"]
        holo_pdb = cfg["holo_pdb"]

        apo_fasta = os.path.join(FASTA_DIR, f"{apo_pdb}.fasta")
        holo_fasta = os.path.join(FASTA_DIR, f"{holo_pdb}.fasta")

        if not os.path.exists(apo_fasta) or not os.path.exists(holo_fasta):
            print(f"[error] {name}: missing FASTA file(s)")
            continue

        apo_seq = extract_chain_sequence(apo_fasta, apo_pdb, chain)
        holo_seq = extract_chain_sequence(holo_fasta, holo_pdb, chain)

        if not apo_seq or not holo_seq:
            print(f"[error] {name}: could not extract chain {chain} sequence")
            continue

        canonical, identity, strategy = check_and_standardize(apo_seq, holo_seq, name)

        canonical_path = os.path.join(FASTA_DIR, f"{name}_canonical.fasta")
        with open(canonical_path, "w") as f:
            f.write(f">{name}_canonical chain_{chain}\n")
            for i in range(0, len(canonical), 80):
                f.write(canonical[i : i + 80] + "\n")

        print(f"[ok] {name}: identity={identity:.1%}, strategy={strategy}, "
              f"apo={len(apo_seq)}, holo={len(holo_seq)}, canonical={len(canonical)}")

        report_rows.append({
            "protein": name,
            "apo_pdb": apo_pdb,
            "holo_pdb": holo_pdb,
            "apo_length": len(apo_seq),
            "holo_length": len(holo_seq),
            "identity_pct": round(identity * 100, 1),
            "strategy": strategy,
        })

    with open(REPORT_CSV, "w", newline="") as f:
        fieldnames = ["protein", "apo_pdb", "holo_pdb", "apo_length", "holo_length",
                       "identity_pct", "strategy"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(report_rows)

    print(f"\nSequence report saved to {REPORT_CSV}")


if __name__ == "__main__":
    main()
