"""Generate Cα contact maps and compute overlap metrics."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import PDB
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from utils import load_config, append_metric

CONTACT_THRESHOLD = 8.0  # Angstroms, Cα-Cα
METRICS_CSV = "results/metrics.csv"
os.makedirs("results/contact_maps", exist_ok=True)


def get_ca_coords(pdb_path: str, chain_id: str = "A") -> np.ndarray:
    """Extract Cα coordinates from a PDB file for a given chain."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("X", pdb_path)
    coords = []
    for model in structure:
        for chain in model:
            if chain.id != chain_id:
                continue
            for residue in chain:
                if PDB.is_aa(residue) and "CA" in residue:
                    coords.append(residue["CA"].get_vector().get_array())
        break  # first model only
    return np.array(coords)


def compute_contact_map(coords: np.ndarray, threshold: float = CONTACT_THRESHOLD) -> np.ndarray:
    """Compute binary Cα contact map."""
    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt((diff ** 2).sum(axis=-1))
    return (dist < threshold).astype(float)


def contact_map_overlap(map1: np.ndarray, map2: np.ndarray,
                        label: str = "") -> dict:
    """Compute precision, recall, F1 of map1 vs map2 (map2 = ground truth).

    Warns if sizes differ significantly — results may not be comparable.
    """
    n1, n2 = map1.shape[0], map2.shape[0]

    if abs(n1 - n2) > 10:
        print(f"[warn] {label}: residue count mismatch — pred={n1}, true={n2}. "
              f"Comparing only first {min(n1, n2)} residues. "
              f"Results are NOT directly comparable.")

    n = min(n1, n2)
    m1, m2 = map1[:n, :n], map2[:n, :n]

    mask = np.ones_like(m1, dtype=bool)
    for i in range(-2, 3):
        np.fill_diagonal(mask[max(0, -i):, max(0, i):], False)

    tp = ((m1 == 1) & (m2 == 1) & mask).sum()
    fp = ((m1 == 1) & (m2 == 0) & mask).sum()
    fn = ((m1 == 0) & (m2 == 1) & mask).sum()

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    return {
        "precision": round(float(precision), 3),
        "recall": round(float(recall), 3),
        "f1": round(float(f1), 3),
        "n_residues_compared": n,
        "pred_total": n1,
        "true_total": n2,
    }


def plot_contact_map(cmap: np.ndarray, title: str, out_path: str):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.imshow(cmap, cmap="Blues", origin="lower", interpolation="nearest")
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("Residue index")
    ax.set_ylabel("Residue index")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def plot_comparison(pred_map, true_map, title, out_path):
    """Overlay: TP=blue, FP=red, FN=gray."""
    n = min(pred_map.shape[0], true_map.shape[0])
    p, t = pred_map[:n, :n], true_map[:n, :n]

    rgb = np.ones((n, n, 3))
    rgb[(p == 1) & (t == 1)] = [0.2, 0.4, 0.8]  # TP: blue
    rgb[(p == 1) & (t == 0)] = [0.9, 0.2, 0.2]  # FP: red
    rgb[(p == 0) & (t == 1)] = [0.7, 0.7, 0.7]  # FN: gray

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.imshow(rgb, origin="lower")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Residue index")
    ax.set_ylabel("Residue index")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def main():
    config = load_config()

    for name, cfg in config["proteins"].items():
        if cfg.get("skip_by_default"):
            continue

        chain = cfg.get("chain", "A")
        holo_chain = cfg.get("holo_chain", chain)

        files = {
            "boltz_apo": f"data/boltz_output/apo/{name}_apo_model_0.pdb",
            "boltz_holo": f"data/boltz_output/holo/{name}_holo_model_0.pdb",
            "true_apo": f"data/pdb_cif/{cfg['apo_pdb']}_true.pdb",
            "true_holo": f"data/pdb_cif/{cfg['holo_pdb']}_true.pdb",
        }

        chain_for_label = {
            "boltz_apo": chain,
            "boltz_holo": chain,
            "true_apo": chain,
            "true_holo": holo_chain,
        }

        maps = {}
        for label, path in files.items():
            if not os.path.exists(path):
                print(f"[skip] {name}/{label}: file not found")
                continue
            coords = get_ca_coords(path, chain_id=chain_for_label[label])
            if len(coords) < 10:
                print(f"[warn] {name}/{label}: only {len(coords)} Cα atoms found, skipping")
                continue
            maps[label] = compute_contact_map(coords)
            plot_contact_map(
                maps[label], f"{name} — {label}",
                f"results/contact_maps/{name}_{label}.png",
            )
            print(f"[ok] {name}/{label}: contact map saved ({len(coords)} residues)")

        for pred_key, true_key in [("boltz_apo", "true_apo"), ("boltz_holo", "true_holo")]:
            if pred_key in maps and true_key in maps:
                overlap = contact_map_overlap(maps[pred_key], maps[true_key],
                                              label=f"{name} {pred_key}_vs_{true_key}")
                plot_comparison(
                    maps[pred_key], maps[true_key],
                    f"{name}: {pred_key} vs {true_key}",
                    f"results/contact_maps/{name}_{pred_key}_vs_{true_key}_comparison.png",
                )
                append_metric(METRICS_CSV, {
                    "protein": name,
                    "comparison": f"{pred_key}_vs_{true_key}_contacts",
                    "contact_precision": overlap["precision"],
                    "contact_recall": overlap["recall"],
                    "contact_f1": overlap["f1"],
                })
                print(f"[ok] {name} {pred_key} vs {true_key}: F1={overlap['f1']}")

    print("Contact maps complete.")


if __name__ == "__main__":
    main()
