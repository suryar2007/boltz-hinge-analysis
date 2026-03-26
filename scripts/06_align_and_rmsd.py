"""Structural alignment and RMSD computation using PyMOL.

Run as: pymol -cq scripts/06_align_and_rmsd.py
"""

from pymol import cmd
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from utils import load_config, append_metric

config = load_config()
METRICS_CSV = "results/metrics.csv"

os.makedirs("results/aligned_pdbs", exist_ok=True)
os.makedirs("results/images", exist_ok=True)


def load_and_clean(path, obj_name, chain="A", remove_ligands=True):
    cmd.load(path)
    raw_name = cmd.get_object_list()[-1]
    cmd.set_name(raw_name, obj_name)
    cmd.remove(f"{obj_name} and solvent")
    if remove_ligands:
        cmd.remove(f"{obj_name} and not polymer")
    cmd.remove(f"{obj_name} and not chain {chain}")


def align_pair(mobile, target, label):
    try:
        result = cmd.cealign(target, mobile)
        rmsd = round(result["RMSD"], 3)
        n = result["alignment_length"]
        print(f"[ok] {label}: RMSD = {rmsd} Å ({n} atoms)")
        return rmsd, n
    except Exception as e:
        print(f"[error] {label}: {e}")
        return None, None


def render(img_path):
    cmd.show("cartoon", "all")
    cmd.hide("lines", "all")
    cmd.bg_color("white")
    cmd.set("ray_shadows", 0)
    cmd.zoom("all", 5)
    cmd.ray(1920, 1080)
    cmd.png(img_path, dpi=300)


for name, cfg in config["proteins"].items():
    if cfg.get("skip_by_default"):
        continue

    apo_pdb = cfg["apo_pdb"]
    holo_pdb = cfg["holo_pdb"]
    chain = cfg.get("chain", "A")

    boltz_apo = f"data/boltz_output/apo/{name}_apo_model_0.pdb"
    boltz_holo = f"data/boltz_output/holo/{name}_holo_model_0.pdb"
    true_apo = f"data/pdb_cif/{apo_pdb}_true.pdb"
    true_holo = f"data/pdb_cif/{holo_pdb}_true.pdb"

    # --- Alignment 1: Boltz apo vs True apo ---
    if os.path.exists(boltz_apo) and os.path.exists(true_apo):
        cmd.reinitialize()
        load_and_clean(boltz_apo, "pred_apo", chain)
        load_and_clean(true_apo, "true_apo", chain)
        cmd.color("cyan", "pred_apo")
        cmd.color("magenta", "true_apo")
        rmsd, n = align_pair("pred_apo", "true_apo", f"{name} boltz_apo vs true_apo")
        render(f"results/images/{name}_boltz_apo_vs_true_apo.png")
        cmd.save(f"results/aligned_pdbs/{name}_boltz_apo_aligned.pdb", "pred_apo")
        append_metric(METRICS_CSV, {
            "protein": name, "comparison": "boltz_apo_vs_true_apo",
            "rmsd": rmsd, "n_atoms": n,
        })

    # --- Alignment 2: Boltz holo vs True holo ---
    if os.path.exists(boltz_holo) and os.path.exists(true_holo):
        cmd.reinitialize()
        load_and_clean(boltz_holo, "pred_holo", chain)
        load_and_clean(true_holo, "true_holo", chain)
        cmd.color("cyan", "pred_holo")
        cmd.color("magenta", "true_holo")
        rmsd, n = align_pair("pred_holo", "true_holo", f"{name} boltz_holo vs true_holo")
        render(f"results/images/{name}_boltz_holo_vs_true_holo.png")
        cmd.save(f"results/aligned_pdbs/{name}_boltz_holo_aligned.pdb", "pred_holo")
        append_metric(METRICS_CSV, {
            "protein": name, "comparison": "boltz_holo_vs_true_holo",
            "rmsd": rmsd, "n_atoms": n,
        })

    # --- Alignment 3: True apo vs True holo (ground truth conformational change) ---
    if os.path.exists(true_apo) and os.path.exists(true_holo):
        cmd.reinitialize()
        load_and_clean(true_apo, "true_apo", chain)
        load_and_clean(true_holo, "true_holo", chain)
        cmd.color("cyan", "true_apo")
        cmd.color("magenta", "true_holo")
        rmsd, n = align_pair("true_apo", "true_holo", f"{name} true_apo vs true_holo")
        render(f"results/images/{name}_true_apo_vs_true_holo.png")
        append_metric(METRICS_CSV, {
            "protein": name, "comparison": "true_apo_vs_true_holo",
            "rmsd": rmsd, "n_atoms": n,
        })

    # --- Alignment 4: Boltz apo vs Boltz holo (did Boltz predict the motion?) ---
    if os.path.exists(boltz_apo) and os.path.exists(boltz_holo):
        cmd.reinitialize()
        load_and_clean(boltz_apo, "pred_apo", chain)
        load_and_clean(boltz_holo, "pred_holo", chain)
        cmd.color("cyan", "pred_apo")
        cmd.color("magenta", "pred_holo")
        rmsd, n = align_pair("pred_apo", "pred_holo", f"{name} boltz_apo vs boltz_holo")
        render(f"results/images/{name}_boltz_apo_vs_boltz_holo.png")
        append_metric(METRICS_CSV, {
            "protein": name, "comparison": "boltz_apo_vs_boltz_holo",
            "rmsd": rmsd, "n_atoms": n,
        })

print("Alignment complete. Results in results/metrics.csv")
