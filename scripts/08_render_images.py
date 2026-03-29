"""Render additional PyMOL images: pLDDT coloring and holo with ligand.

Run as: pymol -cq scripts/08_render_images.py
"""

from pymol import cmd
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from utils import load_config

config = load_config()
os.makedirs("results/images", exist_ok=True)


def render(path, width=1920, height=1080):
    cmd.bg_color("white")
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)
    cmd.zoom("all", 5)
    cmd.ray(width, height)
    cmd.png(path, dpi=300)


for name, cfg in config["proteins"].items():
    if cfg.get("skip_by_default"):
        continue
    chain = cfg.get("chain", "A")

    # pLDDT coloring on Boltz apo
    boltz_apo = f"data/boltz_output/apo/{name}_apo_model_0.pdb"
    if os.path.exists(boltz_apo):
        cmd.reinitialize()
        cmd.load(boltz_apo)
        raw = cmd.get_object_list()[0]
        cmd.set_name(raw, "pred")
        cmd.remove("pred and solvent")
        cmd.remove("pred and not chain " + chain)
        cmd.show("cartoon", "pred")
        cmd.spectrum("b", "blue_white_red", "pred")
        render(f"results/images/{name}_boltz_apo_plddt.png")
        print(f"[ok] {name} pLDDT image saved")

    # Boltz holo with ligand
    boltz_holo = f"data/boltz_output/holo/{name}_holo_model_0.pdb"
    if os.path.exists(boltz_holo):
        cmd.reinitialize()
        cmd.load(boltz_holo)
        raw = cmd.get_object_list()[0]
        cmd.set_name(raw, "pred_holo")
        cmd.remove("pred_holo and solvent")
        cmd.show("cartoon", "pred_holo and polymer")
        cmd.show("sticks", "pred_holo and organic")
        cmd.color("cyan", "pred_holo and polymer")
        cmd.color("yellow", "pred_holo and organic")
        render(f"results/images/{name}_boltz_holo_with_ligand.png")
        print(f"[ok] {name} holo+ligand image saved")

print("Image rendering complete.")
