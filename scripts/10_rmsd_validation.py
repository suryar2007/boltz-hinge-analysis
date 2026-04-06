"""RMSD validation using PyMOL align and rms_cur.

Computes three RMSD measurements for every protein comparison:
  1. `align`   – sequence-based Kabsch fit with iterative outlier rejection
  2. `align_all` – same fit, but RMSD reported over ALL matched atoms (cycles=0)
  3. `rms_cur` – RMSD of the pre-aligned coordinates saved by step 06 (no re-fitting)

Run as:  pymol -cq scripts/10_rmsd_validation.py
"""

from pymol import cmd, stored
import os, sys, csv, math

for _p in [
    os.path.join(os.getcwd(), "scripts"),
    os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else "",
]:
    if _p and _p not in sys.path:
        sys.path.insert(0, _p)

from utils import load_config

config = load_config()
OUT_CSV = "results/rmsd_validation.csv"
os.makedirs("results", exist_ok=True)

FIELDNAMES = [
    "protein", "comparison",
    "align_rmsd", "align_n_atoms", "align_rejected",
    "align_all_rmsd", "align_all_n_atoms",
    "rms_cur_rmsd", "rms_cur_n_atoms",
    "cealign_rmsd", "cealign_n_atoms",
]


def load_and_clean(path, obj_name, chain="A"):
    cmd.load(path)
    raw = cmd.get_object_list()[-1]
    if raw != obj_name:
        cmd.set_name(raw, obj_name)
    cmd.remove(f"{obj_name} and solvent")
    cmd.remove(f"{obj_name} and not polymer")
    cmd.remove(f"{obj_name} and not chain {chain}")


def ca_coords(obj_name):
    stored.coords = []
    cmd.iterate_state(1, f"{obj_name} and name CA", "stored.coords.append((x, y, z))")
    return stored.coords


def compute_rmsd_from_coords(c1, c2):
    n = min(len(c1), len(c2))
    if n == 0:
        return None, 0
    ssq = sum(
        (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2
        for a, b in zip(c1[:n], c2[:n])
    )
    return round(math.sqrt(ssq / n), 3), n


def run_align(mobile, target, cycles=5):
    """Run cmd.align and return (rmsd, n_aligned, n_total)."""
    try:
        result = cmd.align(f"{mobile} and name CA", f"{target} and name CA",
                           cycles=cycles, transform=1)
        rmsd_after = round(result[0], 3)
        n_after = result[1]
        n_before = result[4]
        return rmsd_after, n_after, n_before
    except Exception as e:
        print(f"  [error] align failed: {e}")
        return None, None, None


def run_cealign(mobile, target):
    """Run cmd.cealign and return (rmsd, n_aligned)."""
    try:
        result = cmd.cealign(target, mobile)
        return round(result["RMSD"], 3), result["alignment_length"]
    except Exception as e:
        print(f"  [error] cealign failed: {e}")
        return None, None


rows = []

for name, cfg in config["proteins"].items():
    if cfg.get("skip_by_default"):
        continue

    apo_pdb = cfg["apo_pdb"]
    holo_pdb = cfg["holo_pdb"]
    chain = cfg.get("chain", "A")
    holo_chain = cfg.get("holo_chain", chain)

    pairs = [
        ("boltz_apo_vs_true_apo",
         f"data/boltz_output/apo/{name}_apo_model_0.pdb", chain,
         f"data/pdb_cif/{apo_pdb}_true.pdb", chain,
         f"results/aligned_pdbs/{name}_boltz_apo_aligned.pdb"),
        ("boltz_holo_vs_true_holo",
         f"data/boltz_output/holo/{name}_holo_model_0.pdb", chain,
         f"data/pdb_cif/{holo_pdb}_true.pdb", holo_chain,
         f"results/aligned_pdbs/{name}_boltz_holo_aligned.pdb"),
        ("true_apo_vs_true_holo",
         f"data/pdb_cif/{apo_pdb}_true.pdb", chain,
         f"data/pdb_cif/{holo_pdb}_true.pdb", holo_chain,
         None),
        ("boltz_apo_vs_boltz_holo",
         f"data/boltz_output/apo/{name}_apo_model_0.pdb", chain,
         f"data/boltz_output/holo/{name}_holo_model_0.pdb", chain,
         None),
    ]

    for comp_label, mob_path, mob_chain, tgt_path, tgt_chain, aligned_path in pairs:
        if not os.path.exists(mob_path) or not os.path.exists(tgt_path):
            print(f"[skip] {name} {comp_label}: missing input files")
            continue

        print(f"\n{'='*60}")
        print(f"{name} — {comp_label}")
        print(f"{'='*60}")

        row = {"protein": name, "comparison": comp_label}

        # --- align with iterative rejection (default cycles=5, cutoff=2.0) ---
        cmd.reinitialize()
        load_and_clean(mob_path, "mobile", mob_chain)
        load_and_clean(tgt_path, "target", tgt_chain)
        n_mob = cmd.count_atoms("mobile and name CA")
        n_tgt = cmd.count_atoms("target and name CA")
        print(f"  CA atoms — mobile: {n_mob}, target: {n_tgt}")

        rmsd_ref, n_ref, n_total = run_align("mobile", "target", cycles=5)
        if rmsd_ref is not None:
            rejected = n_total - n_ref
            print(f"  align (refined):  RMSD = {rmsd_ref} A  ({n_ref} atoms, {rejected} rejected)")
            row["align_rmsd"] = rmsd_ref
            row["align_n_atoms"] = n_ref
            row["align_rejected"] = rejected

        # --- align with cycles=0 (no rejection — RMSD over all matched atoms) ---
        cmd.reinitialize()
        load_and_clean(mob_path, "mobile", mob_chain)
        load_and_clean(tgt_path, "target", tgt_chain)
        rmsd_all, n_all, _ = run_align("mobile", "target", cycles=0)
        if rmsd_all is not None:
            print(f"  align (all):      RMSD = {rmsd_all} A  ({n_all} atoms, no rejection)")
            row["align_all_rmsd"] = rmsd_all
            row["align_all_n_atoms"] = n_all

        # --- cealign for direct comparison with step 06 ---
        cmd.reinitialize()
        load_and_clean(mob_path, "mobile", mob_chain)
        load_and_clean(tgt_path, "target", tgt_chain)
        ce_rmsd, ce_n = run_cealign("mobile", "target")
        if ce_rmsd is not None:
            print(f"  cealign:          RMSD = {ce_rmsd} A  ({ce_n} atoms)")
            row["cealign_rmsd"] = ce_rmsd
            row["cealign_n_atoms"] = ce_n

        # --- rms_cur on pre-aligned coordinates (no re-fitting) ---
        if aligned_path and os.path.exists(aligned_path):
            cmd.reinitialize()
            load_and_clean(aligned_path, "aligned", "A")
            load_and_clean(tgt_path, "target", tgt_chain)
            c_aligned = ca_coords("aligned")
            c_target = ca_coords("target")
            rc_rmsd, rc_n = compute_rmsd_from_coords(c_aligned, c_target)
            if rc_rmsd is not None:
                print(f"  rms_cur (saved):  RMSD = {rc_rmsd} A  ({rc_n} CA atoms, no fitting)")
                row["rms_cur_rmsd"] = rc_rmsd
                row["rms_cur_n_atoms"] = rc_n
        else:
            print(f"  rms_cur: N/A (no pre-aligned PDB)")

        rows.append(row)

with open(OUT_CSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=FIELDNAMES, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)

print(f"\n{'='*60}")
print(f"Results written to {OUT_CSV}")
print(f"{'='*60}")

print("\n\n===== SUMMARY TABLE =====\n")
fmt = "{:<14s} {:<28s} {:>10s} {:>8s} {:>12s} {:>8s} {:>12s} {:>8s}"
print(fmt.format("Protein", "Comparison", "align", "n_atom", "align_all", "n_atom", "rms_cur", "n_atom"))
print("-" * 110)
for r in rows:
    print(fmt.format(
        r.get("protein", ""),
        r.get("comparison", ""),
        str(r.get("align_rmsd", "")),
        str(r.get("align_n_atoms", "")),
        str(r.get("align_all_rmsd", "")),
        str(r.get("align_all_n_atoms", "")),
        str(r.get("rms_cur_rmsd", "")),
        str(r.get("rms_cur_n_atoms", "")),
    ))
