#!/bin/bash
# DynDom analysis: compare apo/holo pairs for domain motion detection.
# Requires PDB format — Boltz outputs are already PDB.

set -e

if ! command -v dyndom &> /dev/null; then
    echo "[warn] DynDom not found. Install from https://dyndom.cmp.uea.ac.uk/dyndom/"
    echo "       Alternatively, use the web server at https://dyndom.cmp.uea.ac.uk/dyndom/domainselector.do"
    echo "       Upload: predicted_apo.pdb and predicted_holo.pdb for each protein"
    exit 0
fi

DYNDOM_OUT="results/dyndom"
mkdir -p "$DYNDOM_OUT"

python3 - <<'EOF'
import yaml
import os
import subprocess

with open("config.yaml") as f:
    config = yaml.safe_load(f)

for name, cfg in config["proteins"].items():
    if cfg.get("skip_by_default"):
        continue

    chain = cfg.get("chain", "A")

    pairs = [
        ("true_apo_vs_true_holo",
         f"data/pdb_cif/{cfg['apo_pdb']}_true.pdb",
         f"data/pdb_cif/{cfg['holo_pdb']}_true.pdb"),
        ("boltz_apo_vs_boltz_holo",
         f"data/boltz_output/apo/{name}_apo_model_0.pdb",
         f"data/boltz_output/holo/{name}_holo_model_0.pdb"),
    ]

    for label, pdb1, pdb2 in pairs:
        out_dir = f"results/dyndom/{name}_{label}"
        os.makedirs(out_dir, exist_ok=True)

        if not os.path.exists(pdb1) or not os.path.exists(pdb2):
            print(f"[skip] {name}/{label}: missing PDB file(s)")
            continue

        print(f"[dyndom] {name} — {label}")
        result = subprocess.run(
            ["dyndom", pdb1, pdb2],
            capture_output=True, text=True, cwd=out_dir
        )

        with open(f"{out_dir}/output.txt", "w") as f:
            f.write(result.stdout)
        if result.stderr:
            with open(f"{out_dir}/stderr.txt", "w") as f:
                f.write(result.stderr)

        for line in result.stdout.split("\n"):
            if any(k in line.lower() for k in ["hinge", "rotation", "domain", "angle"]):
                print(f"  {line.strip()}")

        print(f"[ok] {out_dir}/output.txt")
EOF
