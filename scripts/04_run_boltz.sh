#!/bin/bash
set -e

DEVICE="${BOLTZ_DEVICE:-gpu}"

run_boltz() {
    local yaml="$1"
    local out_dir="$2"
    local name=$(basename "$yaml" .yaml)
    local expected="${out_dir}/${name}_model_0.pdb"

    if [ -f "$expected" ]; then
        echo "[skip] $name — output already exists"
        return
    fi

    echo "[run] Boltz: $name"
    boltz predict --no_kernels "$yaml" \
        --out_dir "$out_dir" \
        --accelerator "$DEVICE" \
        --output_format pdb

    # Boltz writes to: {out_dir}/boltz_results_{name}/predictions/{name}/{name}_model_0.pdb
    local boltz_pdb="${out_dir}/boltz_results_${name}/predictions/${name}/${name}_model_0.pdb"
    if [ -f "$boltz_pdb" ]; then
        cp "$boltz_pdb" "$expected"
        echo "[ok] saved to $expected"
    else
        local fallback=$(find "${out_dir}/boltz_results_${name}" -name "*.pdb" 2>/dev/null | head -1)
        if [ -n "$fallback" ]; then
            cp "$fallback" "$expected"
            echo "[ok] saved to $expected (fallback)"
        else
            echo "[error] No PDB found for $name in ${out_dir}/boltz_results_${name}"
        fi
    fi
}

# Run apo predictions
for yaml in data/yaml/apo/*.yaml; do
    [ -f "$yaml" ] || continue
    run_boltz "$yaml" "data/boltz_output/apo"
done

# Run holo predictions
for yaml in data/yaml/holo/*.yaml; do
    [ -f "$yaml" ] || continue
    run_boltz "$yaml" "data/boltz_output/holo"
done

echo "Boltz predictions complete."
