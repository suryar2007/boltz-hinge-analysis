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
    boltz predict "$yaml" \
        --out_dir "$out_dir" \
        --accelerator "$DEVICE" \
        --output_format pdb

    # Rename to clean convention
    local generated=$(find "$out_dir" -name "*.pdb" -newer "$yaml" | head -1)
    if [ -n "$generated" ]; then
        mv "$generated" "$expected"
        echo "[ok] saved to $expected"
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
