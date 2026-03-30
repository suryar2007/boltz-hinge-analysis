#!/bin/bash
set -e

echo "=== 1. Fetch sequences ==="
python scripts/01_fetch_sequences.py

echo "=== 2. Check sequence matches ==="
python scripts/02_check_sequence_match.py

echo "=== 3. Generate Boltz YAMLs ==="
python scripts/03_make_yamls.py

echo "=== 4. Run Boltz predictions ==="
bash scripts/04_run_boltz.sh

echo "=== 5. Fetch true structures ==="
python scripts/05_fetch_true_structures.py

echo "=== 6. Align and compute RMSD ==="
pymol -cq scripts/06_align_and_rmsd.py

echo "=== 7. Contact maps ==="
python scripts/07_contact_maps.py

echo "=== 8. Render images ==="
pymol -cq scripts/08_render_images.py

echo "=== 9. DynDom ==="
bash scripts/09_run_dyndom.sh

echo ""
echo "=== Done. Check results/ ==="
