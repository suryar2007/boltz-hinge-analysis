#!/bin/bash
#SBATCH --job-name=boltz_hinge
#SBATCH --partition=mit_normal_gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/boltz_%j.out
#SBATCH --error=logs/boltz_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=suryar@mit.edu

set -euo pipefail

PROJECT_DIR="$HOME/repos/boltz-hinge-analysis"
START_STEP="${1:-1}"

mkdir -p "$PROJECT_DIR/logs"
cd "$PROJECT_DIR"

echo "============================================"
echo " Boltz Hinge Pipeline (steps 1-8)"
echo " Job ID:     $SLURM_JOB_ID"
echo " Node:       $SLURMD_NODENAME"
echo " Date:       $(date)"
echo " Start step: $START_STEP"
echo "============================================"

# ── GPU — let SLURM assign, don't hardcode ───────────────────────────────
# SLURM sets CUDA_VISIBLE_DEVICES automatically — do not override it
export BOLTZ_DEVICE="gpu"

# ── Conda ─────────────────────────────────────────────────────────────────
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate hinge_test
export PATH="$HOME/miniconda3/envs/hinge_test/bin:$PATH"
echo "[✓] Python:    $(python --version)"
echo "[✓] Conda env: hinge_test"

# ── Dependency check ──────────────────────────────────────────────────────
echo "[*] Checking dependencies..."
MISSING=0
check_cmd() { command -v "$1" &>/dev/null && echo "  [✓] $1" || { echo "  [✗] $1 (step $2)"; MISSING=1; }; }
check_py()  { python -c "import $1" 2>/dev/null && echo "  [✓] $1" || { echo "  [✗] $1 (step $2)"; MISSING=1; }; }

check_cmd pymol      "6 and 8"
check_py  boltz      "4"
check_py  Bio        "2 and 7"
check_py  yaml       "1 2 3"
check_py  numpy      "7"
check_py  matplotlib "7"
check_py  requests   "1 and 5"

[ "$MISSING" -eq 1 ] && { echo "[error] Fix deps above."; exit 1; }
echo "[✓] All dependencies found"

# ── GPU check ─────────────────────────────────────────────────────────────
echo "[*] GPU assigned by SLURM:"
nvidia-smi --query-gpu=index,name,memory.free \
    --format=csv,noheader,nounits | \
    awk -F',' '{printf "  GPU %s (%s): %s MiB free\n",$1,$2,$3}'
echo "[*] CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"

# ── Verify --no_kernels is set in 04_run_boltz.sh ─────────────────
if ! grep -q "no_kernels" scripts/04_run_boltz.sh; then
    echo "[warn] --no_kernels not found in scripts/04_run_boltz.sh"
    echo "       Patching now..."
    sed -i 's/boltz predict/boltz predict --no_kernels/' scripts/04_run_boltz.sh
    echo "[✓] Patched: $(grep 'boltz predict' scripts/04_run_boltz.sh)"
else
    echo "[✓] --no_kernels confirmed in 04_run_boltz.sh"
fi

# ── PyMOL wrapper (pymol -cq returns 0 even on script crash) ─────────────
run_pymol() {
    local script="$1"
    local tmpout
    tmpout=$(mktemp)
    PYTHONPATH="$PROJECT_DIR/scripts:${PYTHONPATH:-}" pymol -cq "$script" 2>&1 | tee "$tmpout"
    if grep -q "Traceback\|ModuleNotFoundError" "$tmpout"; then
        echo "[error] PyMOL script $script failed — see output above"
        rm -f "$tmpout"
        return 1
    fi
    rm -f "$tmpout"
    return 0
}

# ── Step runner ───────────────────────────────────────────────────────────
step() {
    local N=$1 LABEL=$2; shift 2
    if [ "$START_STEP" -le "$N" ]; then
        echo ""
        echo "--------------------------------------------"
        echo "=== Step $N: $LABEL"
        echo "    Started: $(date)"
        echo "--------------------------------------------"
        "$@"
        echo "    Finished: $(date)"
    else
        echo "[skip] Step $N ($LABEL)"
    fi
}

# ── Pipeline ──────────────────────────────────────────────────────────────
step 1 "Fetch sequences"         python scripts/01_fetch_sequences.py
step 2 "Check sequence matches"  python scripts/02_check_sequence_match.py
step 3 "Generate Boltz YAMLs"    python scripts/03_make_yamls.py
step 4 "Run Boltz predictions"   bash   scripts/04_run_boltz.sh
step 5 "Fetch true structures"   python scripts/05_fetch_true_structures.py
step 6 "Align and compute RMSD"  run_pymol scripts/06_align_and_rmsd.py
step 7 "Contact maps"            python scripts/07_contact_maps.py
step 8 "Render images"           run_pymol scripts/08_render_images.py

echo ""
echo "============================================"
echo " DONE — $(date)"
echo " Results: $PROJECT_DIR/results/"
echo " Logs:    $PROJECT_DIR/logs/"
echo "============================================"
