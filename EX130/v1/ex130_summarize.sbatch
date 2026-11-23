#!/bin/bash
#SBATCH --job-name=EX130_merge
#SBATCH --partition=high
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=ex130_merge_%j.out
#SBATCH --error=ex130_merge_%j.err
set -euo pipefail
EX130_DIR="${EX130_DIR:-$HOME/IHEC_Sept_2025/EX130}"
PYTHON_BIN="${PYTHON_BIN:-$HOME/miniconda3/envs/DGL/bin/python3}"
cd "$EX130_DIR"
"$PYTHON_BIN" summarize_attribution.py --outdir "$EX130_DIR/results/fixed_source"

# Generate one 4-panel PNG per candidate for manual screening.
"$PYTHON_BIN" plot_candidate_panels.py \
    --outdir "$EX130_DIR/results/fixed_source" --max-cases 25
