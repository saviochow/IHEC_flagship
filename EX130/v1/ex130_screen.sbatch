#!/bin/bash
#SBATCH --job-name=EX130_screen
#SBATCH --partition=high
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=10:00:00
#SBATCH --output=ex130_screen_%j.out
#SBATCH --error=ex130_screen_%j.err
set -euo pipefail

EX130_DIR="${EX130_DIR:-$HOME/IHEC_Sept_2025/EX130}"
EX51_BASE="${EX51_BASE:-$HOME/IHEC_Mar_2024}"
PYTHON_BIN="${PYTHON_BIN:-$HOME/miniconda3/envs/DGL/bin/python3}"
FEATURE_ROOT="${FEATURE_ROOT:-/home/shared/IHEC/baseline_expt/features}"
SOURCE_MODEL="${SOURCE_MODEL:-IHECRE00001900_mRNA-Seq}"
cd "$EX130_DIR"
"$PYTHON_BIN" -V
mkdir -p results/diagonal results/fixed_source

# Stage 1a: all available held-out diagonal predictions; model varies by sample.
"$PYTHON_BIN" screen_ihec_cases.py screen --base "$EX51_BASE" \
    --outdir "$EX130_DIR/results/diagonal" --top-per-type 25
"$PYTHON_BIN" screen_ihec_cases.py features \
    --outdir "$EX130_DIR/results/diagonal" --features-root "$FEATURE_ROOT"

# Stage 1b: hold source GCN weights constant across evaluated contexts.
"$PYTHON_BIN" screen_ihec_cases.py screen --base "$EX51_BASE" \
    --trained-sample "$SOURCE_MODEL" \
    --outdir "$EX130_DIR/results/fixed_source" --top-per-type 25
"$PYTHON_BIN" screen_ihec_cases.py features \
    --outdir "$EX130_DIR/results/fixed_source" --features-root "$FEATURE_ROOT"

# Include representative auto candidates AND the user-approved EX51 pilot IDs.
"$PYTHON_BIN" make_attribution_manifest.py \
    --screen-dir "$EX130_DIR/results/fixed_source" \
    --logits-dir "$EX51_BASE/EX51/logits_debug_pc1" \
    --seed-cases "$EX130_DIR/seed_cases.csv" \
    --max-auto-candidates 20 \
    --outdir "$EX130_DIR/results/fixed_source"
"$PYTHON_BIN" export_attribution_features.py \
    --outdir "$EX130_DIR/results/fixed_source" --features-root "$FEATURE_ROOT"
"$PYTHON_BIN" feature_schema.py --out "$EX130_DIR/results/feature_index_schema.csv"

echo "READY: $(($(wc -l < results/fixed_source/attribution_tasks.csv) - 1)) attribution tasks"
echo "Review outputs and submit ex130_attribution.sbatch using submit_ex130_attribution.sh"
