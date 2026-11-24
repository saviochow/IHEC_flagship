#!/bin/bash
#SBATCH --job-name=EX130_attr
#SBATCH --partition=high
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=ex130_attr_%A_%a.out
#SBATCH --error=ex130_attr_%A_%a.err
# Array range is set by submit_ex130_attribution.sh after stage 1 completes.
set -euo pipefail

EX130_DIR="${EX130_DIR:-$HOME/IHEC_Sept_2025/EX130}"
EX51_BASE="${EX51_BASE:-$HOME/IHEC_Mar_2024}"
PYTHON_BIN="${PYTHON_BIN:-$HOME/miniconda3/envs/DGL/bin/python3}"
FEATURE_ROOT="${FEATURE_ROOT:-/home/shared/IHEC/baseline_expt/features}"
: "${SLURM_ARRAY_TASK_ID:?Must submit with sbatch --array=0-<N-1>%<limit>}"
cd "$EX130_DIR"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS=1
"$PYTHON_BIN" transcript_attribution.py \
    --base "$EX51_BASE" \
    --features-root "$FEATURE_ROOT" \
    --outdir "$EX130_DIR/results/fixed_source" \
    --task-index "$SLURM_ARRAY_TASK_ID" \
    --cpu-threads "${SLURM_CPUS_PER_TASK:-8}" \
    --skip-completed
