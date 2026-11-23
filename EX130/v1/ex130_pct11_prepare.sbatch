#!/usr/bin/env bash
#SBATCH --job-name=EX130_pct11_prep
#SBATCH --partition=normal
#SBATCH --nodelist=n001
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=ex130_pct11_prep_%j.out
#SBATCH --error=ex130_pct11_prep_%j.err
set -euo pipefail
EX130_DIR="${EX130_DIR:-$HOME/IHEC_Sept_2025/EX130}"
TRACK_PYTHON_BIN="${TRACK_PYTHON_BIN:-/home/hchow/miniconda3/envs/DGL/bin/python3}"
: "${PCT11_OUT_DIR:?must set PCT11_OUT_DIR in submission wrapper}"
cd "$EX130_DIR"
mkdir -p "$PCT11_OUT_DIR/shared"
"$TRACK_PYTHON_BIN" verify_ex130_pct11_inputs.py \
  --ex130-dir "$EX130_DIR" \
  --outdir "$PCT11_OUT_DIR/shared" \
  ${EX130_PCT11_BUILD_INDEX:+--build-index-cache}
"$TRACK_PYTHON_BIN" ex130_raw_tracks_pct11.py prepare \
  --marks H3K36me3 \
  --wgbs-beta-column 11 --wgbs-scale percent \
  --outdir "$PCT11_OUT_DIR/shared" \
  --annotation-cache-dir "$EX130_DIR/results/raw_genomic_tracks/shared" \
  --index-cache-dir "$PCT11_OUT_DIR/shared/wgbs_index_cache"
