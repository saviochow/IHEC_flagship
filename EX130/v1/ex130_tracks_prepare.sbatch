#!/bin/bash
#SBATCH --job-name=EX130_track_prep
#SBATCH --partition=high
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=ex130_tracks_prepare_%j.out
#SBATCH --error=ex130_tracks_prepare_%j.err
set -euo pipefail
EX130_DIR="${EX130_DIR:-$HOME/IHEC_Sept_2025/EX130}"
TRACK_PYTHON_BIN="${TRACK_PYTHON_BIN:-$HOME/miniconda3/envs/DGL/bin/python3}"
OUT="$EX130_DIR/results/raw_genomic_tracks"
cd "$EX130_DIR"
mkdir -p "$OUT/shared"
"$TRACK_PYTHON_BIN" - <<'PY'
import numpy, matplotlib, pyBigWig, pysam
print('track plotting dependencies import successfully')
PY
extra=()
if [[ "${BUILD_WGBS_INDEX_CACHE:-0}" == "1" ]]; then
  extra+=(--build-index-cache --tabix-bin "${TABIX_BIN:-tabix}")
fi
"$TRACK_PYTHON_BIN" ex130_raw_tracks.py prepare \
  --marks "${TRACK_MARKS:-H3K36me3}" \
  --outdir "$OUT/shared" \
  --annotation-cache-dir "$OUT/shared" \
  --index-cache-dir "$OUT/shared/wgbs_index_cache" \
  "${extra[@]}"
