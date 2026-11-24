#!/bin/bash
#SBATCH --job-name=EX130_track_case
#SBATCH --partition=high
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=ex130_tracks_%A_%a.out
#SBATCH --error=ex130_tracks_%A_%a.err
set -euo pipefail
EX130_DIR="${EX130_DIR:-$HOME/IHEC_Sept_2025/EX130}"
TRACK_PYTHON_BIN="${TRACK_PYTHON_BIN:-$HOME/miniconda3/envs/DGL/bin/python3}"
: "${SLURM_ARRAY_TASK_ID:?submit using submit_ex130_tracks.sh}"
cd "$EX130_DIR"
CASE_ID=$("$TRACK_PYTHON_BIN" - <<'PY'
import csv,os
with open('seed_cases.csv',newline='') as handle:
    rows=list(csv.DictReader(handle))
i=int(os.environ['SLURM_ARRAY_TASK_ID'])
if i<0 or i>=len(rows):raise SystemExit(f'out of bounds track array index: {i}')
print(rows[i]['case_id'])
PY
)
OUT="$EX130_DIR/results/raw_genomic_tracks"
CASE_OUT="$OUT/cases/$CASE_ID"
mkdir -p "$CASE_OUT"
"$TRACK_PYTHON_BIN" ex130_raw_tracks.py plot \
  --only-cases "$CASE_ID" \
  --select "${TRACK_SELECT:-curated}" \
  --marks "${TRACK_MARKS:-H3K36me3}" \
  --wgbs-beta-column "${WGBS_BETA_COLUMN:-auto}" \
  --wgbs-scale "${WGBS_SCALE:-fraction}" \
  --outdir "$CASE_OUT" \
  --annotation-cache-dir "$OUT/shared" \
  --index-cache-dir "$OUT/shared/wgbs_index_cache"
