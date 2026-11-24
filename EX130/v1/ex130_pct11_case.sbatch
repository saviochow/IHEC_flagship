#!/usr/bin/env bash
#SBATCH --job-name=EX130_pct11_plot
#SBATCH --partition=normal
#SBATCH --nodelist=n001
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=ex130_pct11_%A_%a.out
#SBATCH --error=ex130_pct11_%A_%a.err
set -euo pipefail
EX130_DIR="${EX130_DIR:-$HOME/IHEC_Sept_2025/EX130}"
TRACK_PYTHON_BIN="${TRACK_PYTHON_BIN:-/home/hchow/miniconda3/envs/DGL/bin/python3}"
: "${PCT11_OUT_DIR:?must set PCT11_OUT_DIR in submission wrapper}"
: "${SLURM_ARRAY_TASK_ID:?array index is missing}"
cd "$EX130_DIR"
CASE_ID=$(/usr/bin/python3 - <<'PY'
import csv,os
with open('seed_cases.csv', newline='') as file:
    cases=list(csv.DictReader(file))
i=int(os.environ['SLURM_ARRAY_TASK_ID'])
if i<0 or i>=len(cases):raise SystemExit('Invalid SLURM_ARRAY_TASK_ID {}'.format(i))
print(cases[i]['case_id'])
PY
)
CASE_OUT="$PCT11_OUT_DIR/cases/$CASE_ID"
mkdir -p "$CASE_OUT"
"$TRACK_PYTHON_BIN" ex130_raw_tracks_pct11.py plot \
  --only-cases "$CASE_ID" \
  --select curated \
  --marks H3K36me3 --signal pval \
  --wgbs-beta-column 11 --wgbs-scale percent \
  --outdir "$CASE_OUT" \
  --annotation-cache-dir "$EX130_DIR/results/raw_genomic_tracks/shared" \
  --index-cache-dir "$PCT11_OUT_DIR/shared/wgbs_index_cache"
# Fail the task if it omitted the beta track because of an incorrect BED column.
/usr/bin/python3 - "$CASE_OUT/track_measurement_status.csv" <<'PY'
import csv,sys
with open(sys.argv[1],newline='') as f:rows=list(csv.DictReader(f))
wgbs=[r for r in rows if r['track']=='WGBS CpG beta']
valid=[r for r in wgbs if r['status']=='USER_SPECIFIED_COL_11' and int(r['points'])>0]
if not wgbs or not valid:
    raise SystemExit('Missing verified CpG beta track for the case')
invalid=[r for r in wgbs if 'BETA_COLUMN' in r['status'] or r['status']=='NO_VALID_BETA_POINTS']
if invalid:raise SystemExit('Invalid WGBS CpG beta column conversion: {}'.format(invalid))
print('PCT11_CASE_VALIDATED WGBS_samples={} valid={}'.format(len(wgbs),len(valid)))
PY
