#!/bin/bash
# Schedule one GFF/index-cache preparation, a throttled independent per-case
# high-partition array, then a terminal receipt collection (even on failures).
set -euo pipefail
EX130_DIR="${EX130_DIR:-$HOME/IHEC_Sept_2025/EX130}"
TRACK_PYTHON_BIN="${TRACK_PYTHON_BIN:-$HOME/miniconda3/envs/DGL/bin/python3}"
MAX_CONCURRENT="${MAX_CONCURRENT:-4}"
case "$MAX_CONCURRENT" in ''|*[!0-9]*) echo 'MAX_CONCURRENT must be positive integer' >&2;exit 2;; esac
if ((MAX_CONCURRENT<1 || MAX_CONCURRENT>16));then echo 'MAX_CONCURRENT must be 1..16' >&2;exit 2;fi
cd "$EX130_DIR"
[[ -f seed_cases.csv && -f ex130_raw_tracks.py ]] || { echo 'EX130 plotting files missing' >&2;exit 2; }
N=$("$TRACK_PYTHON_BIN" - <<'PY'
import csv
with open('seed_cases.csv',newline='') as f:
 rows=list(csv.DictReader(f))
if not rows or len({r['case_id'] for r in rows})!=len(rows):raise SystemExit('empty or duplicate seed cases')
print(len(rows))
PY
)
export EX130_DIR TRACK_PYTHON_BIN
export TRACK_MARKS="${TRACK_MARKS:-H3K36me3}" TRACK_SELECT="${TRACK_SELECT:-curated}"
export WGBS_BETA_COLUMN="${WGBS_BETA_COLUMN:-auto}" WGBS_SCALE="${WGBS_SCALE:-fraction}"
export BUILD_WGBS_INDEX_CACHE="${BUILD_WGBS_INDEX_CACHE:-0}" TABIX_BIN="${TABIX_BIN:-tabix}"
prep=$(sbatch --parsable ex130_tracks_prepare.sbatch)
array=$(sbatch --parsable --dependency="afterok:$prep" --array="0-$((N-1))%${MAX_CONCURRENT}" ex130_tracks.sbatch)
collect=$(sbatch --parsable --dependency="afterany:$array" ex130_tracks_collect.sbatch)
printf 'TRACK_PREP_JOB=%s\nTRACK_PARALLEL_ARRAY_JOB=%s\nTRACK_COLLECT_JOB=%s\nCASE_TASKS=%s\nPARALLEL_LIMIT=%s\n' "$prep" "$array" "$collect" "$N" "$MAX_CONCURRENT"
