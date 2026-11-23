#!/usr/bin/env bash
# Submit EX130 analyses across Yipcluster high/n003 and normal/n001+n002.
# Capacity-weighted deterministic assignment: 3:2:2 (192:128:128 cores).
# Does NOT rewrite task manifests or original .sbatch files.
set -euo pipefail

MODE="${1:-all}"          # attribution | tracks | all
ACTION="${2:---dry-run}"  # --dry-run | --submit
case "$MODE" in attribution|tracks|all) ;; *) echo 'Usage: bash submit_ex130_balanced.sh [attribution|tracks|all] [--dry-run|--submit]' >&2; exit 2;; esac
case "$ACTION" in --dry-run|--submit) ;; *) echo "Unknown action: $ACTION" >&2; exit 2;; esac

EX130_DIR="${EX130_DIR:-$HOME/IHEC_Sept_2025/EX130}"
PYTHON_BIN="${PYTHON_BIN:-/home/hchow/miniconda3/envs/DGL/bin/python3}"
TRACK_PYTHON_BIN="${TRACK_PYTHON_BIN:-$PYTHON_BIN}"
ATTR_MAX_CONCURRENT="${ATTR_MAX_CONCURRENT:-${MAX_CONCURRENT:-4}}"
TRACK_MAX_CONCURRENT="${TRACK_MAX_CONCURRENT:-${MAX_CONCURRENT:-4}}"
TRACK_MARKS="${TRACK_MARKS:-H3K36me3}"
TRACK_SELECT="${TRACK_SELECT:-curated}"
WGBS_BETA_COLUMN="${WGBS_BETA_COLUMN:-auto}"
WGBS_SCALE="${WGBS_SCALE:-fraction}"
BUILD_WGBS_INDEX_CACHE="${BUILD_WGBS_INDEX_CACHE:-0}"
TABIX_BIN="${TABIX_BIN:-tabix}"
export EX130_DIR PYTHON_BIN TRACK_PYTHON_BIN TRACK_MARKS TRACK_SELECT
export WGBS_BETA_COLUMN WGBS_SCALE BUILD_WGBS_INDEX_CACHE TABIX_BIN
cd "$EX130_DIR"

[[ "$ATTR_MAX_CONCURRENT" =~ ^[0-9]+$ && "$TRACK_MAX_CONCURRENT" =~ ^[0-9]+$ ]] || {
  echo 'Concurrency limits must be integers >= 3' >&2; exit 2;
}
(( ATTR_MAX_CONCURRENT >= 3 && TRACK_MAX_CONCURRENT >= 3 )) || {
  echo 'Concurrency limits must be >=3 to give each node at least one slot' >&2; exit 2;
}
(( ATTR_MAX_CONCURRENT <= 48 && TRACK_MAX_CONCURRENT <= 48 )) || {
  echo 'Concurrency >48 requires editing and reviewing resource admission explicitly' >&2; exit 2;
}

# CPU-proportional task assignment; each task is assigned exactly once.
# Index residues 0,1,2 -> n003; 3,4 -> n001; 5,6 -> n002.
make_arrays() {
  local n="$1" i residue
  ARR_N003='' ARR_N001='' ARR_N002=''
  for (( i=0; i<n; i++ )); do
    residue=$((i % 7))
    if (( residue < 3 )); then ARR_N003+="${ARR_N003:+,}${i}"
    elif (( residue < 5 )); then ARR_N001+="${ARR_N001:+,}${i}"
    else ARR_N002+="${ARR_N002:+,}${i}"
    fi
  done
  [[ -n "$ARR_N003" && -n "$ARR_N001" && -n "$ARR_N002" ]] || {
    echo "Too few tasks ($n) for three node pools" >&2; exit 2;
  }
  COUNT_N003=$(awk -F, '{print NF}' <<< "$ARR_N003")
  COUNT_N001=$(awk -F, '{print NF}' <<< "$ARR_N001")
  COUNT_N002=$(awk -F, '{print NF}' <<< "$ARR_N002")
  (( COUNT_N003 + COUNT_N001 + COUNT_N002 == n )) || {
    echo 'Internal assignment error' >&2; exit 2;
  }
}

# Spread total in-flight tasks over node pools using target 3:2:2 ratio.
# The slot allocator assigns each next slot to the currently most under-served pool.
make_caps() {
  local total="$1" s j best gap bestgap
  local -a weights=(3 2 2) slots=(1 1 1)
  for ((s=4; s<=total; s++)); do
    best=0
    bestgap=-999999
    for ((j=0; j<3; j++)); do
      gap=$((s * weights[j] - 7 * slots[j]))
      if (( gap > bestgap )); then bestgap=$gap; best=$j; fi
    done
    slots[best]=$((slots[best]+1))
  done
  CAP_N003="${slots[0]}" CAP_N001="${slots[1]}" CAP_N002="${slots[2]}"
}

# Verify that explicitly pinned nodes belong to their requested Slurm partitions.
check_partitions() {
  command -v sinfo >/dev/null || { echo 'sinfo unavailable: cannot validate partition/node membership' >&2; exit 3; }
  command -v sbatch >/dev/null || { echo 'sbatch unavailable' >&2; exit 3; }
  local partition node actual
  for pair in 'high:n003' 'normal:n001' 'normal:n002'; do
    partition="${pair%%:*}"; node="${pair#*:}"
    actual=$(sinfo -h -N -p "$partition" -o '%N')
    if ! grep -Fxq -- "$node" <<< "$actual"; then
      echo "Refusing submission: $node is not listed as a $partition node" >&2
      echo "Inspect: sinfo -N -p high,normal -o '%P %N %c %m %t'" >&2
      exit 3
    fi
  done
  echo 'SLURM_NODE_PARTITION_MEMBERSHIP_VERIFIED=true'
  sinfo -N -p high,normal -o '%P %N %c %m %t'
}

RECEIPT=''
log() { printf '%s\n' "$*"; if [[ -n "$RECEIPT" ]]; then printf '%s\n' "$*" >> "$RECEIPT"; fi; }
submit_job() {
  # The global LAST_JOB is assigned without subshell, preserving state and receipts.
  local label="$1" response; shift
  if [[ "$ACTION" == '--dry-run' ]]; then
    LAST_JOB="DRYRUN_${label}"
    log "DRYRUN $label: sbatch --parsable $*"
  else
    response=$(sbatch --parsable "$@")
    LAST_JOB="${response%%;*}"
    [[ "$LAST_JOB" =~ ^[0-9]+$ ]] || { echo "Unexpected sbatch response for $label: $response" >&2; exit 4; }
    log "$label=$LAST_JOB"
  fi
}

if [[ "$ACTION" == '--submit' ]]; then
  check_partitions
  [[ -x "$PYTHON_BIN" ]] || { echo "DGL interpreter unavailable: $PYTHON_BIN" >&2; exit 3; }
  [[ -x "$TRACK_PYTHON_BIN" ]] || { echo "Plot interpreter unavailable: $TRACK_PYTHON_BIN" >&2; exit 3; }
  # The earlier Python 3.8 error must be fixed BEFORE submitting all 69 tasks.
  if [[ "$MODE" == 'attribution' || "$MODE" == 'all' ]]; then
    if grep -Fq 'action=argparse.BooleanOptionalAction' transcript_attribution.py; then
      echo 'Python 3.8 incompatibility remains: run /usr/bin/python3 patch_ex130_python38.py first.' >&2
      exit 3
    fi
  fi
  mkdir -p results/balanced_submissions
  RECEIPT="results/balanced_submissions/$(date -u +'%Y%m%dT%H%M%SZ')_${MODE}_jobs.txt"
  : > "$RECEIPT"
fi

# Validate every requested workflow before the first sbatch submission.
if [[ "$MODE" == 'attribution' || "$MODE" == 'all' ]]; then
  [[ -s results/fixed_source/attribution_tasks.csv && -s results/fixed_source/attribution_targets.csv && -f ex130_attribution.sbatch && -f ex130_summarize.sbatch ]] || {
    echo 'Attribution input/worker files missing; refusing partial submission' >&2; exit 3;
  }
fi
if [[ "$MODE" == 'tracks' || "$MODE" == 'all' ]]; then
  [[ -s seed_cases.csv && -f ex130_raw_tracks.py && -f ex130_tracks_prepare.sbatch && -f ex130_tracks.sbatch && -f ex130_tracks_collect.sbatch ]] || {
    echo 'Genomic-track input/worker files missing; refusing partial submission' >&2; exit 3;
  }
fi

log "EX130_BALANCED_MODE=$MODE ACTION=$ACTION"
log 'WEIGHTED_NODE_CAPACITY=n003:192,n001:128,n002:128; job partition map high:n003,normal:n001,normal:n002'
log "ATTR_MAX_CONCURRENT=$ATTR_MAX_CONCURRENT TRACK_MAX_CONCURRENT=$TRACK_MAX_CONCURRENT"

if [[ "$MODE" == 'attribution' || "$MODE" == 'all' ]]; then
  [[ -s results/fixed_source/attribution_tasks.csv && -s results/fixed_source/attribution_targets.csv ]] || {
    echo 'Missing attribution manifests: do not rerun EX130 screening unless needed' >&2; exit 3;
  }
  [[ -f ex130_attribution.sbatch && -f ex130_summarize.sbatch ]] || { echo 'Missing attribution batch scripts' >&2; exit 3; }
  N=$(( $(wc -l < results/fixed_source/attribution_tasks.csv) - 1 ))
  (( N > 2 )) || { echo 'Attribution task manifest has fewer than 3 data rows' >&2; exit 3; }
  make_arrays "$N"; make_caps "$ATTR_MAX_CONCURRENT"
  log "ATTR_DISTRIBUTION total=$N high:n003=$COUNT_N003 normal:n001=$COUNT_N001 normal:n002=$COUNT_N002 in_flight=$CAP_N003,$CAP_N001,$CAP_N002"
  submit_job ATTR_HIGH --partition=high --nodelist=n003 --nodes=1 --job-name=EX130_A_H --array="${ARR_N003}%${CAP_N003}" --export=ALL ex130_attribution.sbatch
  attr_high="$LAST_JOB"
  submit_job ATTR_NORMAL_N001 --partition=normal --nodelist=n001 --nodes=1 --job-name=EX130_A_N1 --array="${ARR_N001}%${CAP_N001}" --export=ALL ex130_attribution.sbatch
  attr_n001="$LAST_JOB"
  submit_job ATTR_NORMAL_N002 --partition=normal --nodelist=n002 --nodes=1 --job-name=EX130_A_N2 --array="${ARR_N002}%${CAP_N002}" --export=ALL ex130_attribution.sbatch
  attr_n002="$LAST_JOB"
  # Summary runs only once, after ALL three attribution arrays succeed.
  submit_job ATTR_SUMMARY --partition=normal --nodelist=n001 --nodes=1 --job-name=EX130_A_SUM --kill-on-invalid-dep=yes --dependency="afterok:${attr_high}:${attr_n001}:${attr_n002}" --export=ALL ex130_summarize.sbatch
  log "ATTR_SUMMARY_DEPENDS_ON=$attr_high,$attr_n001,$attr_n002"
fi

if [[ "$MODE" == 'tracks' || "$MODE" == 'all' ]]; then
  [[ -s seed_cases.csv && -f ex130_raw_tracks.py && -f ex130_tracks_prepare.sbatch && -f ex130_tracks.sbatch && -f ex130_tracks_collect.sbatch ]] || {
    echo 'Missing raw-track seed cases or batch scripts' >&2; exit 3;
  }
  # seed_cases.csv is generated with a header and one row per case; no DGL imports on login.
  N=$(( $(wc -l < seed_cases.csv) - 1 ))
  (( N > 2 )) || { echo 'Track cases have fewer than 3 data rows' >&2; exit 3; }
  make_arrays "$N"; make_caps "$TRACK_MAX_CONCURRENT"
  log "TRACK_DISTRIBUTION total=$N high:n003=$COUNT_N003 normal:n001=$COUNT_N001 normal:n002=$COUNT_N002 in_flight=$CAP_N003,$CAP_N001,$CAP_N002"
  # Exactly one preparation job writes the shared GFF cache and optional tabix sidecars.
  submit_job TRACK_PREP --partition=normal --nodelist=n001 --nodes=1 --job-name=EX130_T_PREP --export=ALL ex130_tracks_prepare.sbatch
  track_prep="$LAST_JOB"
  submit_job TRACK_HIGH --partition=high --nodelist=n003 --nodes=1 --job-name=EX130_T_H --kill-on-invalid-dep=yes --array="${ARR_N003}%${CAP_N003}" --dependency="afterok:${track_prep}" --export=ALL ex130_tracks.sbatch
  track_high="$LAST_JOB"
  submit_job TRACK_NORMAL_N001 --partition=normal --nodelist=n001 --nodes=1 --job-name=EX130_T_N1 --kill-on-invalid-dep=yes --array="${ARR_N001}%${CAP_N001}" --dependency="afterok:${track_prep}" --export=ALL ex130_tracks.sbatch
  track_n001="$LAST_JOB"
  submit_job TRACK_NORMAL_N002 --partition=normal --nodelist=n002 --nodes=1 --job-name=EX130_T_N2 --kill-on-invalid-dep=yes --array="${ARR_N002}%${CAP_N002}" --dependency="afterok:${track_prep}" --export=ALL ex130_tracks.sbatch
  track_n002="$LAST_JOB"
  # Collector executes after all three arrays have reached terminal state.
  submit_job TRACK_COLLECT --partition=normal --nodelist=n002 --nodes=1 --job-name=EX130_T_COL --kill-on-invalid-dep=yes --dependency="afterany:${track_high}:${track_n001}:${track_n002}" --export=ALL ex130_tracks_collect.sbatch
  log "TRACK_COLLECT_DEPENDS_ON=$track_high,$track_n001,$track_n002"
fi

log "SUBMISSION_RECEIPT=${RECEIPT:-DRYRUN_ONLY}"
log 'NOTE: CPU capacity weights determine task placement, but Slurm memory/time admission still governs simultaneous execution.'
