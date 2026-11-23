#!/bin/bash
# Submit after ex130_screen.sbatch has finished successfully and task manifest exists.
set -euo pipefail
EX130_DIR="${EX130_DIR:-$HOME/IHEC_Sept_2025/EX130}"
MAX_CONCURRENT="${MAX_CONCURRENT:-4}"
[[ "$MAX_CONCURRENT" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid MAX_CONCURRENT" >&2; exit 2; }
cd "$EX130_DIR"
manifest="results/fixed_source/attribution_tasks.csv"
targets="results/fixed_source/attribution_targets.csv"
[[ -s "$manifest" && -s "$targets" ]] || { echo "Missing stage-1 manifests: run ex130_screen.sbatch first" >&2; exit 3; }
N="$(( $(wc -l < "$manifest") - 1 ))"
(( N > 0 )) || { echo "No attribution tasks" >&2; exit 4; }
# Do not silently launch an unbounded array: explicit concurrency throttle.
attr_job="$(sbatch --parsable --array="0-$((N-1))%${MAX_CONCURRENT}" ex130_attribution.sbatch)"
merge_job="$(sbatch --parsable --dependency="afterok:${attr_job}" ex130_summarize.sbatch)"
printf 'Submitted %s tasks, max %s concurrent\nAttribution array job: %s\nSummary depends on successful array: %s\n' "$N" "$MAX_CONCURRENT" "$attr_job" "$merge_job"
