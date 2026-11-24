#!/usr/bin/env bash
# Run from Yipcluster EX130. Never overwrites original track PNGs/SVGs.
set -euo pipefail
EX130_DIR="${EX130_DIR:-$HOME/IHEC_Sept_2025/EX130}"
TRACK_PYTHON_BIN="${TRACK_PYTHON_BIN:-/home/hchow/miniconda3/envs/DGL/bin/python3}"
PCT11_MAX_CONCURRENT="${PCT11_MAX_CONCURRENT:-4}"
EX130_PCT11_BUILD_INDEX="${EX130_PCT11_BUILD_INDEX:-1}"
[[ "$PCT11_MAX_CONCURRENT" =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid PCT11_MAX_CONCURRENT' >&2;exit 2; }
((PCT11_MAX_CONCURRENT<=8)) || { echo 'PCT11_MAX_CONCURRENT must be <=8' >&2;exit 2; }
[[ "$EX130_PCT11_BUILD_INDEX" == 0 || "$EX130_PCT11_BUILD_INDEX" == 1 ]] || { echo 'Index cache flag must be 0 or 1' >&2;exit 2; }
cd "$EX130_DIR"
for f in seed_cases.csv track_sample_manifest.csv ex130_raw_tracks.py ex130_raw_tracks_pct11.py \
         verify_ex130_pct11_inputs.py ex130_pct11_prepare.sbatch \
         ex130_pct11_case.sbatch ex130_pct11_collect.sbatch collect_ex130_pct11.py;do
  [[ -s "$f" ]] || { echo "Missing $f" >&2;exit 2; }
done
[[ -x "$TRACK_PYTHON_BIN" ]] || { echo "DGL Python unavailable: $TRACK_PYTHON_BIN" >&2;exit 2; }
[[ -s results/raw_genomic_tracks/shared/gencode_selected_transcripts.json ]] || { echo 'Prepared GENCODE cache missing' >&2;exit 2; }
N=$(( $(wc -l < seed_cases.csv)-1 ))
((N==13)) || { echo "Expected 13 PILOT cases, found $N" >&2;exit 2; }
if [[ "$EX130_PCT11_BUILD_INDEX" == 0 ]];then unset EX130_PCT11_BUILD_INDEX; fi
RUN_LABEL="${EX130_PCT11_RUN_LABEL:-$(date -u +%Y%m%dT%H%M%SZ)}"
[[ "$RUN_LABEL" =~ ^[A-Za-z0-9_-]+$ ]] || { echo 'Unsafe run label' >&2;exit 2; }
PCT11_OUT_DIR="$EX130_DIR/results/raw_genomic_tracks_pct11/$RUN_LABEL"
if [[ -e "$PCT11_OUT_DIR" ]];then echo "Refusing overwrite: $PCT11_OUT_DIR" >&2;exit 2;fi
mkdir -p "$PCT11_OUT_DIR"
export EX130_DIR TRACK_PYTHON_BIN PCT11_OUT_DIR EX130_PCT11_BUILD_INDEX
# Use available n001 on normal. Existing attribution arrays and original plots are untouched.
prep=$(sbatch --parsable --partition=normal --nodelist=n001 --export=ALL ex130_pct11_prepare.sbatch)
array=$(sbatch --parsable --partition=normal --nodelist=n001 --dependency="afterok:$prep" \
    --array="0-$((N-1))%${PCT11_MAX_CONCURRENT}" --export=ALL ex130_pct11_case.sbatch)
collect=$(sbatch --parsable --partition=normal --nodelist=n001 --dependency="afterany:$array" \
    --export=ALL ex130_pct11_collect.sbatch)
cat > "$PCT11_OUT_DIR/submission_receipt.txt" <<EOF
EX130_PCT11_ROOT=$PCT11_OUT_DIR
EXPECTED_PILOTS=$N
WGBS_PERCENTAGE_COLUMN_ONE_BASED=11
WGBS_BETA_FORMULA=percentage/100
TRACK_PYTHON_BIN=$TRACK_PYTHON_BIN
NODE=n001
PARTITION=normal
PREP=$prep
ARRAY=$array
COLLECT=$collect
EOF
cat "$PCT11_OUT_DIR/submission_receipt.txt"
