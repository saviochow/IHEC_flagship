#!/bin/bash
#SBATCH --job-name=EX130_track_collect
#SBATCH --partition=high
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=ex130_tracks_collect_%j.out
#SBATCH --error=ex130_tracks_collect_%j.err
set -euo pipefail
EX130_DIR="${EX130_DIR:-$HOME/IHEC_Sept_2025/EX130}"
TRACK_PYTHON_BIN="${TRACK_PYTHON_BIN:-$HOME/miniconda3/envs/DGL/bin/python3}"
cd "$EX130_DIR"
"$TRACK_PYTHON_BIN" collect_ex130_tracks.py --ex130-dir "$EX130_DIR"
