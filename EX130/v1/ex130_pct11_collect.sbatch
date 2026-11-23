#!/usr/bin/env bash
#SBATCH --job-name=EX130_pct11_collect
#SBATCH --partition=normal
#SBATCH --nodelist=n001
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=ex130_pct11_collect_%j.out
#SBATCH --error=ex130_pct11_collect_%j.err
set -euo pipefail
EX130_DIR="${EX130_DIR:-$HOME/IHEC_Sept_2025/EX130}"
: "${PCT11_OUT_DIR:?must set PCT11_OUT_DIR in submission wrapper}"
cd "$EX130_DIR"
/usr/bin/python3 collect_ex130_pct11.py --ex130-dir "$EX130_DIR" --outdir "$PCT11_OUT_DIR"
