#!/usr/bin/env bash
# Stage exact code only from saviocluster to a new Yipcluster scheduler work-root folder.
# SSH and rsync are staging channels, NOT ClusterMCP job submission.
set -Eeuo pipefail
: "${IHEC_SAVIO_CHECKOUT:?Set to exact clean checkout produced by bootstrap_savio_ex130_v2.sh}"
: "${IHEC_PINNED_SHA:?Set to exact 40-character GitHub commit SHA}"
[[ "$IHEC_PINNED_SHA" =~ ^[a-f0-9]{40}$ ]] || exit 2
[[ -d "$IHEC_SAVIO_CHECKOUT/.git" ]] || { echo 'Missing git checkout' >&2; exit 2; }
[[ "$(git -C "$IHEC_SAVIO_CHECKOUT" rev-parse HEAD)" == "$IHEC_PINNED_SHA" ]] || {
  echo 'Savio source is not the requested pinned commit' >&2; exit 2;
}
ARCHIVE="$IHEC_SAVIO_CHECKOUT/EX130/v2/isoform_switch_v2.tar.gz"
[[ -f "$ARCHIVE" ]] || exit 2
printf '%s  %s\n' '0aa9ed114e6f0c89077d8f2490f776c8f10b03152158cc93a6825edaa59bcdea' "$ARCHIVE" | sha256sum -c -
YIP_HOST="${YIP_HOST:-hchow@yipcluster}"
[[ "$YIP_HOST" =~ ^[a-zA-Z0-9_.-]+@[a-zA-Z0-9_.-]+$ ]] || { echo 'Unsafe YIP_HOST' >&2; exit 2; }
BASE='/autofs/login001_ssd/hchow/scratch/clustermcp-runs'
DEST="$BASE/IHEC_flagship_source_${IHEC_PINNED_SHA}"
# Refuse existing source tree; never adopt or overwrite a previous run or staging attempt.
ssh "$YIP_HOST" "umask 077; mkdir -p '$BASE'; test ! -e '$DEST'; mkdir '$DEST'; mkdir '$DEST/EX130'"
rsync -a -- "$IHEC_SAVIO_CHECKOUT/EX130/v2/" "$YIP_HOST:$DEST/EX130/v2/"
# Verify original GitHub-pinned archive on Yipcluster as well.
ssh "$YIP_HOST" "printf '%s  %s\\n' '0aa9ed114e6f0c89077d8f2490f776c8f10b03152158cc93a6825edaa59bcdea' '$DEST/EX130/v2/isoform_switch_v2.tar.gz' | sha256sum -c -"
printf 'YIP_VERSIONED_SOURCE_ROOT=%s\nIHEC_SOURCE_SHA=%s\n' "$DEST" "$IHEC_PINNED_SHA"
echo 'STAGED_CODE_ONLY=true; NO_SLURM_SUBMISSION_PERFORMED=true'
