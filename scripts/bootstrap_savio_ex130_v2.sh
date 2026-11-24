#!/usr/bin/env bash
# Run on saviocluster. Installs no packages, modifies no existing checkout and submits no jobs.
set -Eeuo pipefail
REPO='https://github.com/saviochow/IHEC_flagship.git'
: "${IHEC_PINNED_SHA:?Set IHEC_PINNED_SHA to the exact published review-branch commit}"
[[ "$IHEC_PINNED_SHA" =~ ^[a-f0-9]{40}$ ]] || { echo 'IHEC_PINNED_SHA must be 40 hex digits' >&2; exit 2; }
DEST="${IHEC_DEST:-$HOME/Agent_projects/IHEC_flagship_v2_${IHEC_PINNED_SHA:0:12}}"
[[ ! -e "$DEST" ]] || { echo "Destination already exists: $DEST" >&2; exit 2; }
mkdir -p "$(dirname "$DEST")"
git clone --no-checkout "$REPO" "$DEST"
git -C "$DEST" fetch --depth 1 origin "$IHEC_PINNED_SHA"
[[ "$(git -C "$DEST" rev-parse FETCH_HEAD)" == "$IHEC_PINNED_SHA" ]] || { echo 'GitHub source identity mismatch' >&2; exit 3; }
git -C "$DEST" checkout --detach "$IHEC_PINNED_SHA"
ARCHIVE="$DEST/EX130/v2/isoform_switch_v2.tar.gz"
BUNDLES=("$DEST"/EX130/v2/bundle/[0-9][0-9].b64)
[[ ${#BUNDLES[@]} -eq 11 ]] || { echo "Expected 11 bundle chunks, got ${#BUNDLES[@]}" >&2; exit 3; }
cat "${BUNDLES[@]}" | base64 -d > "$ARCHIVE"
echo '0aa9ed114e6f0c89077d8f2490f776c8f10b03152158cc93a6825edaa59bcdea  '"$ARCHIVE" | sha256sum -c -
# Verify archive file names before extracting; never execute shell code from the archive at this step.
tar -tzf "$ARCHIVE" | while IFS= read -r path; do
  case "$path" in EX130/v2/*) ;; *) echo "Unexpected archive path: $path" >&2; exit 4 ;; esac
done
tar -xzf "$ARCHIVE" --no-same-owner -C "$DEST"
python3 -m compileall -q "$DEST/EX130/v2"
python3 -m unittest discover -s "$DEST/EX130/v2/tests" -v
printf '\nIHEC_EX130_V2_SOURCE_SHA=%s\nSAVIO_CHECKOUT=%s\n' "$IHEC_PINNED_SHA" "$DEST"
# MCP workflow is separately approval-gated. This is strictly a read-only capability probe.
if [[ "${IHEC_MCP_STATUS_CHECK:-0}" == 1 ]]; then
  python3 "$DEST/EX130/v2/mcp_dispatch.py" status
fi
