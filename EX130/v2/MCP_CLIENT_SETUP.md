# SavioCluster MCP client compatibility and execution

**Repaired on 2026-11-23 following live operator evidence and GitHub CI failure.**

The original immutable EX130 v2 scientific archive was built with the old MCP `streamable_http_client(url, headers=...)` interface. As of `mcp==1.30.0`, this call is no longer accepted, even though the package retains the three-item return tuple. Installing the Python package alone would therefore not repair the archived client.

**Use the separately versioned, source-reviewed replacement `scripts/mcp_dispatch_fixed.py`** committed in this GitHub branch. It passes authentication through a caller-owned `httpx.AsyncClient(headers=..., timeout=...)` and passes that client as `streamable_http_client(url, http_client=http_client)`, compatible with `mcp==1.30.0`. The archived scientific code, the exact SHA used for the Yipcluster staging, the existing output directories and Slurm scripts are unchanged.

## SavioCluster: isolated MCP client environment

The following has no scheduler side effects, and does not modify your pinned scientific checkout.

```bash
set -euo pipefail
VENV="$HOME/.venvs/ihec-ex130-mcp-1.30"
python3 -m venv "$VENV"
"$VENV/bin/python" -m pip install --disable-pip-version-check 'mcp==1.30.0'
"$VENV/bin/python" - <<'PY'
from importlib.metadata import version
from inspect import signature
from mcp import ClientSession
from mcp.client.streamable_http import streamable_http_client
assert version("mcp") == "1.30.0"
assert "http_client" in signature(streamable_http_client).parameters
assert "headers" not in signature(streamable_http_client).parameters
print("EX130_MCP_1_30_TRANSPORT_COMPATIBLE=true")
PY
```

From the already checked out pinned research source, fetch the new review-branch fix by exact commit SHA and export *only* the client script to a standalone path. Do not reset/rebase the original science worktree.

```bash
export IHEC_SAVIO_CHECKOUT="$HOME/Agent_projects/IHEC_flagship_v2_176a7aeb49ba"
export EX130_MCP_FIX_SHA='<40-hex commit of the verified compatibility fix>'
git -C "$IHEC_SAVIO_CHECKOUT" fetch --depth=1 origin "$EX130_MCP_FIX_SHA"
test "$(git -C "$IHEC_SAVIO_CHECKOUT" rev-parse FETCH_HEAD)" = "$EX130_MCP_FIX_SHA"
mkdir -p "$HOME/.local/libexec"
git -C "$IHEC_SAVIO_CHECKOUT" show "$EX130_MCP_FIX_SHA:scripts/mcp_dispatch_fixed.py" > "$HOME/.local/libexec/ihec_ex130_mcp_dispatch.py"
"$VENV/bin/python" "$HOME/.local/libexec/ihec_ex130_mcp_dispatch.py" --help
```

Configure an actual scoped bearer token and MCP Streamable-HTTP endpoint through existing operator-approved authentication. A ChatGPT connected plugin's token is not automatically shared with the Savio Python environment. Keep secrets out of logs and shell history.

```bash
: "${CLUSTERMCP_URL:?Set actual authenticated MCP URL}"
: "${CLUSTERMCP_TOKEN:?Set scoped token through operator vault}"
"$VENV/bin/python" "$HOME/.local/libexec/ihec_ex130_mcp_dispatch.py" status
```

**Current independent readback:** ClusterMCP's Yipcluster Slurm gateway is still `enabled=false` and `reachable=false`, and the Researchstack board contains no EX130 v2 approved plan. The fix permits client import and authentication but **cannot** enable the broker or approve/dispatch science jobs.

Do not call propose/prepare/submit with dummy placeholder IDs. Use the deployed Researchstack plan schema, versioned code under the authorized Yipcluster work root, APPROVED ledger state, and one-use operator approval. On uncertain dispatch, inspect the original operation UUID; never retry with a new one.
