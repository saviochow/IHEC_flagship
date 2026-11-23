# EX130 v2 MCP client: SavioCluster setup and actual dispatch gates

**Updated 2026-11-23** after the real SavioCluster output reported `ModuleNotFoundError: No module named 'mcp'`. No Slurm job or Researchstack EX130 proposal was submitted.

The pinned EX130 v2 archive's `mcp_dispatch.py` uses the **1.x** low-level MCP Streamable-HTTP client interface (`streamable_http_client(..., headers=...)` and a 3-item transport tuple). MCP SDK **2.x** removed these interfaces, so installing `mcp==2.3.0` without migrating source would create a second failure. Use the exact maintenance-line client `mcp==1.30.0` from `requirements-mcp-client.txt` for this immutable archived script. GitHub CI now checks the imports, function signature, CLI startup, and credential failure behavior. The library version is distinct from the ClusterMCP server version.

On **saviocluster** using your already source-verified checkout:

```bash
set -euo pipefail
export IHEC_SAVIO_CHECKOUT="$HOME/Agent_projects/IHEC_flagship_v2_176a7aeb49ba"
python3 --version  # use Python >=3.10 for pinned SDK
CLIENT_ENV="$HOME/.venvs/ihec-ex130-mcp-1.30"
python3 -m venv "$CLIENT_ENV"
"$CLIENT_ENV/bin/python" -m pip install --disable-pip-version-check 'mcp==1.30.0'
"$CLIENT_ENV/bin/python" - <<'PY'
from importlib.metadata import version
from inspect import signature
from mcp import ClientSession
from mcp.client.streamable_http import streamable_http_client
assert version("mcp") == "1.30.0"
assert "headers" in signature(streamable_http_client).parameters
print("MCP_V1_CLIENT_IMPORT_OK")
PY
```

Next, perform a *read-only* status check with a trusted, explicitly configured `CLUSTERMCP_URL` and `CLUSTERMCP_TOKEN`:

```bash
: "${CLUSTERMCP_URL:?Configure actual MCP Streamable-HTTP endpoint}"
: "${CLUSTERMCP_TOKEN:?Acquire a valid scoped token through existing operator authentication}"
"$CLIENT_ENV/bin/python" "$IHEC_SAVIO_CHECKOUT/EX130/v2/mcp_dispatch.py" status
```

Do not paste tokens into GitHub, chat, shell history, or scripts. A ChatGPT connected plugin session does not give a local Python process an interchangeable bearer token.

**Separate production gate:** On 2026-11-23, direct ClusterMCP `slurm_gateway_status` returns `enabled=false`, `reachable=false`; the Researchstack run board contains only previously PROPOSED, non-executing Pines acceptance fixtures. No v2 EX130 plan is approved. Installing the SDK will not enable this gateway.

Do not run the previous examples with literal placeholders (`researchstack_plan.json`, `RUN_ID`, `PLAN_SHA256`, `ONE_USE_APPROVAL_ID`). A real plan must match the deployed Researchstack schema, refer to the already staged source under the Yipcluster allowlisted work root and a fresh isolated output root, pass independent review and be in APPROVED state. Then `slurm_prepare` returns the operator approval digest and one-shot approval is required for `slurm_submit`. Preserve the original dispatch UUID if the response is uncertain; never retry with a new UUID.

If normal operator policy cannot enable the gateway, do not bypass it by running unreviewed direct `sbatch` or by modifying the ClusterMCP policy in a research script. Independent direct Slurm can be used only as a separately approved operator fallback, with its own provenance.
