#!/usr/bin/env python3
"""EX130 v2 Savio-side MCP client for pinned mcp==1.30.0.

Independent Researchstack approval and ClusterMCP one-use Slurm approval
are required. This client never fabricates approval, retries a dispatch,
or changes the Yipcluster research source snapshot.
"""
from __future__ import annotations

import argparse
import asyncio
import json
import os
import sys
import uuid
from pathlib import Path


def unpack(result):
    if getattr(result, "isError", False):
        raise RuntimeError("ClusterMCP rejected tool call: " + str(result))
    if getattr(result, "structuredContent", None) is not None:
        return result.structuredContent
    pieces = []
    for part in getattr(result, "content", []):
        if getattr(part, "text", None):
            try:
                pieces.append(json.loads(part.text))
            except ValueError:
                pieces.append(part.text)
    return pieces[0] if len(pieces) == 1 else pieces


async def invoke(args):
    url = os.environ.get("CLUSTERMCP_URL", "")
    token = os.environ.get("CLUSTERMCP_TOKEN", "")
    if not (url.startswith("https://") or url.startswith("http://127.0.0.1:")) or not token:
        raise ValueError("Set CLUSTERMCP_URL (HTTPS or localhost) and CLUSTERMCP_TOKEN; no embedded credentials")

    import httpx
    from mcp import ClientSession
    from mcp.client.streamable_http import streamable_http_client

    context = {
        "purpose": "EX130 v2 fixed-source within-gene isoform analysis on Yipcluster with independent scheduler approval",
        "turn_label": "ihec-flagship-ex130-v2",
    }
    timeout = httpx.Timeout(30.0, read=300.0)
    headers = {"Authorization": f"Bearer {token}"}
    # MCP SDK 1.30 removed headers= from streamable_http_client.
    # Authentication lives on the caller-owned HTTP client.
    async with httpx.AsyncClient(headers=headers, timeout=timeout) as http_client:
        async with streamable_http_client(url, http_client=http_client) as (rx, tx, _get_session_id):
            async with ClientSession(rx, tx) as session:
                await session.initialize()

                async def call(name, data):
                    answer = unpack(await session.call_tool(name, data))
                    print(json.dumps({"tool": name, "answer": answer}, indent=2))
                    return answer

                status = await call("slurm_gateway_status", {"context": context})
                if args.mode == "status":
                    return 0
                if args.mode == "propose":
                    if not args.plan_json or not args.script_file:
                        raise ValueError("--plan-json and --script-file are required for an operator-reviewed plan")
                    plan = json.loads(Path(args.plan_json).read_text())
                    if not isinstance(plan, dict):
                        raise ValueError("Researchstack plan must be a JSON object")
                    script = Path(args.script_file).read_text()
                    if not script.strip():
                        raise ValueError("Submitted script must not be empty")
                    await call("research_propose", {
                        "plan": plan, "script": script,
                        "request_id": str(uuid.uuid4()), "context": context,
                    })
                    return 0
                if not isinstance(status, dict) or not status.get("enabled") or not status.get("reachable"):
                    raise RuntimeError("Slurm broker disabled/unreachable. No submission attempted.")
                if not args.run_id or not args.plan_sha256:
                    raise ValueError("--run-id and --plan-sha256 required for an existing APPROVED plan")
                await call("slurm_prepare", {
                    "run_id": args.run_id, "plan_sha256": args.plan_sha256, "context": context,
                })
                if args.mode == "prepare":
                    return 0
                if not args.approval_id:
                    raise ValueError("--approval-id required; this script does not generate approvals")
                request_id = str(uuid.uuid4())
                print(f"DISPATCH_OPERATION_ID={request_id}", flush=True)
                try:
                    await call("slurm_submit", {
                        "run_id": args.run_id, "plan_sha256": args.plan_sha256,
                        "approval_id": args.approval_id, "request_id": request_id, "context": context,
                    })
                except Exception:
                    print(
                        f"UNCERTAIN_DISPATCH: inspect slurm_operation_status for {request_id}; "
                        "NEVER resubmit with a different UUID",
                        file=sys.stderr,
                    )
                    raise
                return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", choices=("status", "propose", "prepare", "submit-approved"))
    ap.add_argument("--plan-json")
    ap.add_argument("--script-file")
    ap.add_argument("--run-id")
    ap.add_argument("--plan-sha256")
    ap.add_argument("--approval-id")
    args = ap.parse_args()
    try:
        return asyncio.run(invoke(args))
    except Exception as e:
        print(type(e).__name__ + ": " + str(e), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
