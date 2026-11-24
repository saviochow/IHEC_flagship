#!/usr/bin/env python3
"""True per-transcript input gradients for the existing EX51 7-layer DGL GCN.

For a held-out transcript t, compute d(model(X)[t,0])/dX, where X contains
*every* graph node's features. This is fundamentally different from EX88's
Captum Saliency(model).attribute(X) followed by mean over masked nodes.

Records two distinct kinds of sensitivity:
  self: d y_t / d X[t, feature] (the transcript node's own input)
  all_nodes: sum_n |d y_t / d X[n, feature]| (including graph context)
Both are descriptive local model sensitivities, NOT contributions or causality.

Only writes under --outdir. Loads exact EX51 model/graph and verifies prediction,
label, chromosome masking, and WGBS-first feature order for identifiable cases.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from feature_schema import ASSAYS, REGIONS, mapping, write_mapping


def sha256_stream(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            b = fh.read(1024 * 1024)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def target_node_gradient(model, x, node_index: int):
    """Differentiate a single scalar node output, never sum output nodes.

    Works with torch tensors and models whose output shape is (N,) or (N,1).
    Exposed separately for dependency-light regression testing.
    """
    import torch
    if not x.requires_grad:
        raise ValueError("Input features must require gradients")
    output = model(x).reshape(-1)
    if node_index < 0 or node_index >= output.numel():
        raise IndexError(f"Invalid output node index {node_index}")
    grad = torch.autograd.grad(output[node_index], x, retain_graph=False, create_graph=False)[0]
    return output[node_index].detach(), grad.detach()


def find_wgbs(features_root: Path, base_id: str) -> Path:
    matches = sorted((features_root / base_id).glob(f"{base_id}.[0-9]*.standard-WGBS.cpg.all.txt"))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected ONE original WGBS file for {base_id}; found {len(matches)}")
    return matches[0]


def verify_wgbs_first(features_root: Path, base_id: str, txid_to_node: dict, x, *, atol: float = 1e-4):
    """Audit model tensor columns 0..20 against raw WGBS for selected txids.

    EX18 applies no transformation to WGBS; missing beta is represented by -1.
    A matching row establishes the first 21 graph inputs belong to DNAm,
    independently of the external RF reference used by EX88 plotting.
    """
    path = find_wgbs(features_root, base_id)
    cols = ["beta_hg38" + v for v in REGIONS]
    header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    missing = set(cols + ["transcript.id"]) - set(header)
    if missing:
        raise ValueError(f"Missing raw WGBS fields in {path}: {sorted(missing)}")
    waiting = set(txid_to_node)
    hits, checked = 0, []
    for chunk in pd.read_csv(path, sep="\t", usecols=["transcript.id"] + cols,
                             dtype={"transcript.id": str}, chunksize=25000):
        seen = chunk.loc[chunk["transcript.id"].isin(waiting)]
        for r in seen.itertuples(index=False, name=None):
            record = dict(zip(["transcript.id"] + cols, r))
            txid = record.pop("transcript.id")
            node = int(txid_to_node[txid])
            original = pd.to_numeric(pd.Series(record), errors="coerce").fillna(-1).to_numpy(dtype="float32")
            stored = x[node, :len(REGIONS)].detach().cpu().numpy()
            if not np.allclose(original, stored, atol=atol, rtol=1e-4):
                diff = np.abs(original - stored)
                k = int(np.argmax(diff))
                raise ValueError(
                    f"EX18 FEATURE ORDER MISMATCH: {txid}, feature index {k}, "
                    f"raw {cols[k]}={original[k]}, graph={stored[k]}, difference={diff[k]}. "
                    "Stop: do not label saliency assays without resolving graph provenance."
                )
            hits += 1
            checked.append(txid)
            waiting.discard(txid)
        if not waiting:
            break
    if hits != len(txid_to_node):
        raise ValueError(f"Incomplete raw WGBS audit for {base_id}: matched {hits}/{len(txid_to_node)} targets. "
                         f"Missing: {sorted(waiting)[:10]}. Exclude unmeasured candidates before claiming the feature interpretation.")
    return {"wgbs_file": str(path), "n_verified_transcripts": hits,
            "verified_transcripts": checked, "not_in_raw_wgbs": sorted(waiting),
            "result": "WGBS_COLUMNS_0_THROUGH_20_MATCH_SAVED_GRAPH"}


def sum_feature_values(grad, x, node: int, feature_map, torch):
    self_signed = grad[node].detach().cpu().numpy()
    self_abs = np.abs(self_signed)
    # Full receptive-field input sensitivity, not just the node itself.
    graph_abs_by_feat = grad.detach().abs().sum(dim=0).cpu().numpy()
    input_values = x[node].detach().cpu().numpy()
    local_total = float(self_abs.sum())
    graph_total = float(graph_abs_by_feat.sum())
    if not np.all(np.isfinite(self_signed)) or not np.all(np.isfinite(graph_abs_by_feat)):
        raise ValueError("Non-finite transcript attribution")
    features = []
    for m in feature_map:
        j = m["feature_index"]
        features.append({
            "feature_index": j, "assay": m["assay"], "region": m["region"],
            "raw_feature_column": m["raw_feature_column"],
            "self_input_value": float(input_values[j]),
            "self_input_missing_sentinel": bool(np.isclose(input_values[j], -1)),
            "self_gradient_signed": float(self_signed[j]),
            "self_gradient_abs": float(self_abs[j]),
            "self_gradient_abs_share": float(self_abs[j] / local_total) if local_total > 0 else np.nan,
            "all_nodes_gradient_abs_sum": float(graph_abs_by_feat[j]),
            "all_nodes_gradient_abs_share": float(graph_abs_by_feat[j] / graph_total) if graph_total > 0 else np.nan,
        })
    summaries = []
    for assay in ASSAYS:
        ids = [m["feature_index"] for m in feature_map if m["assay"] == assay]
        ii = next(m["feature_index"] for m in feature_map if m["assay"] == assay and m["region"] == "intnintrons")
        summaries.append({
            "assay": assay,
            "n_features": len(ids),
            "self_gradient_abs_sum": float(self_abs[ids].sum()),
            "self_gradient_abs_mean": float(self_abs[ids].mean()),
            "self_gradient_abs_share": float(self_abs[ids].sum() / local_total) if local_total > 0 else np.nan,
            "all_nodes_gradient_abs_sum": float(graph_abs_by_feat[ids].sum()),
            "all_nodes_gradient_abs_share": float(graph_abs_by_feat[ids].sum() / graph_total) if graph_total > 0 else np.nan,
            "internal_intron_self_gradient_signed": float(self_signed[ii]),
            "internal_intron_self_gradient_abs": float(self_abs[ii]),
            "local_abs_mass_all_features": local_total,
            "graph_abs_mass_all_features": graph_total,
            "self_mass_fraction_of_graph_mass": float(local_total / graph_total) if graph_total > 0 else np.nan,
        })
    return features, summaries


def atomic_csv(df: pd.DataFrame, path: Path):
    temp = path.with_name(path.name + f".tmp-{os.getpid()}")
    try:
        df.to_csv(temp, index=False)
        temp.replace(path)
    finally:
        if temp.exists():
            temp.unlink()


def perform(args):
    import torch
    import dgl
    from dgl.data import load_graphs

    if args.task_index < 0:
        raise ValueError("Negative task index")
    tasks = pd.read_csv(args.outdir / "attribution_tasks.csv", dtype={"eval_sample": str, "source_model": str})
    target_file = args.outdir / "attribution_targets.csv"
    targets = pd.read_csv(target_file, dtype={"eval_sample": str, "source_model": str, "txid": str})
    if args.task_index >= len(tasks) or int(tasks.iloc[args.task_index].task_index) != args.task_index:
        raise ValueError(f"Invalid task_index {args.task_index}: manifest has {len(tasks)} rows")
    task = tasks.iloc[args.task_index]
    sample, source = str(task.eval_sample), str(task.source_model)
    selected = targets.loc[(targets.eval_sample == sample) & (targets.source_model == source)].copy()
    if selected.empty:
        raise ValueError(f"No target rows for {source} -> {sample}")
    expected = selected.drop_duplicates("txid").reset_index(drop=True)
    if len(expected) != int(task.n_unique_transcripts):
        raise ValueError("Manifest unique transcript count mismatch")
    # Node- and graph-level files are in EX130, with no writes to historic EX18/EX51/EX88.
    details_dir = args.outdir / "attribution_by_sample"
    details_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = f"{source}__to__{sample}"
    attr_file = details_dir / (output_prefix + ".features.csv")
    summary_file = details_dir / (output_prefix + ".summary.csv")
    proof_file = details_dir / (output_prefix + ".provenance.json")
    model_dir = args.base / "EX51"
    graph_path = (args.base / "EX18" / "graphs_tpmpc1e-3_hmpc1e-9_zscale_clip5_WGBS+HM"
                  / "10sample_+thyroid" / f"{sample}.all.dgl")
    checkpoint_path = model_dir / "models" / f"{source}_128_neurons.model"
    if not graph_path.is_file() or not checkpoint_path.is_file():
        raise FileNotFoundError(f"Missing DGL graph or model: {graph_path}, {checkpoint_path}")
    if attr_file.exists() or summary_file.exists() or proof_file.exists():
        if args.skip_completed and all(p.is_file() and p.stat().st_size > 0 for p in (attr_file, summary_file, proof_file)):
            proof = json.loads(proof_file.read_text())
            current_checkpoint_sha = sha256_stream(checkpoint_path)
            if (proof.get("n_selected_transcripts") == len(expected)
                and proof.get("target_test_ids") == [int(x) for x in expected.test_id.tolist()]
                and proof.get("source_model") == source
                and proof.get("evaluated_sample") == sample
                and proof.get("checkpoint_sha256") == current_checkpoint_sha
                and proof.get("graph_size_bytes") == graph_path.stat().st_size
                and proof.get("graph_mtime_ns") == graph_path.stat().st_mtime_ns
                and proof.get("attribution_method") == "scalar_target_input_gradient"
                and proof.get("wgbs_first_feature_audit", {}).get("result") == "WGBS_COLUMNS_0_THROUGH_20_MATCH_SAVED_GRAPH"):
                print(f"[SKIP] exact checkpoint/graph/target identity matches prior verified result: {output_prefix}")
                return
        raise FileExistsError(f"Preexisting attribution output for {source} to {sample}; preserving instead of overwriting")

    sys.path.insert(0, str(model_dir))
    from gcn_regression import GCN
    torch.set_num_threads(max(1, args.cpu_threads))
    g = load_graphs(str(graph_path))[0][0].int()
    if args.self_loop:
        g = dgl.remove_self_loop(g)
        g = dgl.add_self_loop(g)
    degs = g.in_degrees().float()
    norm = degs.pow(-0.5)
    norm[torch.isinf(norm)] = 0
    g.ndata["norm"] = norm.unsqueeze(1)
    x0 = g.ndata["feat"]
    label = g.ndata["exp_label"]
    testid = g.ndata["test_id"]
    testmask = g.ndata["test_mask"].bool()
    if x0.ndim != 2 or x0.shape[1] != len(mapping()):
        raise ValueError(f"Unexpected graph feature width: {tuple(x0.shape)} (expected 147)")
    model = GCN(g, x0.shape[1], 128, 1, args.n_layers, torch.nn.functional.relu, args.dropout)
    try:
        weights = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
    except TypeError:
        # Old cluster PyTorch predates weights_only. Use ONLY known local trusted checkpoints.
        weights = torch.load(checkpoint_path, map_location="cpu")
    model.load_state_dict(weights, strict=True)
    model.eval()
    x0 = x0.detach()
    # Use actual test_id rather than assuming graph node index equals CSV test_id.
    matches = torch.nonzero(testmask, as_tuple=False).reshape(-1)
    testids = testid[matches].detach().cpu().numpy()
    id_to_node = {}
    for n, k in zip(matches.detach().cpu().tolist(), testids):
        ki = int(k)
        if ki in id_to_node:
            raise ValueError(f"Duplicate test_id={ki} in graph")
        id_to_node[ki] = int(n)
    txid_nodes = {}
    for row in expected.itertuples(index=False):
        k = int(row.test_id)
        if k not in id_to_node:
            raise ValueError(f"test_id {k} is absent from held-out mask in graph for {row.txid}")
        n = id_to_node[k]
        if not np.isclose(float(label[n]), float(row.expected_label), atol=args.label_tolerance, rtol=0):
            raise ValueError(f"CSV vs DGL label mismatch for {row.txid}: {float(label[n])} vs {row.expected_label}")
        txid_nodes[str(row.txid)] = n

    audit = {"result": "NOT_RUN_BY_EXPLICIT_OVERRIDE"} if args.skip_feature_audit else verify_wgbs_first(
        args.features_root, sample.split("_")[0], txid_nodes, x0)
    if (not args.skip_feature_audit) and audit["n_verified_transcripts"] < len(txid_nodes):
        print(f"[WARN] WGBS feature audit covered {audit['n_verified_transcripts']} of {len(txid_nodes)} targets", file=sys.stderr)

    with torch.no_grad():
        observed_pred = model(x0).reshape(-1).detach().cpu()
    for row in expected.itertuples(index=False):
        n = txid_nodes[str(row.txid)]
        actual = float(observed_pred[n])
        if abs(actual - float(row.expected_logit)) > args.logit_tolerance:
            raise ValueError(
                f"CHECKPOINT/GRAPH REPLAY FAILURE {row.txid} in {sample}: "
                f"EX51={row.expected_logit}, reloaded={actual}, tolerance={args.logit_tolerance}. "
                "Stop: likely model/graph/architecture mismatch."
            )

    x = x0.detach().clone().requires_grad_(True)
    fmap = mapping()
    features, summaries = [], []
    for i, row in enumerate(expected.itertuples(index=False), 1):
        n = txid_nodes[str(row.txid)]
        y, grad = target_node_gradient(model, x, n)
        detail, smry = sum_feature_values(grad, x, n, fmap, torch)
        common = {"source_model": source, "eval_sample": sample,
                  "txid": str(row.txid), "test_id": int(row.test_id), "chr": str(row.chr),
                  "observed_log2_tpm_pc1e3": float(row.expected_label),
                  "predicted_log2_tpm_pc1e3": float(row.expected_logit),
                  "predicted_replayed": float(y), "prediction_absolute_error":
                      abs(float(row.expected_logit) - float(row.expected_label))}
        features.extend([{**common, **z} for z in detail])
        summaries.extend([{**common, **z} for z in smry])
        if i % 5 == 0 or i == len(expected):
            print(f"[{sample}] attributed {i}/{len(expected)} distinct held-out transcript outputs", flush=True)
    atomic_csv(pd.DataFrame(features), attr_file)
    atomic_csv(pd.DataFrame(summaries), summary_file)
    proof = {
        "attribution_method": "scalar_target_input_gradient",
        "source_model": source, "evaluated_sample": sample,
        "n_selected_transcripts": len(expected), "target_test_ids": expected.test_id.tolist(),
        "feature_map": "EX18 WGBS FIRST; 21 features per mark, 147 total",
        "graph_num_nodes": int(g.num_nodes()), "graph_num_edges": int(g.num_edges()),
        "graph_path": str(graph_path), "checkpoint_path": str(checkpoint_path),
        "graph_size_bytes": graph_path.stat().st_size,
        "graph_mtime_ns": graph_path.stat().st_mtime_ns,
        "checkpoint_sha256": sha256_stream(checkpoint_path),
        "torch_version": torch.__version__, "dgl_version": dgl.__version__,
        "n_layers": args.n_layers, "n_hidden": 128, "dropout_in_eval": args.dropout,
        "self_loop": args.self_loop, "logit_tolerance": args.logit_tolerance,
        "wgbs_first_feature_audit": audit,
        "scope": "Target-specific input sensitivity, not biological causal attribution",
        "n_feature_rows": len(features), "n_assay_summary_rows": len(summaries),
    }
    proof_file.write_text(json.dumps(proof, indent=2))
    print(f"[DONE] {sample}: {len(expected)} transcript targets, checkpoint replay PASS, feature audit={audit['result']}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base", type=Path, default=Path.home() / "IHEC_Mar_2024")
    ap.add_argument("--outdir", type=Path, required=True,
                    help="Directory containing attribution_tasks.csv and attribution_targets.csv")
    ap.add_argument("--features-root", type=Path, default=Path("/home/shared/IHEC/baseline_expt/features"))
    ap.add_argument("--task-index", type=int, required=True)
    ap.add_argument("--n-layers", type=int, default=7)
    ap.add_argument("--dropout", type=float, default=0.45)
    ap.add_argument("--logit-tolerance", type=float, default=0.05)
    ap.add_argument("--label-tolerance", type=float, default=1e-4)
    ap.add_argument("--cpu-threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "4")))
    ap.add_argument("--self-loop", dest="self_loop", action="store_true", default=True)
    ap.add_argument("--no-self-loop", dest="self_loop", action="store_false")
    ap.add_argument("--skip-feature-audit", action="store_true", help="Only for diagnosis; results are not figure-ready")
    ap.add_argument("--skip-completed", action="store_true", help="Do not overwrite validated existing results")
    args = ap.parse_args()
    if args.cpu_threads < 1:
        ap.error("cpu-threads must be >=1")
    perform(args)


if __name__ == "__main__":
    main()
