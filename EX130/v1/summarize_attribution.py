#!/usr/bin/env python3
"""Join transcript-targeted GCN gradients to candidate contrasts for manual review.

A strict merger: fails if any scheduled attribution file or provenance receipt is
missing, or if EX18 WGBS feature-order audit was not verified. Output is
exploratory and should not be described as causal or directionally regulatory.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def require(df, columns, name):
    absent = set(columns) - set(df.columns)
    if absent:
        raise ValueError(f"{name}: missing {sorted(absent)}")


def summarize(args):
    d = args.outdir
    task_df = pd.read_csv(d / "attribution_tasks.csv", dtype={"source_model": str, "eval_sample": str})
    target_df = pd.read_csv(d / "attribution_targets.csv", dtype={"source_model": str, "eval_sample": str, "txid": str, "candidate_id": str})
    require(target_df, ["source_model", "eval_sample", "txid", "candidate_id", "role", "tissue"], "targets")
    all_summaries, failures = [], []
    for task in task_df.itertuples(index=False):
        prefix = f"{task.source_model}__to__{task.eval_sample}"
        folder = d / "attribution_by_sample"
        f = folder / (prefix + ".summary.csv")
        p = folder / (prefix + ".provenance.json")
        if not f.is_file() or not p.is_file():
            failures.append(prefix + " missing summary or provenance")
            continue
        proof = json.loads(p.read_text())
        audit = proof.get("wgbs_first_feature_audit", {})
        if not args.allow_unaudited and audit.get("result") != "WGBS_COLUMNS_0_THROUGH_20_MATCH_SAVED_GRAPH":
            failures.append(prefix + " missing WGBS-first graph verification")
            continue
        if proof.get("source_model") != task.source_model or proof.get("evaluated_sample") != task.eval_sample:
            failures.append(prefix + " source/target provenance mismatch")
            continue
        frame = pd.read_csv(f, dtype={"txid": str, "source_model": str, "eval_sample": str})
        if frame.txid.nunique() != int(task.n_unique_transcripts) or len(frame) != int(task.n_unique_transcripts) * 7:
            failures.append(prefix + " wrong attribution row count")
            continue
        all_summaries.append(frame)
    (d / "attribution_merge_issues.json").write_text(json.dumps(failures, indent=2))
    if failures:
        raise RuntimeError(f"Incomplete or unverified attribution batch: {len(failures)} tasks. Examples: {failures[:4]}")
    data = pd.concat(all_summaries, ignore_index=True)
    require(data, ["self_gradient_abs_share", "all_nodes_gradient_abs_share", "assay"], "source attribution")
    for (s, e, t), frame in data.groupby(["source_model", "eval_sample", "txid"]):
        if set(frame.assay) != {"DNAm", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3"}:
            raise ValueError(f"Missing assay(s): {s}, {e}, {t}")
        if float(frame.self_gradient_abs_share.sum()) > 0:
            if not np.isclose(frame.self_gradient_abs_share.sum(), 1.0, atol=0.001):
                raise ValueError(f"Self-gradients shares do not sum to 1: {s},{e},{t}")
    joined = target_df.merge(data, on=["source_model", "eval_sample", "txid"], how="left", validate="many_to_many")
    if joined.assay.isna().any():
        raise ValueError("Missing attribution data for requested targets")
    joined.to_csv(d / "candidate_per_sample_assay_attribution.csv", index=False)
    summaries = joined.groupby(["candidate_id", "origin", "txid", "role", "tissue", "assay"], dropna=False).agg(
        n_samples=("eval_sample", "nunique"),
        mean_self_gradient_share=("self_gradient_abs_share", "mean"),
        median_self_gradient_share=("self_gradient_abs_share", "median"),
        mean_graph_gradient_share=("all_nodes_gradient_abs_share", "mean"),
        mean_internal_intron_self_abs=("internal_intron_self_gradient_abs", "mean"),
        mean_internal_intron_self_signed=("internal_intron_self_gradient_signed", "mean"),
        mean_abs_expression_prediction_error=("prediction_absolute_error", "mean"),
    ).reset_index()
    summaries.to_csv(d / "candidate_role_assay_summary.csv", index=False)
    # Pivot into one row per candidate+role, with the two marks beside each other.
    sub = summaries.loc[summaries.assay.isin(["DNAm", "H3K36me3"])].copy()
    pivot = sub.pivot_table(
        index=["candidate_id", "origin", "txid", "role", "tissue"],
        columns="assay", values=["mean_self_gradient_share", "mean_graph_gradient_share",
                                "mean_internal_intron_self_abs", "n_samples"],
        aggfunc="first",
    )
    pivot.columns = [f"{metric}_{assay}" for metric, assay in pivot.columns]
    pivot = pivot.reset_index()
    pivot["combined_self_share_H3K36me3_plus_DNAm"] = (
        pivot["mean_self_gradient_share_H3K36me3"] + pivot["mean_self_gradient_share_DNAm"]
    )
    pivot.to_csv(d / "candidate_review_attribution.csv", index=False)
    info = {
        "n_scheduled_tasks": len(task_df), "n_target_rows": len(target_df),
        "n_unique_transcripts": target_df.txid.nunique(),
        "n_cases": target_df.candidate_id.nunique(),
        "selection_scope": "chr18/chr19 held-out EX51; one fixed source model for compared contexts",
        "attribution": "target-specific scalar output input gradients, both self-row and all-nodes",
        "interpretation": "predictive sensitivity only; neither causality nor directionality of biology",
        "pass": True,
    }
    (d / "attribution_analysis_summary.json").write_text(json.dumps(info, indent=2))
    print(json.dumps(info, indent=2))
    print(pivot.sort_values("combined_self_share_H3K36me3_plus_DNAm", ascending=False).head(15).to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--allow-unaudited", action="store_true", help="Diagnostics only, not publication-ready")
    summarize(ap.parse_args())


if __name__ == "__main__":
    main()
