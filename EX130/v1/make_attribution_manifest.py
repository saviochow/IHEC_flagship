#!/usr/bin/env python3
"""Freeze a bounded manifest of (source model, evaluated sample, transcript) targets.

Auto candidates are taken from screen_ihec_cases.py's representative_samples.csv.
Known pilot cases are added from seed_cases.csv using real held-out EX51 logits.
Only metadata-matched contexts and exact transcript IDs are admitted.
Seed cases may use different explicitly named fixed source models; each CASE is
still evaluated with a single frozen GCN across its compared biological contexts.
No graph/model loads or input mutations occur here.
"""
from __future__ import annotations

import argparse
import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED = {"test_id", "logits", "labels", "txid", "chr"}


def require_cols(df, required, filename):
    missing = set(required) - set(df)
    if missing:
        raise ValueError(f"{filename}: missing {sorted(missing)}")


def build_manifest(args):
    args.outdir.mkdir(parents=True, exist_ok=True)
    screen = args.screen_dir
    rs_path = screen / "run_summary.json"
    run_summary = json.loads(rs_path.read_text())
    if run_summary.get("model_scope") != "fixed_cross_sample":
        raise ValueError("Attribution case contrasts require a FIXED source model; rerun screen with --trained-sample")
    source = run_summary["trained_sample"]
    if not source:
        raise ValueError("Missing fixed model ID")
    sample_manifest = pd.read_csv(screen / "sample_manifest.csv", dtype=str).fillna("")
    require_cols(sample_manifest, ("sample_id", "base_id", "rna_type", "tissue", "disease"), "sample_manifest")
    if sample_manifest.sample_id.duplicated().any():
        raise ValueError("Duplicate evaluated sample IDs")
    allowed = set(sample_manifest.sample_id)
    meta = sample_manifest.set_index("sample_id")

    @lru_cache(maxsize=1000)
    def get_logits(model_id: str, eval_id: str):
        if eval_id not in allowed or model_id not in allowed:
            raise ValueError(f"Out-of-manifest source or sample: {model_id} -> {eval_id}")
        path = args.logits_dir / f"10sample_+thyroid_{model_id}_{eval_id}_128_test_logits.csv"
        if not path.is_file():
            return None
        df = pd.read_csv(path, dtype={"txid": str}, usecols=list(REQUIRED))
        df = df.loc[df["chr"].isin(["chr18", "chr19"]) & df.txid.notna()].copy()
        if df.txid.duplicated().any():
            raise ValueError(f"Duplicate held-out transcripts in {path}")
        return df.set_index("txid", drop=False)

    auto = pd.read_csv(screen / "representative_samples.csv", dtype=str).fillna("")
    require_cols(auto, ("candidate_id", "txid", "sample_id", "role", "rna_type", "tissue"), "representative_samples")
    candidate_rank = screen / "candidate_ranked_with_features.csv"
    if candidate_rank.is_file():
        ranks = pd.read_csv(candidate_rank, dtype={"candidate_id": str}).head(args.max_auto_candidates)
        allowed_candidates = set(ranks.candidate_id)
    else:
        ranked_pairs = pd.read_csv(screen / "candidate_pairs.csv", dtype={"candidate_id": str}).head(args.max_auto_candidates)
        allowed_candidates = set(ranked_pairs.candidate_id)
    auto = auto.loc[auto.candidate_id.isin(allowed_candidates)]
    requested, warnings = [], []

    def add(case_id, txid, eval_id, role, origin, model_id):
        if eval_id not in allowed:
            warnings.append(f"unrecognized sample {eval_id} for {case_id}")
            return
        df = get_logits(model_id, eval_id)
        if df is None or txid not in df.index:
            warnings.append(f"missing held-out transcript {txid} in {eval_id} / {case_id}")
            return
        x = df.loc[txid]
        pred, obs = float(x["logits"]), float(x["labels"])
        if not np.isfinite(pred) or not np.isfinite(obs):
            warnings.append(f"nonfinite logits/label {txid} in {eval_id}")
            return
        sm = meta.loc[eval_id]
        requested.append({
            "candidate_id": case_id, "origin": origin, "source_model": model_id,
            "eval_sample": eval_id, "txid": txid, "test_id": int(x["test_id"]),
            "chr": str(x["chr"]), "role": role, "rna_type": sm.rna_type,
            "tissue": sm.tissue, "base_id": sm.base_id,
            "expected_logit": pred, "expected_label": obs,
            "absolute_prediction_error": abs(pred - obs),
        })

    for r in auto.itertuples(index=False):
        add(str(r.candidate_id), str(r.txid), str(r.sample_id), str(r.role), "automated", source)

    if args.seed_cases.is_file():
        cases = pd.read_csv(args.seed_cases, dtype=str).fillna("")
        require_cols(cases, ("case_id", "txid", "rna_type", "tissue_high", "tissue_low", "source_model"), "seed_cases")
        for r in cases.itertuples(index=False):
            if r.source_model not in allowed:
                warnings.append(f"skipped seed {r.case_id}: trained source model {r.source_model} absent from sample manifest")
                continue
            for role, tissue in (("high", r.tissue_high), ("low", r.tissue_low)):
                options = sample_manifest.loc[
                    (sample_manifest.rna_type == r.rna_type) &
                    (sample_manifest.tissue == tissue) &
                    (sample_manifest.disease == "non-disease")
                ].sample_id.tolist()
                eligible = []
                for sid in options:
                    df = get_logits(r.source_model, sid)
                    if df is None or r.txid not in df.index:
                        continue
                    z = df.loc[r.txid]
                    if np.isfinite(z["logits"]) and np.isfinite(z["labels"]):
                        eligible.append((sid, float(z["labels"]), float(z["logits"])))
                if len(eligible) < args.min_seed_group:
                    warnings.append(f"seed {r.case_id} role={role} has only {len(eligible)} available datasets")
                    continue
                med_obs = float(np.median([v[1] for v in eligible]))
                med_pred = float(np.median([v[2] for v in eligible]))
                eligible.sort(key=lambda z: (abs(z[1]-med_obs) + abs(z[2]-med_pred) + abs(z[1]-z[2]), z[0]))
                for sid, *_ in eligible[:args.seed_reps]:
                    add(r.case_id, r.txid, sid, role, "pilot_seed", r.source_model)

    if not requested:
        raise ValueError(f"No valid target rows: {warnings[:10]}")
    targets = pd.DataFrame(requested).drop_duplicates(["candidate_id", "source_model", "eval_sample", "txid", "role"])
    targets = targets.sort_values(["source_model", "eval_sample", "txid", "candidate_id"]).reset_index(drop=True)
    targets.to_csv(args.outdir / "attribution_targets.csv", index=False)
    tasks = (targets.groupby(["source_model", "eval_sample"], sort=True)
             .agg(n_target_rows=("txid", "size"), n_unique_transcripts=("txid", "nunique"))
             .reset_index())
    tasks.insert(0, "task_index", list(range(len(tasks))))
    tasks.to_csv(args.outdir / "attribution_tasks.csv", index=False)
    manifest = {
        "primary_screen_source_model": source,
        "included_source_models": sorted(set(tasks.source_model)),
        "outdir": str(args.outdir), "tasks": len(tasks),
        "candidate_ids": int(targets.candidate_id.nunique()),
        "unique_sample_transcript_pairs": int(targets[["eval_sample", "txid"]].drop_duplicates().shape[0]),
        "target_rows": len(targets), "warnings": warnings,
        "feature_order": "EX18: DNAm (0-20), six HM blocks (21-146); verify against DGL tensor",
        "selection_note": "exploratory held-out sample selection; does not imply causal feature importance",
    }
    (args.outdir / "attribution_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps(manifest, indent=2))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--screen-dir", type=Path, required=True)
    ap.add_argument("--logits-dir", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--seed-cases", type=Path, default=Path(__file__).with_name("seed_cases.csv"))
    ap.add_argument("--seed-reps", type=int, default=2)
    ap.add_argument("--min-seed-group", type=int, default=3)
    ap.add_argument("--max-auto-candidates", type=int, default=20)
    args = ap.parse_args()
    if args.seed_reps < 1 or args.max_auto_candidates < 0:
        ap.error("Invalid positive sample/candidate limits")
    build_manifest(args)


if __name__ == "__main__":
    main()
