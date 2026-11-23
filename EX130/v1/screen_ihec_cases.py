#!/usr/bin/env python3
"""Two-stage, read-only IHEC transcript case-study candidate screen.

Stage 1 uses all available DIAGONAL EX51 held-out prediction CSVs (chr18/19).
Stage 2 selectively joins original WGBS + H3K36me3 features for shortlisted
transcript/sample pairs. No existing input is modified and no model is retrained.

Requires: Python >=3.9, numpy, pandas.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

SAMPLE_RE = re.compile(r"^(IHECRE\d+)(?:_(mRNA-Seq|total-RNA-Seq))?$")
CHRS = ("chr18", "chr19")
H3_COLUMNS = [
    "hg38H3K36me3intnintrons", "hg38H3K36me3gene_body",
    "hg38H3K36me3firstintron", "hg38H3K36me3lastintron",
]
DNAM_COLUMNS = [
    "beta_hg38intnintrons", "beta_hg38gene_body",
    "beta_hg38firstintron", "beta_hg38promoter", "beta_hg38TSS1",
]


def metadata_frame(path: Path) -> pd.DataFrame:
    with path.open(encoding="utf-8-sig") as f:
        first = f.readline()
    # Note: the observed .tsv file is actually COMMA-delimited.
    sep = "," if first.count(",") > first.count("\t") else "\t"
    df = pd.read_csv(path, sep=sep, dtype=str).fillna("")
    required = ["Sample ID", "Tissue origin (simplified)", "Biomaterial type", "Disease status", "project"]
    for col in required:
        if col not in df:
            raise ValueError(f"Missing metadata column: {col!r}. Found {list(df.columns)}")
    mrna_cols = [x for x in df.columns if x.lower().startswith("mrna") and "availability" in x.lower()]
    total_cols = [x for x in df.columns if x.lower().startswith("total") and "availability" in x.lower()]
    if len(mrna_cols) != 1 or len(total_cols) != 1:
        raise ValueError("Cannot identify mRNA and total-RNA availability columns")
    if df["Sample ID"].duplicated().any():
        raise ValueError("Duplicate Sample ID in metadata")
    df = df.rename(columns={mrna_cols[0]: "mrna_available", total_cols[0]: "total_available"})
    return df.set_index("Sample ID", drop=False)


def build_sample_manifest(args, md: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    if not args.samples_dir.is_dir():
        raise FileNotFoundError(args.samples_dir)
    names = sorted(x.name for x in args.samples_dir.iterdir() if x.is_dir() and SAMPLE_RE.fullmatch(x.name))
    rows, errors = [], []
    for sid in names:
        match = SAMPLE_RE.fullmatch(sid)
        base, suffix = match.groups()
        if base not in md.index:
            errors.append(f"No metadata for {sid}")
            continue
        meta = md.loc[base]
        if suffix == "mRNA-Seq":
            rna = "mRNA"
        elif suffix == "total-RNA-Seq":
            rna = "total"
        else:
            mrna, total = str(meta.mrna_available).upper() == "T", str(meta.total_available).upper() == "T"
            if mrna == total:
                errors.append(f"Ambiguous RNA protocol for unsuffixed sample {sid}")
                continue
            rna = "mRNA" if mrna else "total"
        rows.append({
            "sample_id": sid, "base_id": base, "rna_type": rna,
            "tissue": str(meta["Tissue origin (simplified)"]).strip(),
            "biomaterial": str(meta["Biomaterial type"]).strip(),
            "disease": str(meta["Disease status"]).strip(),
            "project": str(meta["project"]).strip(),
        })
    sample_df = pd.DataFrame(rows)
    if sample_df.empty:
        raise ValueError("No valid sample entries in manifest")
    if sample_df["sample_id"].duplicated().any():
        raise ValueError("Duplicate sample IDs")
    return sample_df, errors


def observed_tpm(y):
    return np.maximum(0, np.exp2(y) - 0.001)


def screen(args):
    args.outdir.mkdir(parents=True, exist_ok=True)
    md = metadata_frame(args.metadata)
    samples, issues = build_sample_manifest(args, md)
    if args.trained_sample and args.trained_sample not in set(samples.sample_id):
        raise ValueError(f"Fixed training model not in sample manifest: {args.trained_sample}")
    if args.max_samples:
        samples = samples.head(args.max_samples)
    samples.to_csv(args.outdir / "sample_manifest.csv", index=False)
    frames, qc = [], []
    for k, s in enumerate(samples.to_dict("records"), 1):
        sid = s["sample_id"]
        model_sample = args.trained_sample or sid
        path = args.logits_dir / f"10sample_+thyroid_{model_sample}_{sid}_128_test_logits.csv"
        if not path.is_file():
            issues.append(f"Missing diagonal test file: {path}")
            continue
        df = pd.read_csv(
            path, usecols=["test_id", "logits", "labels", "txid", "chr"],
            dtype={"test_id": "Int64", "logits": "float32", "labels": "float32", "txid": str, "chr": str},
        )
        df = df.loc[df.chr.isin(CHRS)].copy()
        df = df.loc[np.isfinite(df.logits) & np.isfinite(df.labels) & df.txid.notna()].copy()
        if df.txid.duplicated().any():
            raise ValueError(f"Transcript IDs duplicated within {path}")
        if not len(df):
            issues.append(f"No usable test rows in {path}")
            continue
        qc.append({"sample_id": sid, "n_test_transcripts": len(df),
                   "mae_log2": float(np.mean(np.abs(df.logits - df.labels))),
                   "zero_like_fraction": float(np.mean(df.labels <= np.log2(0.0011)))})
        df = df.rename(columns={"logits": "pred_log2", "labels": "obs_log2"})
        df["sample_id"] = sid
        for attr in ("base_id", "rna_type", "tissue", "biomaterial", "disease", "project"):
            df[attr] = s[attr]
        df["abs_err"] = (df.pred_log2 - df.obs_log2).abs()
        frames.append(df.drop(columns=["test_id", "chr"]))
        if k % 30 == 0:
            print(f"Read {k}/{len(samples)} diagonal CSVs", file=sys.stderr, flush=True)
    if not frames:
        raise ValueError("No held-out prediction rows loaded")
    pd.DataFrame(qc).to_csv(args.outdir / "sample_qc.csv", index=False)
    obs = pd.concat(frames, ignore_index=True)
    del frames
    for col in ("txid", "sample_id", "base_id", "rna_type", "tissue", "biomaterial", "disease", "project"):
        obs[col] = obs[col].astype("category")
    if not args.include_disease:
        obs = obs.loc[obs.disease == "non-disease"].copy()
    if obs.empty:
        raise ValueError("No eligible rows after metadata filters")
    keys = ["txid", "rna_type", "tissue", "biomaterial"]
    grp = obs.groupby(keys, observed=True, sort=False).agg(
        n=("sample_id", "nunique"), obs_median=("obs_log2", "median"),
        pred_median=("pred_log2", "median"), mae=("abs_err", "mean")
    ).reset_index()
    grp = grp.loc[(grp.n >= args.min_group_n) & (grp.mae <= args.max_group_mae)].copy()
    grp["median_tpm"] = observed_tpm(grp.obs_median)
    grp.to_csv(args.outdir / "tissue_group_summary.csv", index=False)

    highs = grp.loc[grp.median_tpm >= args.high_tpm].copy()
    lows = grp.loc[grp.median_tpm <= args.low_tpm].copy()
    # A cross-join within a transcript and RNA protocol, never across protocols.
    pairs = highs.merge(lows, on=["txid", "rna_type"], suffixes=("_high", "_low"))
    pairs = pairs.loc[pairs.tissue_high != pairs.tissue_low].copy()
    if not args.allow_mixed_biomaterial:
        pairs = pairs.loc[pairs.biomaterial_high == pairs.biomaterial_low].copy()
    pairs["observed_delta_log2"] = pairs.obs_median_high - pairs.obs_median_low
    pairs["predicted_delta_log2"] = pairs.pred_median_high - pairs.pred_median_low
    pairs["delta_discrepancy_log2"] = (pairs.predicted_delta_log2 - pairs.observed_delta_log2).abs()
    pairs = pairs.loc[
        (pairs.observed_delta_log2 >= args.min_observed_delta) &
        (pairs.predicted_delta_log2 >= args.min_predicted_delta) &
        (pairs.delta_discrepancy_log2 <= args.max_delta_discrepancy)
    ].copy()
    if pairs.empty:
        issues.append("No contrast pairs passed filters; relax thresholds with CLI flags")
        (args.outdir / "issues.json").write_text(json.dumps(issues, indent=2))
        print("No candidate pairs passed. Wrote sample QC and tissue groups.")
        return
    pairs["case_type"] = np.where(pairs.median_tpm_low >= args.graded_low_tpm, "graded", "on_off")
    # Screening heuristic, not a biological statistic or p-value.
    pairs["prediction_screen_score"] = (
        pairs.observed_delta_log2.clip(upper=14) -
        1.2 * pairs.delta_discrepancy_log2 -
        0.5 * (pairs.mae_high + pairs.mae_low) -
        0.75 * (pairs.case_type == "on_off").astype(int)
    )
    pairs = pairs.sort_values("prediction_screen_score", ascending=False).reset_index(drop=True)
    pairs.insert(0, "candidate_id", [f"C{i:05d}" for i in range(1, len(pairs) + 1)])
    pairs.to_csv(args.outdir / "candidate_pairs.csv", index=False)
    # Unique transcript/protocol representative for screening. Preserve both classes.
    unique = pairs.drop_duplicates(["txid", "rna_type", "case_type"])
    unique.to_csv(args.outdir / "candidate_shortlist.csv", index=False)

    selected = (pd.concat([
        unique.loc[unique.case_type == "graded"].head(args.top_per_type),
        unique.loc[unique.case_type == "on_off"].head(args.top_per_type),
    ], ignore_index=True))
    chosen_txs = set(selected.txid.astype(str))
    small = obs.loc[obs.txid.astype(str).isin(chosen_txs)].copy()
    by_tx = {str(t): d for t, d in small.groupby("txid", observed=True)}
    reps = []
    for p in selected.to_dict("records"):
        pool = by_tx.get(str(p["txid"]))
        if pool is None:
            continue
        for role in ("high", "low"):
            tissue = p[f"tissue_{role}"]
            bm = p[f"biomaterial_{role}"]
            rs = pool.loc[(pool.rna_type == p["rna_type"]) &
                          (pool.tissue == tissue) & (pool.biomaterial == bm)].copy()
            rs["pick_distance"] = ((rs.obs_log2 - p[f"obs_median_{role}"]).abs() +
                                    (rs.pred_log2 - p[f"pred_median_{role}"]).abs() + rs.abs_err)
            rs = rs.sort_values("pick_distance").drop_duplicates("base_id").head(args.rep_samples)
            for _, z in rs.iterrows():
                reps.append({"candidate_id": p["candidate_id"], "txid": str(p["txid"]),
                             "case_type": p["case_type"], "role": role,
                             "sample_id": str(z.sample_id), "base_id": str(z.base_id),
                             "rna_type": str(z.rna_type), "tissue": str(z.tissue),
                             "biomaterial": str(z.biomaterial), "disease": str(z.disease),
                             "project": str(z.project), "observed_log2": float(z.obs_log2),
                             "predicted_log2": float(z.pred_log2),
                             "observed_tpm": float(observed_tpm(z.obs_log2)),
                             "prediction_abs_error_log2": float(z.abs_err)})
    pd.DataFrame(reps).to_csv(args.outdir / "representative_samples.csv", index=False)
    summary = {
        "sample_manifest_count": int(len(samples)), "loaded_sample_count": int(len(qc)),
        "unique_biospecimens_in_manifest": int(samples.base_id.nunique()),
        "test_transcript_ids": int(obs.txid.nunique()),
        "eligible_observation_rows": int(len(obs)), "group_rows_after_filters": int(len(grp)),
        "candidate_pairs": int(len(pairs)), "distinct_candidate_ids": int(len(unique)),
        "representative_records": int(len(reps)), "issues": issues,
        "analysis_scope": "Held-out chr18/chr19 transcript predictions; descriptive cross-context expression contrasts",
        "model_scope": "fixed_cross_sample" if args.trained_sample else "diagonal_sample_specific_models",
        "trained_sample": args.trained_sample,
        "causal_inference": False, "transcript_attributions_computed": False,
        "feature_index_source": "EX18 graph tensor order, NOT EX88 RF reference label order",
    }
    (args.outdir / "run_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


def extract_feature_rows(path: Path, keep_ids: set[str], feature_cols: list[str], prefix: str) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    header = list(pd.read_csv(path, sep="\t", nrows=0).columns)
    keep = ["transcript.id", "gene.id", "symbol"] + feature_cols
    if "transcript.id" not in header:
        raise ValueError(f"Missing transcript.id in {path}")
    actual = [x for x in keep if x in header]
    if len(actual) < len(keep):
        print(f"WARN: missing columns in {path.name}: {sorted(set(keep) - set(actual))}", file=sys.stderr)
    hits = []
    for chunk in pd.read_csv(path, sep="\t", usecols=actual,
                             chunksize=60000, dtype={"transcript.id": str}):
        sub = chunk.loc[chunk["transcript.id"].isin(keep_ids)]
        if len(sub):
            hits.append(sub.copy())
    if not hits:
        return pd.DataFrame(columns=["txid", "gene.id", "symbol"] + [prefix + col for col in feature_cols])
    df = pd.concat(hits, ignore_index=True).rename(columns={"transcript.id": "txid"})
    if df.txid.duplicated().any():
        raise ValueError(f"Duplicate transcript.id in {path}")
    return df.rename(columns={x: prefix + x for x in feature_cols if x in df.columns})


def features(args):
    if not args.outdir.is_dir():
        raise FileNotFoundError(f"First run screen to create {args.outdir}")
    pair_path = args.outdir / "candidate_pairs.csv"
    reps_path = args.outdir / "representative_samples.csv"
    pairs = pd.read_csv(pair_path)
    reps = pd.read_csv(reps_path)
    if len(reps) == 0:
        raise ValueError("No representative rows from stage 1")
    by_base = defaultdict(set)
    for _, row in reps.iterrows():
        by_base[str(row.base_id)].add(str(row.txid))
    all_features, issues = [], []
    for i, (base, txids) in enumerate(sorted(by_base.items()), 1):
        directory = args.features_root / base
        hpaths = sorted(directory.glob(f"{base}.[0-9]*.H3K36me3.pval.all.txt"))
        mpaths = sorted(directory.glob(f"{base}.[0-9]*.standard-WGBS.cpg.all.txt"))
        # A second naming convention can exist in archived datasets; avoid silently using it.
        if len(hpaths) != 1 or len(mpaths) != 1:
            issues.append({"base_id": base, "H3K36me3_files": len(hpaths), "WGBS_files": len(mpaths)})
            continue
        h3 = extract_feature_rows(hpaths[0], txids, H3_COLUMNS, "h3_")
        met = extract_feature_rows(mpaths[0], txids, DNAM_COLUMNS, "dna_")
        if "gene.id" in met:
            met = met.drop(columns=["gene.id", "symbol"], errors="ignore")
        frame = h3.merge(met, on="txid", how="outer", validate="one_to_one")
        frame["base_id"] = base
        all_features.append(frame)
        if i % 20 == 0:
            print(f"Feature-read {i}/{len(by_base)} biospecimens", file=sys.stderr, flush=True)
    if not all_features:
        raise ValueError(f"No usable feature files. Example errors: {issues[:3]}")
    combined = pd.concat(all_features, ignore_index=True)
    for col in H3_COLUMNS:
        name = "h3_" + col
        if name not in combined:
            combined[name] = np.nan
    for col in DNAM_COLUMNS:
        name = "dna_" + col
        if name not in combined:
            combined[name] = np.nan
    for col in ("h3_" + x for x in H3_COLUMNS):
        combined[col] = pd.to_numeric(combined[col], errors="coerce")
        combined.loc[combined[col] < 0, col] = np.nan
    for col in ("dna_" + x for x in DNAM_COLUMNS):
        combined[col] = pd.to_numeric(combined[col], errors="coerce")
        combined.loc[~combined[col].between(0, 1), col] = np.nan
    annotated = reps.merge(combined, on=["base_id", "txid"], how="left", validate="many_to_one")
    annotated.to_csv(args.outdir / "representative_feature_values.csv", index=False)

    h3_col = "h3_hg38H3K36me3intnintrons"
    dna_col = "dna_beta_hg38intnintrons"
    checks = []
    for candidate, group in annotated.groupby("candidate_id", sort=False):
        p = pairs.loc[pairs.candidate_id == candidate]
        if p.empty:
            continue
        p = p.iloc[0]
        hi = group[group.role == "high"]
        lo = group[group.role == "low"]
        def safe_med(frame, col):
            v = pd.to_numeric(frame[col], errors="coerce").dropna()
            return float(v.median()) if len(v) else np.nan
        h_hi, h_lo = safe_med(hi, h3_col), safe_med(lo, h3_col)
        d_hi, d_lo = safe_med(hi, dna_col), safe_med(lo, dna_col)
        hdelta = float(np.log2(1 + h_hi) - np.log2(1 + h_lo)) if np.isfinite(h_hi) and np.isfinite(h_lo) else np.nan
        ddelta = float(d_hi - d_lo) if np.isfinite(d_hi) and np.isfinite(d_lo) else np.nan
        h_support = np.isfinite(hdelta) and hdelta >= args.h3_delta_threshold
        d_support = np.isfinite(ddelta) and abs(ddelta) >= args.dna_delta_threshold
        cls = "both" if h_support and d_support else "H3K36me3" if h_support else "DNAme" if d_support else "neither_or_missing"
        gene = next((str(x) for x in group.get("symbol", pd.Series(dtype=str)).dropna() if str(x)), "")
        checks.append({"candidate_id": candidate, "txid": str(p.txid), "gene_symbol": gene,
                       "rna_type": p.rna_type, "tissue_high": p.tissue_high, "tissue_low": p.tissue_low,
                       "obs_high_tpm": p.median_tpm_high, "obs_low_tpm": p.median_tpm_low,
                       "observed_delta_log2": p.observed_delta_log2,
                       "predicted_delta_log2": p.predicted_delta_log2,
                       "prediction_screen_score": p.prediction_screen_score,
                       "h3_high_internal_introns": h_hi, "h3_low_internal_introns": h_lo,
                       "h3_delta_log2_1plus": hdelta,
                       "dna_high_internal_intron_beta": d_hi, "dna_low_internal_intron_beta": d_lo,
                       "dna_beta_difference": ddelta, "n_high_reps": len(hi), "n_low_reps": len(lo),
                       "feature_support_class": cls, "case_type": p.case_type})
    table = pd.DataFrame(checks)
    if table.empty:
        raise ValueError("No candidate feature checks generated")
    class_rank = {"both": 0, "H3K36me3": 1, "DNAme": 2, "neither_or_missing": 3}
    table["review_priority"] = table.feature_support_class.map(class_rank)
    table = table.sort_values(["review_priority", "prediction_screen_score"], ascending=[True, False])
    table.to_csv(args.outdir / "candidate_ranked_with_features.csv", index=False)
    (args.outdir / "feature_issues.json").write_text(json.dumps(issues, indent=2))
    print(f"Annotated {len(table)} candidates, {len(issues)} missing/ambiguous feature-file pairs")
    print("These are cross-context associations, NOT individual-transcript model attribution or causality.")
    print(table.head(15).to_string(index=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("stage", choices=("screen", "features"))
    parser.add_argument("--base", type=Path, default=Path.home() / "IHEC_Mar_2024")
    parser.add_argument("--logits-dir", type=Path)
    parser.add_argument("--samples-dir", type=Path)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--features-root", type=Path, default=Path("/home/shared/IHEC/baseline_expt/features"))
    parser.add_argument("--outdir", type=Path, default=Path.home() / "IHEC_case_screen")
    parser.add_argument("--trained-sample", type=str, default=None,
                        help="Optional single fixed source model used for predictions in every evaluated sample; default is diagonal sample-specific models")
    parser.add_argument("--include-disease", action="store_true")
    parser.add_argument("--allow-mixed-biomaterial", action="store_true")
    parser.add_argument("--min-group-n", type=int, default=3)
    parser.add_argument("--high-tpm", type=float, default=5)
    parser.add_argument("--low-tpm", type=float, default=1)
    parser.add_argument("--graded-low-tpm", type=float, default=0.1)
    parser.add_argument("--max-group-mae", type=float, default=2)
    parser.add_argument("--min-observed-delta", type=float, default=3)
    parser.add_argument("--min-predicted-delta", type=float, default=1.5)
    parser.add_argument("--max-delta-discrepancy", type=float, default=3)
    parser.add_argument("--top-per-type", type=int, default=30)
    parser.add_argument("--rep-samples", type=int, default=2)
    parser.add_argument("--h3-delta-threshold", type=float, default=0.7)
    parser.add_argument("--dna-delta-threshold", type=float, default=0.2)
    parser.add_argument("--max-samples", type=int, default=0,
                        help="Optional sorted manifest prefix for smoke-test only; 0 reads ALL 358")
    args = parser.parse_args()
    args.logits_dir = args.logits_dir or args.base / "EX51" / "logits_debug_pc1"
    args.samples_dir = args.samples_dir or args.base / "EX18" / "samples_needed"
    args.metadata = args.metadata or args.base / "metadata" / "SuppTable1_with_project.tsv"
    if args.stage == "screen":
        screen(args)
    else:
        features(args)


if __name__ == "__main__":
    main()
