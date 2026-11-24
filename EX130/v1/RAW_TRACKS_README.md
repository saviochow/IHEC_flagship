# EX130 revised pilot provenance and parallel genomic tracks (v3)

## Status and raw-file verification

ClusterMCP confirmed that the source data for **IHECRE00001892** are under `~/data_store/IHEC_data_processing_Aug_2023_partial_epigenomes/IHECRE00001892/`, with `ChIP-Seq/H3K36me3.pval.signal.bigwig` and indexed `WGBS/*.standard-WGBS.cpg.all.bed.gz`. The main raw directory holds only RNA-seq for that sample. The track reader now searches **both** raw trees without mixing files between source directories.

ClusterMCP also confirmed the raw H3K36me3 BigWig and raw whole-region WGBS BED for **IHECRE00003822** in the main raw directory. Its whole-region `*.cpg.all.bed.gz` lacks a `.tbi` index in the inspected directory, while the exon-specific `*.cpg.4exon.bed.gz.tbi` exists. This is an **index availability issue, not a missing biological sample or raw WGBS data**. Sample 03822 is kept in all sample groups, expression screens and H3K36me3 panels. Whole-region CpG extraction requires indexing or an explicitly limited alternative; the plotter does not pretend the exon-specific file covers all introns.

In the revised 35-entry track sample manifest, the original 33 curated samples retain their original identifiers and both previous alternatives remain optional. The v3 manifest adds back `IHECRE00001892_mRNA-Seq` and `IHECRE00003822` and prefers the original three evaluated sample IDs for their respective pilot groups.

File-level remote checks were performed through ClusterMCP, **not** actual BigWig or Tabix interval decoding. A live Yipcluster compute job is still required to verify each transcript region. No Slurm submission was performed through ClusterMCP because its Yipcluster gateway was disabled.

## Exact sample IDs supporting the 17 pilot comparison rows

`pilot_candidates.csv` and `pilot_candidates_with_samples.csv` now contain the following columns in addition to all original prediction metrics:

- `high_sample_ids`, `low_sample_ids`: all **three** sample IDs included in each original *pilot* median and MAE, separated by semicolons. These are not claims that a tissue has only three biospecimens across the full cohort.
- `n_high_samples`, `n_low_samples`: both 3 for the pilot, preserving exactly the calculations previously described.
- `high_trained_sample_ids`, `low_trained_sample_ids`: for `diagonal` cases, the model for each sample is trained on itself; for `fixed_cross_sample` rows both roles use the single checkpoint in `trained_sample`.
- `sample_scope`: provenance clarification to prevent accidental comparison with the full 358-dataset screen.

The `trained_sample` cell in diagonal rows now reads `per_evaluated_sample` rather than remaining blank. **Do not substitute** IHECRE00001905 for IHECRE00001892 or IHECRE00003834 for IHECRE00003822 without recalculating the resulting expression medians, deltas and MAEs.

## Deploy on Yipcluster

Transfer `IHEC_EX130_pilot_ids_parallel_tracks_v3.zip` to Yipcluster and extract into your existing EX130 directory **after confirming any earlier track plotting scripts are not currently running**. The add-on does not replace EX51, EX88, or the ongoing EX130 screening/attribution code.

```bash
cd "$HOME/IHEC_Sept_2025/EX130"
unzip -o /path/to/IHEC_EX130_pilot_ids_parallel_tracks_v3.zip

# Optional: verify sample inventory; this is read-only on raw files.
python ex130_raw_tracks.py audit \
  --outdir "$HOME/IHEC_Sept_2025/EX130/results/raw_genomic_tracks/shared"

# Check in your chosen Yipcluster environment.
TRACK_PYTHON_BIN="$HOME/miniconda3/envs/DGL/bin/python3"
"$TRACK_PYTHON_BIN" -c 'import numpy, matplotlib, pyBigWig, pysam'

# Submit a single preparation job, then a 13-case high-partition
# array with FOUR simultaneously running case jobs, and a collector.
TRACK_PYTHON_BIN="$TRACK_PYTHON_BIN" MAX_CONCURRENT=4 \
  bash submit_ex130_tracks.sh
```

The array is **0-12%4**, with one case per job, 2 CPUs and 24 GB RAM per case, on the `high` partition. Preparation streams the 1.2-GB GENCODE v29 GFF3 **once** into a shared, source-identity-checked cache, and the 13 case workers read that cache. Case outputs and status receipts are kept in separate folders, preventing race conditions or overwrites. The collector waits for array termination (even if some tasks fail) and writes an honest merged status receipt.

### Optional: index the unchanged IHECRE00003822 WGBS BED without touching raw data

If you need the **whole-transcript CpG plot** for IHECRE00003822 and `tabix` is installed, use:

```bash
BUILD_WGBS_INDEX_CACHE=1 MAX_CONCURRENT=4 \
  TRACK_PYTHON_BIN="$HOME/miniconda3/envs/DGL/bin/python3" \
  bash submit_ex130_tracks.sh
```

The preparation task creates an **EX130-output symlink** to the existing BGZF-compressed whole-region BED and runs `tabix -p bed` on the symlink. The `.tbi` sidecar is created in EX130's `results/raw_genomic_tracks/shared/wgbs_index_cache`, *not* beside or inside the original raw source. This is only valid if the raw compressed BED is coordinate-sorted and BGZF, which is checked by the actual tabix exit status. If indexing fails or the binary is absent, the run records `INDEX_FAILED` or `INDEX_TOOL_UNAVAILABLE`, retains the sample and H3K36me3 tracks, and flags the CpG subplot as unindexed. Do not interpret a missing CpG subplot as zero methylation. The code never recompresses or rewrites the raw BED.

**Use the optional index-cache operation at most once per concurrent analysis run**; do not start separate prepare jobs racing to build the same sidecar. Check `results/raw_genomic_tracks/shared/index_cache_status.json` and `track_availability.csv` after preparation.

## Review outputs

- `results/raw_genomic_tracks/shared/gencode_selected_transcripts.json`: shared GENCODE coordinate cache for all cases.
- `results/raw_genomic_tracks/shared/track_availability.csv`: actual source folder, BigWig and indexed BED status for all candidates.
- `results/raw_genomic_tracks/cases/PILOT_01/figures/*.png` and `*.svg`: region tracks, exon coordinates, WGBS CpGs when indexed, and fixed-source observed/predicted expression.
- `results/raw_genomic_tracks/cases/PILOT_XX/PILOT_XX_provenance.json`: exact evaluated sample IDs and source-model ID.
- `results/raw_genomic_tracks/case_plot_status_merged.csv`: case-level status of all 13 independent tasks.
- `results/raw_genomic_tracks/track_measurement_status_merged.csv`: per-sample track decode status.
- `results/raw_genomic_tracks/track_array_receipt.json`: strict completeness receipt, non-success when any case lacks a valid result.

Each source BigWig contains the original histone p-value signal, **not normalized GCN attribution**. WGBS beta measurements are plotted only after confirming the methylation column and its fraction/percent units. The default `--select curated` preserves the three pilot samples per tissue. These tissue comparisons are illustrative, not causal perturbations.

## Validation

Run locally or on Yipcluster:

```bash
python test_raw_tracks.py
python test_parallel_tracks.py
```

The synthetic integration tests cover main/partial root fallback, preserving the unindexed raw WGBS sample, opt-in index-cache alias creation, shared GFF cache reuse, sample ID provenance, independent case plotting, fake SLURM array submission with afterok/afterany dependencies, and strict aggregation. **Actual Yipcluster decoding and running jobs are not yet verified.**
