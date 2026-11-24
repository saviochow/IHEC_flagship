# EX130 WGBS percent-column corrected genomic tracks

This add-on re-renders **all 13** existing pilot examples using methylation beta equal to `percentage / 100`, with the **one-based 11th** BED column named `percentage`, rather than an ambiguous inferred column. The original plots and attribution results are unchanged. The script uses the existing DGL Python 3.8 environment and the previously prepared GENCODE annotation cache.

## Why it is necessary

The raw CpG files have tab-separated columns:

`seqnames start end strand name cov_capped start_thick end_thick color cov percentage ref_gen sample_gen quality transcript.id gene.id symbol context`

The `percentage` field is the percentage of reads methylated at that CpG (0 to 100). The figure's beta values must be `percentage * 0.01` (0 to 1). The earlier EX130 figures displayed `UNKNOWN_BETA_COLUMN` and omitted methylation points. The existing `ex130_raw_tracks_pct11.py` has correct explicit-column conversion, so the original script does **not** need replacement; this package includes a corrected-caption copy.

## Install on Yipcluster and run

Download/copy the ZIP into `~/IHEC_Sept_2025/EX130/`:

```bash
cd ~/IHEC_Sept_2025/EX130
unzip -n /path/to/IHEC_EX130_pct11_replot.zip
bash submit_ex130_pct11_replot.sh
```

This submits a preparation job, **13** per-case independent plotting jobs (array `0-12%4`), and a collector. It uses `normal`/`n001` by default to avoid the unavailable `n002` allocation in earlier runs. The retry preserves all previous figures and uses versioned output:

`results/raw_genomic_tracks_pct11/<UTC_TIMESTAMP>/cases/PILOT_XX/figures/`

A submission receipt lists the new SLURM job IDs and the exact output directory. The GENCODE cache from `results/raw_genomic_tracks/shared/gencode_selected_transcripts.json` is reused.

For all 13 cases, the prep stage checks each selected sample's gzip BED header to confirm **`percentage` is the 11th field** and checks that a real measurement is in the range 0–100. If the format differs for any sample, it refuses to misplot that case.

Optional, enabled by default, is a **local index cache** for an unindexed BED (especially IHECRE00003822). The header is un-commented, so the cache invokes `tabix -S 1 -p bed` on an EX130-local **symlink**. The raw BED is never overwritten. If `tabix` cannot create the index, the sample is retained but its CpG signal may be unavailable and is reported. Disable index cache with `EX130_PCT11_BUILD_INDEX=0`.

## Inspect and transfer results

```bash
cd ~/IHEC_Sept_2025/EX130
cat results/raw_genomic_tracks_pct11/*/submission_receipt.txt
find results/raw_genomic_tracks_pct11 -type f \( -name '*.png' -o -name '*.svg' \) | sort
cat results/raw_genomic_tracks_pct11/*/track_array_receipt.json
```

The collector writes `all_13_figures.csv`, `case_plot_status_merged.csv`, `track_measurement_status_merged.csv`, `track_array_receipt.json` and the schema audit `shared/percentage_column_audit.json`.

Each successful WGBS track should have the status `USER_SPECIFIED_COL_11` and a positive number of CpG points. Every figure is saved in PNG and SVG. An incomplete or failed case results in a non-successful collector receipt, not a false completion claim. Non-methylation warnings such as curated (rather than exact attribution-target) sample selection remain visible in the provenance.

## Not performed automatically in this session

ClusterMCP's Yipcluster scheduler gateway is disabled, so the ZIP cannot itself submit the jobs. The user must invoke the one-line submission command above on Yipcluster. No previous jobs or results are cancelled or overwritten by the new scripts.
