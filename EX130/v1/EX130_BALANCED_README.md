# EX130 Yipcluster balanced Slurm submission (v4)

This add-on submits **existing EX130 workloads** across three explicitly pinned
Yipcluster nodes. It does **not** change the original screens, GCN models,
transcript targets, source graph files, or raw epigenomic tracks. It uses the
existing `/home/hchow/miniconda3/envs/DGL/bin/python3` interpreter.

## Placement and load

Assumed physical capacity, subject to Slurm confirmation at submission:

| Partition | Node | Reported cores | Attribution tasks (of 69) | Track cases (of 13) |
|---|---|---:|---:|---:|
| high | n003 | 192 | 30 | 6 |
| normal | n001 | 128 | 20 | 4 |
| normal | n002 | 128 | 19 | 3 |

Indices are assigned by `index % 7`: residues 0/1/2 to n003, 3/4 to n001,
5/6 to n002. This produces **30 high / 39 normal attribution tasks** and
**6 high / 7 normal track cases**, matching 192:256 total partition cores as
closely as integer task counts allow. Task manifests and their row numbers are
left unchanged.

The default `MAX_CONCURRENT=4` is **per workflow**: each attribution or track
workflow runs at most 2 jobs on n003, 1 on n001 and 1 on n002, subject to Slurm
admission. Submitting `all` can therefore run up to 8 total jobs (4 attribution,
4 plotting) at once, **not** 4 total across both workflows. Each attribution
job still requests the original 8 CPUs/96 GB; each track job requests the
original 2 CPUs/24 GB. Core counts alone do not prove memory headroom.
Increase concurrency only after checking `sinfo` and node memory. You can set
`ATTR_MAX_CONCURRENT=8` or `TRACK_MAX_CONCURRENT=8` separately.

## Installation

Unzip the add-on **directly into** `~/IHEC_Sept_2025/EX130`. These are new
files, not replacements for the old scientific scripts:

```bash
cd ~/IHEC_Sept_2025/EX130
unzip /path/to/IHEC_EX130_balanced_high_normal_v4.zip -d .

# Check the physical node / partition membership first.
sinfo -N -p high,normal -o '%P %N %c %m %t'

# Fix the already diagnosed Python 3.8-only argparse issue (idempotent).
/usr/bin/python3 patch_ex130_python38.py
```

If this check says the source is already patched, no change is made.
A backup of the original `transcript_attribution.py` is created for first-time
patches. The wrapper refuses to schedule attribution tasks when the Python
3.9-only `BooleanOptionalAction` expression remains.

**Track prerequisite:** The existing DGL environment must be able to import
`pyBigWig`, `pysam`, `numpy`, and `matplotlib`. The older preparation job failed
because `pyBigWig` was missing. If required, add only readers to the **existing
DGL** environment (no new Conda env):

```bash
export PYTHON_BIN=/home/hchow/miniconda3/envs/DGL/bin/python3
"$PYTHON_BIN" -m pip install --only-binary=:all: --no-deps \
    'pyBigWig==0.3.22' 'pysam==0.22.1'
```

This does not upgrade Python, DGL, PyTorch, or the NumPy installation.
Do not resubmit the plotting arrays if the preparation imports still fail.

## Dry run, then submission

```bash
cd ~/IHEC_Sept_2025/EX130

# No Slurm state mutation. Check task distribution and dependencies.
MAX_CONCURRENT=4 bash submit_ex130_balanced.sh all --dry-run

# New independent attribution and track jobs; old failed jobs are not retried.
MAX_CONCURRENT=4 bash submit_ex130_balanced.sh all --submit
```

You may submit the two workflows separately:

```bash
MAX_CONCURRENT=4 bash submit_ex130_balanced.sh attribution --submit
MAX_CONCURRENT=4 bash submit_ex130_balanced.sh tracks --submit
```

No preparation or GCN retraining is repeated. The existing manifest drives
69 per-transcript tasks; 13 seed cases drive raw-track figures. The prior
`--skip-completed` attribution behavior remains intact. The script does not
automatically cancel any previous pending jobs. If old jobs `105859`, `105882`,
or `105883` are still stuck with unsatisfiable dependencies, cancel only those
specific jobs after checking with `squeue`.

The wrapper checks live membership of n003 in `high` and n001/n002 in
`normal` before submission. It refuses to guess if your cluster's
partition-to-node mapping differs. Existing `#SBATCH --partition=high`
directives are overridden by the explicit `sbatch --partition` CLI options.
`--nodelist` and `--nodes=1` pin each task to exactly one selected node.

## Dependencies and results

* Attribution: three node-pinned arrays, then one summary under `normal/n001`
  with `afterok` dependencies on all three arrays. A failed task cancels the
  dependent summary rather than leaving an eternally pending job.
* Tracks: one shared GENCODE/index-cache preparation on `normal/n001`; then
  three node-pinned arrays with `afterok` on the preparation; then one
  collection job on `normal/n002` with `afterany` on all three arrays.
  A failed preparation cancels the downstream invalid arrays; the collector
  can report missing receipts.
* Submission IDs are recorded in
  `results/balanced_submissions/<UTC-timestamp>_<mode>_jobs.txt`.

Monitor with:

```bash
squeue -u "$USER" -o '%.19i %.20j %.10P %.10T %.5C %R'
ls -ltr results/balanced_submissions
find results/fixed_source/attribution_by_sample -name '*.provenance.json' 2>/dev/null | wc -l
cat results/fixed_source/attribution_analysis_summary.json
cat results/raw_genomic_tracks/track_array_receipt.json
```

Look for 69 attribution receipts with matching targets/checkpoint provenance,
33 intended candidate cases in the attribution manifest, and complete status
for 13 track figures. **Do not infer successful scientific results merely from
submission or completed Slurm job states.**

## Verification

`test_balanced_submission.py` uses temporary fixtures and mocked `sinfo` and
`sbatch`, never the actual scheduler. It verifies all indices are assigned
once, job counts, throttles, node pinning, both dependency fan-ins, early
Python 3.8 rejection and idempotent patching. This is NOT a live acceptance
result on Yipcluster.
