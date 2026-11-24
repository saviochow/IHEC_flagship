# IHEC_flagship | EX130 reproducible pilot analyses

This repository preserves the original EX130 reviewer pilot code in `EX130/v1/` and adds an exploratory **within-gene isoform-switch** analysis in `EX130/v2/`.

**v1:** transcription prediction screens, per-transcript gradient attribution, raw genomic BigWig and whole-genome WGBS BED plotting, explicitly corrected WGBS 11th-column methylation = percentage / 100, provenance manifests, Slurm wrappers and tests. Previously generated PNGs and private study data are not committed.

**v2:** in each biological sample, compare two transcripts from the same gene. Identify paired tissue groups in which their within-sample expression order reverses; independently assess DNA methylation and H3K36me3 contrasts in matched samples. Include stable-expression negative controls and mark missing assays as missing, not zero. The discovery screen is expression-only; epigenomic measurements are validation, not selection. These analyses are descriptive, not causal.

V2 code is in a **content-hashed deterministic source archive** committed in text-safe chunks (`EX130/v2/bundle/00.b64` ... `10.b64`), unpacked by `scripts/bootstrap_savio_ex130_v2.sh`. This keeps source and tests reproducible through ordinary GitHub without storing private epigenomic datasets. Archive SHA256: `0aa9ed114e6f0c89077d8f2490f776c8f10b03152158cc93a6825edaa59bcdea`.

On saviocluster, supply an exact pinned commit SHA:
```bash
export IHEC_PINNED_SHA='<40-character commit SHA of this branch>'
git clone --branch research/ex130-isoform-switch-v2-20261123 https://github.com/saviochow/IHEC_flagship.git
bash IHEC_flagship/scripts/bootstrap_savio_ex130_v2.sh
```

The bootstrap creates a new, isolated checkout; checks source identity and the archive digest; and runs synthetic unit/integration tests. It **does not submit Slurm jobs**. An authenticated ClusterMCP connection and independently APPROVED Researchstack plan with a one-use approval are required for Yipcluster dispatch. The broker was disabled and unreachable at the time of publication; no Yipcluster v2 execution is claimed. Full methods, safeguards, command examples, and output schema are inside the v2 archive (`EX130/v2/README.md`).

**Repository hygiene:** never add dataset-level features, human WGBS BEDs, BigWigs, intermediate expression tables, model checkpoints, lab-private Slurm output, SSH credentials, tokens, or source-machine directories as binaries.

## MCP client compatibility repair (2026-11-23)

The frozen v2 research archive has a historical low-level MCP client whose HTTP arguments are not accepted by the maintained SDK. To preserve the exact staged research source, run the new isolated [`scripts/mcp_dispatch_fixed.py`](scripts/mcp_dispatch_fixed.py) with the pinned `mcp==1.30.0` environment. This fix is independently CI-tested; consult [`EX130/v2/MCP_CLIENT_SETUP.md`](EX130/v2/MCP_CLIENT_SETUP.md). Scientific code and Yipcluster work-root SHA remain unchanged. The scheduler gateway remains disabled until separate operator/Researchstack approval.
