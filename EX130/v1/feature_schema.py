"""Explicit EX18 GCN feature positions; *not* the EX88 RF plot index order.

Source inspected: EX18/data_preproc_regression_all_tpmpc1e-3_hmpc1e-9_zscale_clip5_WGBS+HM.py
Constructed for the test node features in that file, not inferred from the published plot.
The live DGL graph's first 21 columns are independently checked against WGBS
by transcript_attribution.py before attribution is accepted.
"""
from __future__ import annotations

REGIONS = (
    "Up1", "Up2", "Up3", "Up4", "Up5", "TSS1", "TSS2", "TSS3", "TSS4", "TSS5",
    "firstexon", "firstintron", "intnexons", "intnintrons", "lastintron", "lastexon",
    "Dw1", "Dw2", "Dw3", "Dw4", "Dw5",
)
# EX18 graph construction: WGBS first, then histone marks.
ASSAYS = ("DNAm", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3")
# The positional names loaded in the existing EX88 plot, for audit only.
EX88_REFERENCE_ASSAYS = ("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3", "DNAm")


def mapping():
    return [
        {
            "feature_index": ai * len(REGIONS) + ri,
            "assay": assay,
            "region": region,
            "raw_feature_column": ("beta_hg38" if assay == "DNAm" else "hg38" + assay) + region,
            "ex88_reference_assay_at_this_index": EX88_REFERENCE_ASSAYS[ai],
            "index_label_differs_from_ex88_reference": assay != EX88_REFERENCE_ASSAYS[ai],
        }
        for ai, assay in enumerate(ASSAYS)
        for ri, region in enumerate(REGIONS)
    ]


def write_mapping(path):
    import csv
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(mapping()[0]))
        w.writeheader()
        w.writerows(mapping())


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True)
    args = p.parse_args()
    write_mapping(args.out)
