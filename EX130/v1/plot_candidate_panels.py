#!/usr/bin/env python3
"""PNG quick-review panels for fixed-source GCN candidate transcripts.

Requires matplotlib. Outputs are descriptive, not causal/mechanistic.
No plotting is attempted until true transcript-specific attribution exists.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--outdir', type=Path, required=True)
    ap.add_argument('--max-cases', type=int, default=25)
    args = ap.parse_args()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    folder = args.outdir
    target = pd.read_csv(folder/'attribution_target_raw_features.csv', dtype={'txid':str, 'candidate_id':str})
    attrs = pd.read_csv(folder/'candidate_per_sample_assay_attribution.csv', dtype={'txid':str, 'candidate_id':str})
    destination = folder/'candidate_panel_png'
    destination.mkdir(exist_ok=True)
    case_ids = sorted(target.candidate_id.unique(), key=lambda x: (not x.startswith('PILOT_'), x))
    case_ids = case_ids[:args.max_cases]
    for case in case_ids:
        t = target[target.candidate_id == case].drop_duplicates(['txid','eval_sample','role'])
        a = attrs[attrs.candidate_id == case]
        if t.empty or a.empty:
            continue
        txid = str(t.iloc[0].txid)
        gene = t.get('symbol', pd.Series(dtype=str)).dropna().astype(str)
        gene_name = gene.iloc[0] if len(gene) else ''
        roles = ['high','low']
        fig, axs = plt.subplots(2,2, figsize=(11.5,8.2), constrained_layout=True)
        ax = axs[0,0]
        for i, role in enumerate(roles):
            z = t[t.role == role]
            obs = pd.to_numeric(z.expected_label, errors='coerce')
            pred = pd.to_numeric(z.expected_logit, errors='coerce')
            ax.scatter(np.full(len(obs),i-.09),obs,label='Observed' if i==0 else None,s=36,alpha=.85)
            ax.scatter(np.full(len(pred),i+.09),pred,label='Predicted' if i==0 else None,s=36,marker='s',alpha=.85)
        ax.set_xticks([0,1],['Higher-expression context','Lower-expression context'])
        ax.set_ylabel('log2(TPM + 0.001)')
        ax.set_title('Held-out transcript expression')
        ax.legend()
        ax = axs[0,1]
        hcol = 'h3_hg38H3K36me3intnintrons'
        for i, role in enumerate(roles):
            z = t[t.role==role]
            val = pd.to_numeric(z[hcol], errors='coerce')
            ax.scatter(np.full(len(val),i),np.log2(1+val),s=45)
        ax.set_xticks([0,1],['Higher','Lower'])
        ax.set_ylabel('log2(1 + original H3K36me3 signal)')
        ax.set_title('Original internal-intron H3K36me3')
        ax = axs[1,0]
        dcol = 'dna_beta_hg38intnintrons'
        for i, role in enumerate(roles):
            z=t[t.role==role]
            val=pd.to_numeric(z[dcol], errors='coerce')
            ax.scatter(np.full(len(val),i),val,s=45)
        ax.set_xticks([0,1],['Higher','Lower'])
        ax.set_ylim(-.05,1.05)
        ax.set_ylabel('WGBS methylation beta')
        ax.set_title('Original internal-intron DNA methylation')
        ax = axs[1,1]
        assays=['DNAm','H3K27ac','H3K27me3','H3K36me3','H3K4me1','H3K4me3','H3K9me3']
        ys=np.arange(len(assays))
        for role, offset in [('high',-.18),('low',.18)]:
            z=a[a.role==role].groupby('assay').self_gradient_abs_share.mean().reindex(assays)
            ax.barh(ys+offset,z.to_numpy(dtype=float),height=.32,label=role)
        ax.set_yticks(ys,assays)
        ax.invert_yaxis()
        ax.set_xlabel('Share of absolute target-node feature gradient')
        ax.set_title('Transcript-specific GCN input sensitivity')
        ax.legend()
        desc = ', '.join(f'{r}: {", ".join(sorted(set(t[t.role==r].tissue.astype(str))))}' for r in roles)
        fig.suptitle(f'{case}: {gene_name}  {txid}\n{desc}', fontsize=12)
        fig.text(.02,.008,'Exploratory visualization: mark contrasts and model saliency are not causal effects. GCN source model held fixed.',fontsize=8)
        fig.savefig(destination / f'{case}_{txid}.png', dpi=170)
        plt.close(fig)
    print(f'Wrote up to {len(case_ids)} candidate PNG panels to {destination}')

if __name__=='__main__':main()
