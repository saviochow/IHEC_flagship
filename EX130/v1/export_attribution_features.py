#!/usr/bin/env python3
"""Read original WGBS/H3K36me3 values for all targeted EX130 cases, including seeds.

Reads only selected transcript rows from original 60-110 MB tables, in chunks.
This is separate from GCN input attribution: raw marks are descriptive and not
interchangeable with the standardized 147-dimensional DGL input matrix.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import pandas as pd

from screen_ihec_cases import extract_feature_rows, DNAM_COLUMNS, H3_COLUMNS


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--outdir', type=Path, required=True)
    ap.add_argument('--features-root', type=Path, default=Path('/home/shared/IHEC/baseline_expt/features'))
    ap.add_argument('--allow-missing', action='store_true', help='Diagnostics only; missing features must be excluded from publication')
    args = ap.parse_args()
    targets = pd.read_csv(args.outdir / 'attribution_targets.csv', dtype=str)
    grouped = targets.groupby('base_id').txid.agg(lambda x: set(x)).to_dict()
    data, issues = [], []
    for i, (base, txids) in enumerate(sorted(grouped.items()), 1):
        directory = args.features_root / base
        h = list(directory.glob(f'{base}.[0-9]*.H3K36me3.pval.all.txt'))
        m = list(directory.glob(f'{base}.[0-9]*.standard-WGBS.cpg.all.txt'))
        if len(h) != 1 or len(m) != 1:
            issues.append(f'{base}: ambiguous or missing feature files; H3={len(h)}, WGBS={len(m)}')
            continue
        try:
            h3 = extract_feature_rows(h[0], txids, H3_COLUMNS, 'h3_')
            wgbs = extract_feature_rows(m[0], txids, DNAM_COLUMNS, 'dna_')
            wgbs = wgbs.drop(columns=['gene.id', 'symbol'], errors='ignore')
            joined = h3.merge(wgbs, on='txid', how='outer', validate='one_to_one')
            joined['base_id'] = base
            data.append(joined)
        except (OSError, ValueError, KeyError) as exc:
            issues.append(f'{base}: {exc}')
        if i % 20 == 0:
            print(f'Read raw feature tables for {i}/{len(grouped)} biological specimens', flush=True)
    if not data:
        raise RuntimeError(f'No original features retrieved: {issues[:5]}')
    raw = pd.concat(data, ignore_index=True)
    out = targets.merge(raw, on=['base_id', 'txid'], how='left', validate='many_to_one')
    hc='h3_hg38H3K36me3intnintrons'
    dc='dna_beta_hg38intnintrons'
    if hc not in out.columns or dc not in out.columns:
        raise RuntimeError('Missing required internal-intron features')
    for col in (hc,dc):out[col] = pd.to_numeric(out[col], errors='coerce')
    out['has_h3_internal_introns'] = out[hc].ge(0) & out[hc].notna()
    out['has_dna_internal_introns'] = out[dc].between(0,1) & out[dc].notna()
    out['has_both_internal_introns'] = out.has_h3_internal_introns & out.has_dna_internal_introns
    out.to_csv(args.outdir / 'attribution_target_raw_features.csv', index=False)
    if not out.has_both_internal_introns.all():
        issues.append(f'{int((~out.has_both_internal_introns).sum())} target rows lack one or both internal-intron measurements; see attribution_target_raw_features.csv')
    (args.outdir / 'raw_feature_issues.json').write_text(json.dumps(issues, indent=2))
    print(f'Raw features joined: {out.has_both_internal_introns.sum()}/{len(out)} rows have both internal-intron measurements; {len(issues)} issues')
    if issues and not args.allow_missing:
        print('WARNING: missing features in some targets; attribution may be valid for measured transcripts but case panels should exclude missing rows')

if __name__ == '__main__':
    main()
