#!/usr/bin/env python3
"""Consolidate 13 EX130 corrected genomic-track panels and CpG-beta QC receipts."""
import argparse
import csv
import json
from pathlib import Path


def rows(path):
    if not path.is_file():return None
    with path.open(newline='') as f:return list(csv.DictReader(f))


def write_csv(path, records, fields):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore')
        writer.writeheader();writer.writerows(records)


def collect(root, outdir):
    root,outdir=Path(root),Path(outdir)
    seeds=rows(root/'seed_cases.csv')
    if not seeds or len(seeds)!=13:raise RuntimeError('Exactly 13 seed cases are required')
    statuses=[];measures=[];figures=[];issues=[]
    for case in seeds:
        cid,txid=case['case_id'],case['txid']
        folder=outdir/'cases'/cid
        status=rows(folder/'case_plot_status.csv')
        track=rows(folder/'track_measurement_status.csv')
        png=folder/'figures'/'{}_{}.png'.format(cid,txid)
        svg=folder/'figures'/'{}_{}.svg'.format(cid,txid)
        good=status is not None and len(status)==1 and status[0].get('status') in ('PLOTTED','PLOTTED_WITH_WARNINGS')
        valid=[r for r in (track or []) if r['track']=='WGBS CpG beta' and r['status']=='USER_SPECIFIED_COL_11' and int(r['points'])>0]
        invalid=[r for r in (track or []) if r['track']=='WGBS CpG beta' and ('BETA_COLUMN' in r['status'] or r['status']=='NO_VALID_BETA_POINTS')]
        absent=[r for r in (track or []) if r['track']=='WGBS CpG beta' and r['status'] not in ('USER_SPECIFIED_COL_11',)]
        result={'case_id':cid,'txid':txid,'case_status':status[0]['status'] if status else 'MISSING',
                'n_cpg_valid_samples':len(valid),'n_cpg_unavailable_samples':len(absent),
                'n_cpg_invalid_samples':len(invalid),'has_png':png.is_file(),
                'has_svg':svg.is_file(),'png_path':str(png) if png.is_file() else '',
                'svg_path':str(svg) if svg.is_file() else '',
                'issues':status[0].get('issues','') if status else 'Missing status receipt'}
        statuses.append(result)
        if track:measures+=track
        if png.is_file() and svg.is_file():figures.append({'case_id':cid,'txid':txid,'png':str(png),'svg':str(svg)})
        if not(good and png.is_file() and svg.is_file() and valid and not invalid):issues.append(cid)
    write_csv(outdir/'case_plot_status_merged.csv',statuses,
        ['case_id','txid','case_status','n_cpg_valid_samples','n_cpg_unavailable_samples',
         'n_cpg_invalid_samples','has_png','has_svg','png_path','svg_path','issues'])
    write_csv(outdir/'track_measurement_status_merged.csv',measures,
        ['case_id','sample_id','role','tissue','track','status','source','points','annotation'])
    write_csv(outdir/'all_13_figures.csv',figures,['case_id','txid','png','svg'])
    proof={'expected_cases':len(seeds),'cases_with_png_svg':len(figures),
           'cases_with_valid_CpG_beta':sum(int(x['n_cpg_valid_samples'])>0 for x in statuses),
           'n_valid_sample_CpG_tracks':sum(int(x['n_cpg_valid_samples']) for x in statuses),
           'n_unavailable_sample_CpG_tracks':sum(int(x['n_cpg_unavailable_samples']) for x in statuses),
           'failed_cases':issues,'beta_column_one_based':11,
           'beta_definition':'methylated_reads_percent / 100',
           'status':'COMPLETE' if not issues and len(figures)==len(seeds) else 'INCOMPLETE'}
    (outdir/'track_array_receipt.json').write_text(json.dumps(proof,indent=2)+'\n')
    print(json.dumps(proof,indent=2))
    return proof


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ex130-dir',type=Path,default=Path.home()/'IHEC_Sept_2025/EX130')
    parser.add_argument('--outdir',type=Path,required=True)
    a=parser.parse_args()
    if collect(a.ex130_dir,a.outdir)['status']!='COMPLETE':raise SystemExit(2)


if __name__=='__main__':main()
