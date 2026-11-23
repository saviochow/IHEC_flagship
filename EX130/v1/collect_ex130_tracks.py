#!/usr/bin/env python3
"""Strictly collect separate EX130 genomic-track array outputs after completion."""
import argparse,csv,json
from pathlib import Path

def read(path):
    if not path.is_file():return None
    with path.open(newline='') as f:return list(csv.DictReader(f))

def write(path,rows,cols):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore');w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser();p.add_argument('--ex130-dir',type=Path,default=Path.home()/'IHEC_Sept_2025/EX130')
    args=p.parse_args();root=args.ex130_dir
    seed=read(root/'seed_cases.csv')
    if not seed:raise SystemExit('Missing seed_cases.csv')
    out=root/'results/raw_genomic_tracks';summaries=[];measurements=[];missing=[]
    for c in seed:
        cid=c['case_id'];folder=out/'cases'/cid
        result=read(folder/'case_plot_status.csv')
        if result is None or len(result)!=1 or result[0].get('case_id')!=cid:
            missing.append(cid)
            summaries.append({'case_id':cid,'status':'MISSING_CASE_RECEIPT','issues':'no matching case_plot_status.csv'})
            continue
        summaries.append(result[0]);measurements.extend(read(folder/'track_measurement_status.csv') or [])
    write(out/'case_plot_status_merged.csv',summaries,['case_id','status','issues'])
    write(out/'track_measurement_status_merged.csv',measurements,
          ['case_id','sample_id','role','tissue','track','status','source','points','annotation'])
    failures=[r for r in summaries if r['status'] not in ('PLOTTED','PLOTTED_WITH_WARNINGS')]
    receipt={'expected_cases':len(seed),'completed_cases':len(seed)-len(missing),
             'failed_or_incomplete':failures,'status':'COMPLETE' if not failures else 'INCOMPLETE'}
    (out/'track_array_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
    if failures:raise SystemExit(2)

if __name__=='__main__':main()
