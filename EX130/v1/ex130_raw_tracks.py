#!/usr/bin/env python3
"""EX130 genomic-track reviewer panels from GENCODE v29, ChIP BigWig and indexed CpG BED.

No input modifications.  ``audit`` uses stdlib only; ``plot`` additionally uses
numpy/matplotlib, pyBigWig and pysam. Missing/ambiguous evidence is reported,
not substituted with fabricated signal. RNA protocol and GCN source remain fixed
within a case. Coordinates: GFF3 1-based closed -> 0-based half-open.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import unquote

CASE_FIELDS = ('case_id','txid','rna_type','tissue_high','tissue_low','source_model')
MARKS = ('H3K27ac','H3K27me3','H3K36me3','H3K4me1','H3K4me3','H3K9me3')
SAMPLE_RE = re.compile(r'^(IHECRE\d+)(?:_(mRNA-Seq|total-RNA-Seq))?$')


def csv_rows(path, required=()):
    with Path(path).open(newline='',encoding='utf-8-sig') as f:
        reader=csv.DictReader(f)
        if reader.fieldnames is None or not set(required).issubset(reader.fieldnames):
            raise ValueError(f'{path}: missing CSV columns {sorted(set(required)-set(reader.fieldnames or []))}')
        return list(reader)


def write_rows(path, rows, columns):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore');wr.writeheader();wr.writerows(rows)


def base_id(sample):
    m=SAMPLE_RE.fullmatch(str(sample))
    if not m: raise ValueError(f'Unsafe/unrecognized IHEC sample name: {sample!r}')
    return m.group(1)


def single_glob(directory,pat):
    files=sorted(x for x in Path(directory).glob(pat) if x.is_file())
    return (files[0], 'OK') if len(files)==1 else (None, 'MISSING' if len(files)==0 else 'AMBIGUOUS')


def locate_tracks(sample, root, marks, signal='pval', secondary_root=None,
                  index_cache_dir=None):
    """Find original raw tracks without conflating a missing .tbi with missing BED data.

    A specimen can live in the main raw tree or the partial-epigenomes tree.
    Do not mix files from two source roots for a single specimen. Indexed
    aliases under EX130 output are optional and never alter the source BED.
    """
    base=base_id(sample)
    locations=[Path(root)/base]
    if secondary_root is not None:
        other=Path(secondary_root)/base
        if other!=locations[0]:locations.append(other)
    available=[x for x in locations if x.is_dir()]
    d={'sample_id':sample,'base_id':base,'raw_dir':'', 'mark_paths':{},'wgbs_path':'',
       'wgbs_index':'','wgbs_protocol':'','status':'OK','issues':[]}
    if not available:
        d['raw_dir']=';'.join(str(x) for x in locations)
        d['status']='MISSING';d['issues'].append('RAW_SAMPLE_DIRECTORY_ABSENT');return d
    # Prefer a root with both the histone track and WGBS BED. Main wins ties.
    def score(path):
        hits=sum(bool(list((path/'ChIP-Seq').glob(f'{base}.[0-9]*.{m}.{signal}.signal.bigwig'))) for m in marks)
        bed=bool(list((path/'WGBS').glob(f'{base}.[0-9]*.*WGBS.cpg.all.bed.gz')))
        return (int(bed and (hits==len(marks))), hits+int(bed))
    path=max(available,key=score)
    d['raw_dir']=str(path)
    for mark in marks:
        p,ok=single_glob(path/'ChIP-Seq',f'{base}.[0-9]*.{mark}.{signal}.signal.bigwig')
        d['mark_paths'][mark]=str(p) if p else ''
        if ok!='OK': d['issues'].append(f'{mark}_{signal}_BIGWIG_{ok}')
    p,ok=single_glob(path/'WGBS',f'{base}.[0-9]*.*WGBS.cpg.all.bed.gz')
    if p is None:
        d['issues'].append(f'WHOLE_REGION_CPG_BED_{ok}')
    else:
        d['wgbs_path']=str(p)
        d['wgbs_protocol']='PBAT-WGBS' if '.PBAT-WGBS.' in p.name else ('standard-WGBS' if '.standard-WGBS.' in p.name else 'OTHER')
        index_candidates=[Path(str(p)+'.tbi'),Path(str(p)+'.csi')]
        if index_cache_dir is not None:
            alias=Path(index_cache_dir)/base/p.name
            # Only use aliases which point to this original immutable BED.
            if alias.is_symlink() and alias.resolve()==p.resolve():
                index_candidates += [Path(str(alias)+'.tbi'),Path(str(alias)+'.csi')]
                if any(x.is_file() for x in index_candidates[2:]):
                    d['wgbs_path']=str(alias)
        for index in index_candidates:
            if index.is_file():
                d['wgbs_index']=str(index);break
        if not d['wgbs_index']:
            d['issues'].append('WHOLE_REGION_CPG_INDEX_MISSING')
    if not d['issues']:
        d['status']='OK'
    elif d['issues']==['WHOLE_REGION_CPG_INDEX_MISSING']:
        d['status']='RAW_PRESENT_INDEX_PENDING'
    else:
        d['status']='INCOMPLETE'
    return d


def make_sample_audit(samples, root, marks, signal, secondary_root=None,
                      index_cache_dir=None):
    seen={};out=[]
    for row in samples:
        sample=row['sample_id']
        if sample in seen:raise ValueError(f'Duplicate sample: {sample}')
        seen[sample]=True
        found=locate_tracks(sample,root,marks,signal,secondary_root,index_cache_dir)
        out.append({**row,**{k:v for k,v in found.items() if k not in ('mark_paths',)},
            **{f'{m}_{signal}_bigwig':found['mark_paths'].get(m,'') for m in marks},
            'issues':';'.join(found['issues'])})
    return out


def build_cached_indexes(audit, cache_dir, tabix_bin='tabix'):
    """Optionally create index sidecars on EX130 symlinks, never on raw data.

    This relies on existing BED being coordinate-sorted and BGZF compressed.
    A nonzero tabix exit is recorded; the original BED is never recompressed.
    """
    import subprocess
    cache_dir=Path(cache_dir)
    outcomes=[]
    for row in audit:
        if row['status']!='RAW_PRESENT_INDEX_PENDING':continue
        raw=Path(row['wgbs_path'])
        if not raw.is_file():continue
        alias_dir=cache_dir/row['base_id']
        alias_dir.mkdir(parents=True,exist_ok=True)
        alias=alias_dir/raw.name
        if alias.is_symlink() and alias.resolve()!=raw.resolve():
            raise RuntimeError(f'Refusing mismatched existing index-cache symlink: {alias}')
        if alias.exists() and not alias.is_symlink():
            raise RuntimeError(f'Refusing non-symlink index-cache file: {alias}')
        if not alias.is_symlink():alias.symlink_to(raw.resolve())
        if any(Path(str(alias)+suffix).exists() for suffix in ('.tbi','.csi')):
            outcomes.append({'sample_id':row['sample_id'],'status':'CACHE_ALREADY_INDEXED','path':str(alias)})
            continue
        # subprocess.run does not mutate the raw file, only creates alias.tbi.
        try:
            proc=subprocess.run([tabix_bin,'-p','bed',str(alias)],capture_output=True,text=True)
            result='CACHED_INDEX_CREATED' if proc.returncode==0 and Path(str(alias)+'.tbi').is_file() else 'INDEX_FAILED'
            message=proc.stderr.strip()[:400]
        except OSError as exc:
            result='INDEX_TOOL_UNAVAILABLE';message=str(exc)[:400]
        outcomes.append({'sample_id':row['sample_id'],'status':result,'path':str(alias),
                         'message':message})
    return outcomes


def gff_attributes(field):
    result={}
    for part in field.split(';'):
        if '=' not in part:continue
        key,val=part.split('=',1)
        result[key]=unquote(val)
    return result


def find_annotations(path,txids):
    """Single streaming pass through ~1.2 GB GENCODE v29 GFF3 for exact IDs."""
    targets=set(txids);records={};exons=defaultdict(list)
    # Early reject lines without any requested ENST prefix; do not assume fixed exon order.
    roots={t.split('.')[0] for t in targets}
    rx=re.compile('|'.join(map(re.escape,sorted(roots))))
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt',encoding='utf-8',errors='replace') as f:
        for line in f:
            if line.startswith('#') or not rx.search(line):continue
            cols=line.rstrip('\n').split('\t')
            if len(cols)!=9:continue
            chrom,src,feature,start,end,_,strand,_,attrs=cols
            if feature not in ('transcript','mRNA','exon','CDS'):continue
            a=gff_attributes(attrs)
            if feature in ('transcript','mRNA'):
                cand=(a.get('ID') or a.get('transcript_id') or '').split(':')[-1]
                if cand in targets:
                    records[cand]={'txid':cand,'chr':chrom,'start0':int(start)-1,
                        'end0':int(end),'strand':strand,'gene':a.get('gene_name') or a.get('Name') or '',
                        'gene_id':a.get('gene_id') or a.get('Parent') or '', 'exons':[]}
            elif feature=='exon':
                parent=a.get('Parent') or a.get('transcript_id') or ''
                for id0 in parent.split(','):
                    id0=id0.split(':')[-1]
                    if id0 in targets:exons[id0].append([int(start)-1,int(end)])
    for tx,rec in records.items():rec['exons']=sorted(set(tuple(x) for x in exons.get(tx,[])))
    return records


def annotation_cache(gff, txids, outdir):
    cache=Path(outdir)/'gencode_selected_transcripts.json'
    src=Path(gff).stat()
    identity={'file':str(Path(gff).resolve()),'size':src.st_size,'mtime_ns':src.st_mtime_ns,'txids':sorted(txids)}
    if cache.exists():
        old=json.loads(cache.read_text())
        ident=old.get('identity') or {}
        # A prepared all-pilots cache must be reusable by each independent
        # one-case SLURM array task without concurrent 1.2 GB re-scans.
        if (ident.get('file')==identity['file'] and ident.get('size')==identity['size']
            and ident.get('mtime_ns')==identity['mtime_ns']
            and set(txids).issubset(set(ident.get('txids',[])))):
            return {k:v for k,v in old['transcripts'].items() if k in txids}
    print('Scanning GENCODE GFF3 once for selected transcripts (file may exceed 1 GB)',flush=True)
    records=find_annotations(gff,txids)
    cache.parent.mkdir(parents=True,exist_ok=True)
    cache.write_text(json.dumps({'identity':identity,'transcripts':records},indent=2)+'\n')
    return records


def choose_samples(case, samples, targets, mode='auto', max_per_role=3):
    out={};sources={}
    for role in ('high','low'):
        tissue=case['tissue_'+role]
        options=[r for r in samples if r['rna_type']==case['rna_type'] and r['tissue']==tissue]
        if mode!='curated':
            from_target=[r for r in targets if r['candidate_id']==case['case_id'] and
                         r['role']==role and r['txid']==case['txid'] and
                         r['source_model']==case['source_model']]
            ids={r['eval_sample'] for r in from_target}
            chosen=[r for r in options if r['sample_id'] in ids]
        else:chosen=[]
        if chosen:
            sources[role]='EX130_ATTRIBUTION_TARGETS'
        elif mode=='targets':
            raise ValueError(f'No EX130 attribution targets for {case["case_id"]} {role}')
        else:
            chosen=options;sources[role]='CURATED_TRACK_SAMPLE_MANIFEST'
        out[role]=chosen[:max_per_role]
        if not out[role]:raise ValueError(f'No eligible samples for {case["case_id"]}, {tissue} {case["rna_type"]}')
    return out,sources


def read_pred(sample,source,txid,logits_dir):
    p=Path(logits_dir)/f'10sample_+thyroid_{source}_{sample}_128_test_logits.csv'
    if not p.is_file():return {'status':'PREDICTION_CSV_MISSING'}
    with p.open(newline='') as f:
        r=csv.DictReader(f)
        if not {'txid','labels','logits','chr'}.issubset(set(r.fieldnames or [])):
            return {'status':'INVALID_LOGITS_HEADER'}
        vals=[x for x in r if x['txid']==txid and x['chr'] in ('chr18','chr19')]
    if len(vals)!=1:return {'status':f'TRANSCRIPT_ROWS_{len(vals)}'}
    obs,pred=float(vals[0]['labels']),float(vals[0]['logits'])
    return {'status':'OK','obs':obs,'pred':pred,'obs_tpm':max(0,2**obs-0.001)}


def chrom_alias(chrom, options):
    if chrom in options:return chrom
    alt=chrom[3:] if chrom.startswith('chr') else 'chr'+chrom
    if alt in options:return alt
    raise ValueError(f'Chromosome {chrom!r} not found; available sample: {list(options)[:8]}')


def bw_signal(path, chrom, start0, end0, bins):
    import pyBigWig
    bw=pyBigWig.open(str(path))
    if bw is None:raise ValueError(f'Could not open BigWig: {path}')
    try:
        contig=chrom_alias(chrom,bw.chroms())
        vals=bw.stats(contig,start0,end0,nBins=bins,type='mean')
    finally:
        bw.close()
    import numpy as np
    return np.asarray([float('nan') if x is None else float(x) for x in vals],dtype=float)


def wgbs_rows(path, chrom, start0, end0, max_points=50000):
    import pysam
    with pysam.TabixFile(str(path)) as tab:
        chrom=chrom_alias(chrom,tab.contigs)
        records=[]
        for line in tab.fetch(chrom,start0,end0):
            p=line.split('\t')
            if len(p)<4:continue
            try:pos=int(p[1]);end=int(p[2])
            except ValueError:continue
            if end<=start0 or pos>=end0:continue
            records.append((pos,p))
            if len(records)>=max_points:break
    return records


def infer_beta_column(records, selected='auto', scale='fraction'):
    """Do not silently interpret coverage/counts as methylation beta."""
    if not records:return None,'NO_CPG_IN_REGION'
    nfields=min(len(fields) for _,fields in records)
    if selected!='auto':
        col=int(selected)-1
        if col<3 or col>=nfields:return None,'INVALID_EXPLICIT_BETA_COLUMN'
        return col,'USER_SPECIFIED'
    if nfields==4:
        samples=[]
        for _,fields in records[:min(300,len(records))]:
            try: samples.append(float(fields[3]))
            except ValueError: pass
        ceiling=100 if scale=='percent' else 1
        if samples and len(samples)>=max(1,int(.9*min(300,len(records)))) and all(math.isfinite(v) and 0<=v<=ceiling for v in samples):
            return 3,'AUTO_FOUR_COLUMN_BED'
        return None,'INVALID_OR_UNCONFIRMED_BETA_UNITS'
    # For >=5 fields, only infer if exactly one non-coordinate column always
    # falls within beta fraction range (0..1).  Otherwise require an override.
    plausible=[]
    for col in range(3,nfields):
        nums=[]
        for _,fields in records[:min(300,len(records))]:
            try:nums.append(float(fields[col]))
            except (ValueError,IndexError):pass
        if nums and len(nums)>=max(1,int(.9*min(300,len(records)))) and all(math.isfinite(v) and 0<=v<=1 for v in nums):
            plausible.append(col)
    if scale=='fraction' and len(plausible)==1:return plausible[0],'AUTO_UNIQUE_BETA_FRACTION'
    return None,'AMBIGUOUS_BETA_COLUMN' if len(plausible)>1 else 'UNKNOWN_BETA_COLUMN'


def beta_values(records, col, scale):
    import numpy as np
    if col is None:return np.array([]),np.array([])
    xy=[]
    for pos,fields in records:
        try:v=float(fields[col])*(.01 if scale=='percent' else 1)
        except (ValueError,IndexError):continue
        if math.isfinite(v) and 0<=v<=1:xy.append((pos+.5,v))
    if not xy:return np.array([]),np.array([])
    xy=np.asarray(xy,dtype=float)
    return xy[:,0],xy[:,1]


def plot_case(case, samples, source_tags, ann, args, mark_list):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    tx=case['txid'];start0=max(0,ann['start0']-args.flank);end0=ann['end0']+args.flank
    if end0<=start0:raise ValueError('Invalid genomic annotation interval')
    bins=min(args.bins,end0-start0)
    xs=np.linspace(start0,end0,bins,endpoint=False)+(end0-start0)/(2*bins)
    ntracks=len(mark_list)+1
    figh=3.7+2*ntracks*1.35
    fig,axes=plt.subplots(2+2*ntracks,1,figsize=(14,figh),sharex=False,
                          gridspec_kw={'height_ratios':[.7]+[1]*(2*ntracks)+[1.45]},constrained_layout=True)
    for ax in axes[1:-1]:ax.sharex(axes[0])
    axes[0].set_xlim(start0,end0);axes[0].plot([ann['start0'],ann['end0']],[.5,.5],color='0.4',lw=1)
    for a,b in ann['exons']:axes[0].add_patch(plt.Rectangle((a,.36),b-a,.28,color='0.25'))
    axes[0].set_ylim(0,1);axes[0].set_yticks([])
    axes[0].set_ylabel('Gencode v29\ntranscript')
    axes[0].set_title(f"{ann['gene']}  {tx}  ({ann['chr']}:{ann['start0']+1:,}-{ann['end0']:,}, {ann['strand']})",loc='left',fontsize=10)
    colors=plt.rcParams['axes.prop_cycle'].by_key()['color']
    row_reports=[];missing=[];overall_max={mark:0. for mark in mark_list}
    mean_tracks={}
    for rix,role in enumerate(('high','low')):
        grp=samples[role]; tissue=case['tissue_'+role]
        for j,mark in enumerate(mark_list):
            ax=axes[1+rix*ntracks+j];vals_all=[]
            for sid,entry in enumerate(grp):
                sample=entry['sample_id'];found=locate_tracks(sample,args.raw_root,[mark],args.signal,args.partial_root,args.index_cache_dir)
                file=found['mark_paths'].get(mark)
                if not file:missing.append(f'{sample}: {mark} BigWig missing');continue
                try:
                    vals=bw_signal(file,ann['chr'],start0,end0,bins)
                except Exception as ex:
                    missing.append(f'{sample}: BigWig decode failure: {type(ex).__name__}: {ex}');continue
                vals_all.append(vals);ax.plot(xs,vals,alpha=.48,lw=.75,label=sample)
                finite=vals[np.isfinite(vals)]
                if len(finite):overall_max[mark]=max(overall_max[mark],float(np.max(finite)))
                row_reports.append({'case_id':case['case_id'],'sample_id':sample,'role':role,'tissue':tissue,
                    'track':mark+'_'+args.signal,'status':'OK','source':file,
                    'points':int(np.isfinite(vals).sum()),'annotation':ann['txid']})
            if vals_all:
                with np.errstate(all='ignore'):
                    arr=np.asarray(vals_all)
                    num=np.isfinite(arr).sum(axis=0)
                    avg=np.where(num>0,np.nansum(arr,axis=0)/np.maximum(1,num),np.nan)
                ax.plot(xs,avg,lw=1.7,color='black',label='mean (available samples)')
            ax.set_ylabel(f'{tissue}\n{mark}\n{args.signal}.signal',fontsize=8)
            ax.legend(loc='upper right',fontsize=6,ncol=min(4,len(grp)+1))
            mean_tracks[(role,mark)]=ax
        ax=axes[1+rix*ntracks+len(mark_list)]
        for idx,entry in enumerate(grp):
            sample=entry['sample_id'];found=locate_tracks(sample,args.raw_root,[],args.signal,args.partial_root,args.index_cache_dir)
            path=found['wgbs_path']
            if not path or not found['wgbs_index']:
                missing.append(f'{sample}: whole-region CpG BED exists but tabix index is unavailable' if path else
                               f'{sample}: whole-region CpG BED missing')
                row_reports.append({'case_id':case['case_id'],'sample_id':sample,'role':role,
                    'tissue':tissue,'track':'WGBS CpG beta',
                    'status':'RAW_PRESENT_INDEX_PENDING' if path else 'RAW_BED_MISSING',
                    'source':path,'points':0,'annotation':ann['txid']})
                continue
            try:
                recs=wgbs_rows(path,ann['chr'],start0,end0,args.max_cpg)
                beta_col,why=infer_beta_column(recs,args.wgbs_beta_column,args.wgbs_scale)
                if beta_col is None:
                    missing.append(f'{sample}: {why}; rerun with --wgbs-beta-column N after examining BED fields')
                    row_reports.append({'case_id':case['case_id'],'sample_id':sample,'role':role,'tissue':tissue,
                        'track':'WGBS CpG beta','status':why,'source':path,'points':len(recs),'annotation':ann['txid']})
                    continue
                xx,yy=beta_values(recs,beta_col,args.wgbs_scale)
                if not len(xx):
                    missing.append(f'{sample}: no valid CpG beta fractions in selected BED column')
                    row_reports.append({'case_id':case['case_id'],'sample_id':sample,'role':role,'tissue':tissue,
                        'track':'WGBS CpG beta','status':'NO_VALID_BETA_POINTS','source':path,
                        'points':0,'annotation':ann['txid']})
                    continue
                if len(xx)>args.display_cpg:
                    steps=np.linspace(0,len(xx)-1,args.display_cpg,dtype=int)
                    xx=xx[steps];yy=yy[steps]
                ax.scatter(xx,yy,s=3,alpha=.24,color=colors[idx%len(colors)],label=sample)
                row_reports.append({'case_id':case['case_id'],'sample_id':sample,'role':role,'tissue':tissue,
                    'track':'WGBS CpG beta','status':why+f'_COL_{beta_col+1}',
                    'source':path,'points':len(xx),'annotation':ann['txid']})
            except Exception as ex:
                missing.append(f'{sample}: CpG decode failure: {type(ex).__name__}: {ex}')
        ax.set_ylim(-.04,1.04)
        ax.set_ylabel(f'{tissue}\nCpG beta',fontsize=8)
        ax.legend(loc='upper right',fontsize=6,ncol=min(4,len(grp)))
    for mark in mark_list:
        maxy=overall_max[mark]
        for role in ('high','low'):
            mean_tracks[(role,mark)].set_ylim(bottom=0,top=max(1.05*maxy,.01))
    axes[-2].xaxis.set_major_formatter(FuncFormatter(lambda v,_:f'{v/1e6:.3f} Mb'))
    axes[-2].set_xlabel(f'GRCh38 genomic position  |  {ann["chr"]}; same x/y axes for both tissues')
    ax=axes[-1]
    for group,role in enumerate(('high','low')):
        for sample in samples[role]:
            pred=read_pred(sample['sample_id'],case['source_model'],tx,args.logits_dir)
            if pred['status']!='OK':
                missing.append(f"{sample['sample_id']}: {pred['status']}")
                continue
            ax.scatter(group-.07,pred['obs'],s=30,alpha=.8,color='tab:blue',
                       label='Observed' if group==0 and sample==samples[role][0] else None)
            ax.scatter(group+.07,pred['pred'],s=30,alpha=.8,marker='s',color='tab:orange',
                       label='Predicted' if group==0 and sample==samples[role][0] else None)
    ax.set_xticks([0,1],[case['tissue_high'],case['tissue_low']])
    ax.set_xlim(-.5,1.5)
    ax.set_ylabel('log2(TPM + 0.001)')
    ax.set_title('Held-out transcript expression, fixed source GCN',fontsize=9)
    ax.legend(loc='upper right',fontsize=7)
    fig.suptitle(f"{case['case_id']}: {case['tissue_high']} vs {case['tissue_low']}  |  {case['rna_type']}  |  fixed GCN: {case['source_model']}\n"
                 'Raw descriptive genomic signal; no causal inference. CpG values require confirmed BED column interpretation.',fontsize=11)
    for role in ('high','low'):
        if source_tags[role]!='EX130_ATTRIBUTION_TARGETS':missing.append(f'{role} tracks use curated sample manifest, not completed attribution target selection')
    outdir=Path(args.outdir)/'figures';outdir.mkdir(parents=True,exist_ok=True)
    file=outdir/f"{case['case_id']}_{tx}"
    fig.savefig(str(file)+'.png',dpi=170);fig.savefig(str(file)+'.svg');plt.close(fig)
    return row_reports,missing


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('mode',choices=('audit','prepare','plot'))
    root=Path.home()
    ap.add_argument('--raw-root',type=Path,default=root/'data_store/IHEC_data_processing_Aug_2023')
    ap.add_argument('--partial-root',type=Path,default=root/'data_store/IHEC_data_processing_Aug_2023_partial_epigenomes')
    ap.add_argument('--annotation-cache-dir',type=Path,default=Path(__file__).with_name('results')/'raw_genomic_tracks'/'shared')
    ap.add_argument('--index-cache-dir',type=Path,default=Path(__file__).with_name('results')/'raw_genomic_tracks'/'shared'/'wgbs_index_cache')
    ap.add_argument('--build-index-cache',action='store_true',help='Opt-in: tabix index on EX130 symlink, never original WGBS BED')
    ap.add_argument('--tabix-bin',default='tabix')
    ap.add_argument('--gff',type=Path,default=root/'IHEC_Mar_2023/networks/g_b/gencode.v29.annotation.gff3')
    ap.add_argument('--logits-dir',type=Path,default=root/'IHEC_Mar_2024/EX51/logits_debug_pc1')
    ap.add_argument('--samples',type=Path,default=Path(__file__).with_name('track_sample_manifest.csv'))
    ap.add_argument('--seed-cases',type=Path,default=Path(__file__).with_name('seed_cases.csv'))
    ap.add_argument('--targets',type=Path,default=Path(__file__).with_name('results')/'fixed_source'/'attribution_targets.csv')
    ap.add_argument('--outdir',type=Path,default=Path(__file__).with_name('results')/'raw_genomic_tracks')
    ap.add_argument('--marks',default='H3K36me3',help='Comma-separated BigWig histone assays, H3K36me3 by default')
    ap.add_argument('--signal',choices=('pval','fc'),default='pval',help='BigWig assay signal; pval matches EX18 feature table')
    ap.add_argument('--select',choices=('auto','targets','curated'),default='auto',help='Select actual attribution targets if present, else curated samples')
    ap.add_argument('--only-cases',default='',help='Comma-separated PILOT case IDs, empty means all')
    ap.add_argument('--max-per-role',type=int,default=3)
    ap.add_argument('--flank',type=int,default=5000)
    ap.add_argument('--bins',type=int,default=750)
    ap.add_argument('--max-cpg',type=int,default=40000)
    ap.add_argument('--display-cpg',type=int,default=10000)
    ap.add_argument('--wgbs-beta-column',default='auto',help='One-based BED data column, or auto if unambiguous')
    ap.add_argument('--wgbs-scale',choices=('fraction','percent'),default='fraction')
    args=ap.parse_args()
    marks=[x.strip() for x in args.marks.split(',') if x.strip()]
    if not marks or any(m not in MARKS for m in marks):ap.error(f'Unsupported marks: {marks}')
    if args.max_per_role<1 or args.bins<2 or args.flank<0 or args.max_cpg<1:ap.error('Invalid numeric arguments')
    samples=csv_rows(args.samples,('sample_id','tissue','rna_type'))
    cases=csv_rows(args.seed_cases,CASE_FIELDS)
    if args.only_cases:
        keep=set(args.only_cases.split(','));cases=[c for c in cases if c['case_id'] in keep]
    if not cases:raise ValueError('No cases selected')
    args.outdir.mkdir(parents=True,exist_ok=True)
    audit=make_sample_audit(samples,args.raw_root,marks,args.signal,args.partial_root,args.index_cache_dir)
    if args.build_index_cache:
        if args.mode!='prepare':raise ValueError('--build-index-cache is supported only in prepare mode')
        index_results=build_cached_indexes(audit,args.index_cache_dir,args.tabix_bin)
        (args.outdir/'index_cache_status.json').write_text(json.dumps(index_results,indent=2)+'\n')
        audit=make_sample_audit(samples,args.raw_root,marks,args.signal,args.partial_root,args.index_cache_dir)
    write_rows(args.outdir/'track_availability.csv',audit,
        list(audit[0]) if audit else ['sample_id','status','issues'])
    problems=[{'sample_id':r['sample_id'],'issues':r['issues']} for r in audit if r['status']!='OK']
    print(f'Track inventory: {len(audit)-len(problems)}/{len(audit)} complete; report: {args.outdir}/track_availability.csv',flush=True)
    if args.mode=='audit':
        if problems:print(json.dumps(problems,indent=2))
        return
    if args.mode=='prepare':
        if not args.gff.is_file():raise FileNotFoundError(args.gff)
        prepared=annotation_cache(args.gff,{c['txid'] for c in cases},args.annotation_cache_dir)
        print(f'Prepared {len(prepared)} of {len(cases)} selected transcript annotations',flush=True)
        if len(prepared)<len(set(c['txid'] for c in cases)):
            raise RuntimeError('GENCODE annotation missing for one or more seed transcripts')
        return
    # Lazy imports, so audit works with no pyBigWig/pysam installed.
    import pyBigWig,pysam,numpy,matplotlib
    if not args.gff.is_file():raise FileNotFoundError(args.gff)
    annotations=annotation_cache(args.gff,{c['txid'] for c in cases},args.annotation_cache_dir)
    targets=csv_rows(args.targets,('candidate_id','txid','role','eval_sample','source_model')) if args.targets.is_file() else []
    audit_map={r['sample_id']:r for r in audit};outcomes=[];track_rows=[]
    for c in cases:
        case_id=c['case_id']
        if c['txid'] not in annotations:
            outcomes.append({'case_id':case_id,'status':'ANNOTATION_NOT_FOUND','issues':'transcript ID absent from GENCODE v29'})
            continue
        try:
            selection,sources=choose_samples(c,samples,targets,args.select,args.max_per_role)
            # Do not silently use samples without requested original epigenomic tracks.
            for role in ('high','low'):
                allowed=[s for s in selection[role] if audit_map[s['sample_id']]['status'] in ('OK','RAW_PRESENT_INDEX_PENDING')]
                if len(allowed)<len(selection[role]):
                    missing=[s['sample_id'] for s in selection[role] if s not in allowed]
                    print(f'{case_id} {role}: skipping samples lacking raw H3K36me3 and/or WGBS tracks: {missing}',file=sys.stderr)
                selection[role]=allowed
                if not allowed:
                    raise ValueError(f'No complete epigenomic tracks for {role} role')
            reports,warnings=plot_case(c,selection,sources,annotations[c['txid']],args,marks)
            track_rows+=reports
            pred={role:[{'sample':s['sample_id'],**read_pred(s['sample_id'],c['source_model'],c['txid'],args.logits_dir)} for s in selection[role]] for role in ('high','low')}
            meta={'case':c,'annotation':annotations[c['txid']],'sample_selection':sources,
                  'samples':{role:[s['sample_id'] for s in selection[role]] for role in ('high','low')},
                  'predictions':pred,'mark_tracks':marks,'signal':args.signal,
                  'wgbs_beta_column':args.wgbs_beta_column,'wgbs_scale':args.wgbs_scale,
                  'warnings':warnings,'generated_utc':datetime.now(timezone.utc).isoformat()}
            (args.outdir/f'{case_id}_provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
            outcomes.append({'case_id':case_id,'status':'PLOTTED_WITH_WARNINGS' if warnings else 'PLOTTED',
                             'issues':';'.join(warnings)})
            print(f'{case_id}: figure written ({len(warnings)} notices)',flush=True)
        except Exception as ex:
            outcomes.append({'case_id':case_id,'status':'FAILED','issues':f'{type(ex).__name__}: {ex}'})
            print(f'{case_id}: FAILED {ex}',file=sys.stderr,flush=True)
    write_rows(args.outdir/'case_plot_status.csv',outcomes,['case_id','status','issues'])
    write_rows(args.outdir/'track_measurement_status.csv',track_rows,
               ['case_id','sample_id','role','tissue','track','status','source','points','annotation'])
    print(f'Case status: {args.outdir}/case_plot_status.csv',flush=True)
    if any(x['status']=='FAILED' for x in outcomes):sys.exit(2)


if __name__=='__main__':main()
