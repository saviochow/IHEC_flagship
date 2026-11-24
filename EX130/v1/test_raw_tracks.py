"""Synthetic end-to-end sanity test for EX130 raw track reader and figures.

Creates fake BigWig and Tabix backends. Does not validate the true Yipcluster BED
schema or the GENCODE v29 source; live audit and plot remain required.
"""
import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from ex130_raw_tracks import find_annotations,infer_beta_column,locate_tracks,build_cached_indexes,annotation_cache,choose_samples


def test_parsers():
    assert infer_beta_column([(15,['chr19','15','16','0.7'])],'auto','fraction')==(3,'AUTO_FOUR_COLUMN_BED')
    assert infer_beta_column([(15,['chr19','15','16','70'])],'auto','fraction')[0] is None
    assert infer_beta_column([(15,['chr19','15','16','70'])],'auto','percent')[0]==3
    assert infer_beta_column([(15,['chr19','15','16','0.7','0.3'])],'auto','fraction')[0] is None
    with tempfile.TemporaryDirectory() as d:
        p=Path(d)/'a.gff3'
        p.write_text('##gff-version 3\nchr19\tsrc\ttranscript\t101\t301\t.\t-\t.\tID=ENST00000000001.1;gene_name=GENE%2FNAME;gene_id=ENSG01\n'
                     'chr19\tsrc\texon\t101\t151\t.\t-\t.\tParent=ENST00000000001.1\n'
                     'chr19\tsrc\texon\t251\t301\t.\t-\t.\tParent=ENST00000000001.1\n')
        records=find_annotations(p,{'ENST00000000001.1'})
        x=records['ENST00000000001.1'];assert x['start0']==100 and x['end0']==301
        assert x['exons']==[(100,151),(250,301)],x
        assert x['gene']=='GENE/NAME'
        base=Path(d)/'raw'/'IHECRE00000001'
        (base/'ChIP-Seq').mkdir(parents=True);(base/'WGBS').mkdir()
        (base/'ChIP-Seq'/'IHECRE00000001.4.H3K36me3.pval.signal.bigwig').write_bytes(b'fake')
        (base/'WGBS'/'IHECRE00000001.4.PBAT-WGBS.cpg.all.bed.gz').write_bytes(b'fake')
        (base/'WGBS'/'IHECRE00000001.4.PBAT-WGBS.cpg.all.bed.gz.tbi').write_bytes(b'fake')
        found=locate_tracks('IHECRE00000001',Path(d)/'raw',['H3K36me3'])
        assert found['status']=='OK' and found['wgbs_protocol']=='PBAT-WGBS',found
        (base/'WGBS'/'IHECRE00000001.4.PBAT-WGBS.cpg.all.bed.gz.tbi').unlink()
        found=locate_tracks('IHECRE00000001',Path(d)/'raw',['H3K36me3'])
        assert 'WHOLE_REGION_CPG_INDEX_MISSING' in found['issues']


def test_full():
    with tempfile.TemporaryDirectory() as d:
        d=Path(d);mocks=d/'mocks';mocks.mkdir()
        (mocks/'pyBigWig.py').write_text('''
class Big:
 def __init__(self,p):self.p=str(p)
 def __enter__(self):return self
 def __exit__(self,*args):return False
 def close(self):pass
 def chroms(self):return {'chr19':100000}
 def stats(self,chrom,start,end,nBins,type):return [2.0 if '00001' in self.p else .2]*nBins
def open(p):return Big(p)
''')
        (mocks/'pysam.py').write_text('''
class TabixFile:
 def __init__(self,p):self.p=p;self.contigs=['chr19']
 def __enter__(self):return self
 def __exit__(self,*args):return False
 def fetch(self,chrom,start,end):return ['chr19\\t110\\t111\\t0.2','chr19\\t200\\t201\\t0.9','chr19\\t260\\t261\\t0.4']
''')
        raw=d/'raw'
        sample_ids=['IHECRE00000001','IHECRE00000002']
        for sample in sample_ids:
            base=raw/sample
            (base/'ChIP-Seq').mkdir(parents=True);(base/'WGBS').mkdir()
            (base/'ChIP-Seq'/f'{sample}.1.H3K36me3.pval.signal.bigwig').write_bytes(b'fake')
            (base/'WGBS'/f'{sample}.1.PBAT-WGBS.cpg.all.bed.gz').write_bytes(b'fake')
            (base/'WGBS'/f'{sample}.1.PBAT-WGBS.cpg.all.bed.gz.tbi').write_bytes(b'fake')
        annotation=d/'a.gff3'
        annotation.write_text('##gff-version 3\nchr19\tsrc\ttranscript\t101\t301\t.\t+\t.\tID=ENST00000000001.1;gene_name=TESTGENE;gene_id=ENSG01\n'
                              'chr19\tsrc\texon\t101\t151\t.\t+\t.\tParent=ENST00000000001.1\n'
                              'chr19\tsrc\texon\t251\t301\t.\t+\t.\tParent=ENST00000000001.1\n')
        cases=d/'seed.csv'
        cases.write_text('case_id,txid,rna_type,tissue_high,tissue_low,source_model\nPILOT_01,ENST00000000001.1,mRNA,Liver,Colon,IHECRE00000001\n')
        samples=d/'samples.csv'
        samples.write_text('sample_id,tissue,rna_type\nIHECRE00000001,Liver,mRNA\nIHECRE00000002,Colon,mRNA\n')
        l=d/'logits';l.mkdir()
        for sample in sample_ids:
            (l/f'10sample_+thyroid_IHECRE00000001_{sample}_128_test_logits.csv').write_text(
               'test_id,logits,labels,txid,chr\n100,1.0,0.8,ENST00000000001.1,chr19\n')
        output=d/'output';env=dict(os.environ,PYTHONPATH=str(mocks)+os.pathsep+str(ROOT))
        for mode in ('audit','plot'):
            subprocess.run([sys.executable,str(ROOT/'ex130_raw_tracks.py'),mode,
                            '--raw-root',str(raw),'--samples',str(samples),'--seed-cases',str(cases),
                            '--gff',str(annotation),'--logits-dir',str(l),'--outdir',str(output),
                            '--flank','20','--bins','30'],env=env,check=True)
        assert (output/'figures/PILOT_01_ENST00000000001.1.png').stat().st_size>3000
        assert (output/'figures/PILOT_01_ENST00000000001.1.svg').stat().st_size>3000
        assert (output/'case_plot_status.csv').read_text().splitlines()[1].startswith('PILOT_01,PLOTTED_WITH_WARNINGS')
        assert 'AUTO_FOUR_COLUMN_BED_COL_4' in (output/'track_measurement_status.csv').read_text()
        assert 'PREDICTION_CSV' not in (output/'PILOT_01_provenance.json').read_text()
        print('PASS raw-track synthetic integration: indexed PBAT, GENCODE coordinates, BigWig, CpG, files and visual exports')



def test_partial_roots_and_index_cache():
    """Original IHECRE00001892 should use partial raw tree; 03822 is not missing data."""
    with tempfile.TemporaryDirectory() as d:
        d=Path(d);main=d/'main';partial=d/'partial';index=d/'indexes'
        for sid in ('IHECRE00001892','IHECRE00003822'):
            root=partial if sid=='IHECRE00001892' else main
            path=root/sid
            (path/'ChIP-Seq').mkdir(parents=True)
            (path/'WGBS').mkdir()
            (path/'ChIP-Seq'/f'{sid}.1.H3K36me3.pval.signal.bigwig').write_bytes(b'bigwig')
            (path/'WGBS'/f'{sid}.1.standard-WGBS.cpg.all.bed.gz').write_bytes(b'raw-data')
        (main/'IHECRE00001892').mkdir(parents=True)
        (main/'IHECRE00001892'/'RNA-Seq').mkdir()
        p0=partial/'IHECRE00001892'/'WGBS'/'IHECRE00001892.1.standard-WGBS.cpg.all.bed.gz'
        Path(str(p0)+'.tbi').write_bytes(b'index')
        a=locate_tracks('IHECRE00001892_mRNA-Seq',main,['H3K36me3'],secondary_root=partial,index_cache_dir=index)
        assert a['status']=='OK',a
        assert a['raw_dir']==str(partial/'IHECRE00001892')
        b=locate_tracks('IHECRE00003822',main,['H3K36me3'],secondary_root=partial,index_cache_dir=index)
        assert b['status']=='RAW_PRESENT_INDEX_PENDING',b
        assert b['mark_paths']['H3K36me3'] and b['wgbs_path']
        # Stub tabix only writes to the EX130 output symlink, never raw data.
        stub=d/'tabix'
        stub.write_text('#!/bin/sh\ncase "$1" in -p) ;; *) exit 11;; esac\n[ "$2" = bed ] || exit 12\nprintf stub > "$3.tbi"\n')
        stub.chmod(0o700)
        results=build_cached_indexes([b],index,str(stub))
        assert results[0]['status']=='CACHED_INDEX_CREATED',results
        b2=locate_tracks('IHECRE00003822',main,['H3K36me3'],secondary_root=partial,index_cache_dir=index)
        assert b2['status']=='OK',b2
        assert b2['wgbs_path']!=b['wgbs_path']
        assert Path(b['wgbs_path']).read_bytes()==b'raw-data'
        print('PASS partial epigenomes and unindexed raw-bed index caching')


def test_shared_annotation_cache():
    with tempfile.TemporaryDirectory() as d:
        d=Path(d);src=d/'a.gff3';cache=d/'shared';cache.mkdir()
        src.write_text('##gff-version 3\nchr18\tsrc\ttranscript\t101\t130\t.\t+\t.\tID=ENST00000000001.1\n'
                       'chr19\tsrc\ttranscript\t201\t230\t.\t-\t.\tID=ENST00000000002.1\n')
        bulk=annotation_cache(src,{'ENST00000000001.1','ENST00000000002.1'},cache)
        subset=annotation_cache(src,{'ENST00000000002.1'},cache)
        assert len(bulk)==2 and list(subset)==['ENST00000000002.1']
        print('PASS shared multi-transcript GENCODE cache reused by one-case workers')


def test_pilot_sample_manifest():
    import csv
    root=ROOT
    with (root/'pilot_candidates_with_samples.csv').open(newline='') as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==17
    for r in rows:
        hi=r['high_sample_ids'].split(';');lo=r['low_sample_ids'].split(';')
        assert len(hi)==int(r['n_high_samples'])==3
        assert len(lo)==int(r['n_low_samples'])==3
        if r['screen_type']=='diagonal':
            assert r['high_trained_sample_ids']==r['high_sample_ids']
            assert r['low_trained_sample_ids']==r['low_sample_ids']
        else:
            assert r['high_trained_sample_ids']==r['trained_sample']==r['low_trained_sample_ids']
    venous=[r for r in rows if r['lower_tissue']=='Venous Blood']
    assert venous and all('IHECRE00001892_mRNA-Seq' in r['low_sample_ids'] for r in venous)
    cord=[r for r in rows if r['higher_tissue']=='Cord Blood']
    assert cord and all('IHECRE00003822' in r['high_sample_ids'] for r in cord)
    print('PASS every pilot row has exact three sample identities per tissue')

if __name__=='__main__':
    test_parsers();test_full();test_partial_roots_and_index_cache();test_shared_annotation_cache();test_pilot_sample_manifest()
