import math
from pathlib import Path
import pandas as pd
import subprocess
import sys
import tempfile

with tempfile.TemporaryDirectory() as d:
    d = Path(d)
    base = d/'IHEC_Mar_2024'
    out = d/'out'
    feat = d/'features'
    (base/'EX18'/'samples_needed').mkdir(parents=True)
    (base/'EX51'/'logits_debug_pc1').mkdir(parents=True)
    (base/'metadata').mkdir(parents=True)
    lines=['Sample ID,Tissue origin (simplified),Biomaterial type,Total RNA−Seq availability,mRNA−Seq availability,WGBS protocol,Sex,Disease status,Sample pooling status,project']
    for i in range(6):
        sid=f'IHECRE{i+1:08d}'
        tissue='Liver' if i<3 else 'Breast'
        lines.append(f'{sid},{tissue},tissue,F,T,WGBS,M,non-disease,single,CEEHRC')
        (base/'EX18'/'samples_needed'/sid).mkdir()
        yhi=math.log2(10+1e-3); ylo=math.log2(.3+1e-3)
        df=pd.DataFrame([{'test_id':1,'logits':yhi-.12 if i<3 else ylo+.1,'labels':yhi if i<3 else ylo,'txid':'ENST00000012345.2','chr':'chr18'},{'test_id':2,'logits':0,'labels':0,'txid':'ENST00000000555.2','chr':'chr19'}, {'test_id':3,'logits':9,'labels':9,'txid':'ENSTtrain.1','chr':'chr1'}])
        df.to_csv(base/'EX51'/'logits_debug_pc1'/f'10sample_+thyroid_{sid}_{sid}_128_test_logits.csv',index=False)
        fd=feat/sid;fd.mkdir(parents=True)
        h=pd.DataFrame([{'transcript.id':'ENST00000012345.2','gene.id':'ENSG0001.1','symbol':'GENEX','hg38H3K36me3intnintrons':10 if i<3 else .5,'hg38H3K36me3gene_body':10,'hg38H3K36me3firstintron':2,'hg38H3K36me3lastintron':5}])
        h.to_csv(fd/f'{sid}.4.H3K36me3.pval.all.txt',sep='\t',index=False)
        m=pd.DataFrame([{'transcript.id':'ENST00000012345.2','gene.id':'ENSG0001.1','symbol':'GENEX','beta_hg38intnintrons':.8 if i<3 else .2,'beta_hg38gene_body':.5,'beta_hg38firstintron':.3,'beta_hg38promoter':.1,'beta_hg38TSS1':.05}])
        m.to_csv(fd/f'{sid}.4.standard-WGBS.cpg.all.txt',sep='\t',index=False)
    (base/'metadata'/'SuppTable1_with_project.tsv').write_text('\n'.join(lines)+'\n')
    script=Path(__file__).parent/'screen_ihec_cases.py'
    cmd=[sys.executable, str(script), 'screen','--base',str(base),'--outdir',str(out),'--min-group-n','3','--top-per-type','4']
    subprocess.run(cmd,check=True)
    pairs=pd.read_csv(out/'candidate_pairs.csv')
    assert len(pairs) == 1 and pairs.iloc[0].txid=='ENST00000012345.2',pairs
    reps=pd.read_csv(out/'representative_samples.csv')
    assert len(reps)==4,reps
    subprocess.run([sys.executable,str(script),'features','--outdir',str(out),'--features-root',str(feat)],check=True)
    f=pd.read_csv(out/'candidate_ranked_with_features.csv')
    assert f.iloc[0].feature_support_class=='both',f
    assert f.iloc[0].gene_symbol=='GENEX',f
    print('PASS: synthetic end-to-end test, held-out chromosome, replicates, provenance, feature contrasts and gene annotation')
