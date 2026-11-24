"""Synthetic end-to-end EX130 test (fake DGL adapter; NOT live Yipcluster verification)."""
from __future__ import annotations

import math
import os
import pickle
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from feature_schema import ASSAYS, REGIONS, mapping
from transcript_attribution import target_node_gradient, verify_wgbs_first

ROOT = Path(__file__).resolve().parent


def call(args, env):
    subprocess.run([sys.executable, *map(str, args)], check=True, env=env)


def test_scalar_gradient():
    x = torch.tensor([[1., 2.], [3., 4.], [5., 6.]], requires_grad=True)
    A = torch.tensor([[1., 2., 0.], [0., 1., 0.], [2., 0., 3.]])
    class Model:
        def __call__(self, z):
            return A @ z @ torch.tensor([[1.], [2.]])
    y, grad = target_node_gradient(Model(), x, 0)
    assert y.item() == 27
    assert torch.allclose(grad, torch.tensor([[1., 2.], [2., 4.], [0., 0.]]))
    # In particular, node 2 does not affect node 0 in this graph.


def run():
    test_scalar_gradient()
    with tempfile.TemporaryDirectory() as temp:
        temp = Path(temp)
        mocks = temp / "mocks" / "dgl"
        (mocks / "data").mkdir(parents=True)
        (mocks / "__init__.py").write_text('''
import torch
class FakeGraph:
 def __init__(self, nd): self.ndata = nd
 def int(self): return self
 def in_degrees(self): return torch.ones(len(self.ndata["feat"]), dtype=torch.int32)
 def num_nodes(self): return len(self.ndata["feat"])
 def num_edges(self): return len(self.ndata["feat"])
def remove_self_loop(g): return g
def add_self_loop(g): return g
__version__="FAKE_DGL_FOR_TESTS_ONLY"
''')
        (mocks / "data" / "__init__.py").write_text('''
import pickle
from dgl import FakeGraph
def load_graphs(path):
 with open(path, "rb") as fh: d = pickle.load(fh)
 return [FakeGraph(d)], {}
''')
        base = temp / "IHEC_Mar_2024"
        for p in (base / "EX18" / "samples_needed", base / "EX18" / "graphs_tpmpc1e-3_hmpc1e-9_zscale_clip5_WGBS+HM" / "10sample_+thyroid",
                  base / "EX51" / "models", base / "EX51" / "logits_debug_pc1", base / "metadata"):
            p.mkdir(parents=True, exist_ok=True)
        (base / "EX51" / "gcn_regression.py").write_text('''
import torch
from torch import nn
class GCN(nn.Module):
 def __init__(self,g,in_feats,n_hidden,n_classes,n_layers,activation,dropout):
  super().__init__()
  self.l = nn.Linear(in_feats,1)
 def forward(self,x):return self.l(x)
''')
        # Create a deterministic checkpoint compatible with the synthetic GCN.
        sys.path.insert(0, str(base / "EX51"))
        from gcn_regression import GCN
        m = GCN(None, 147, 128, 1, 7, None, .45)
        with torch.no_grad():
            m.l.weight.fill_(0)
            m.l.weight[0, 0] = 6
            m.l.weight[0, 63] = 1
            m.l.bias.fill_(-2.5)
        sample_ids = [f"IHECRE{i+1:08d}" for i in range(6)]
        torch.save(m.state_dict(), base / "EX51" / "models" / f"{sample_ids[0]}_128_neurons.model")
        txid = "ENST00000012345.2"
        features_root = temp / "raw_features"
        metadata = ['Sample ID,Tissue origin (simplified),Biomaterial type,Total RNA−Seq availability,mRNA−Seq availability,WGBS protocol,Sex,Disease status,Sample pooling status,project']
        for i, sid in enumerate(sample_ids):
            high = i < 3
            tpm = 10 if high else .3
            tissue = "Liver" if high else "Breast"
            metadata.append(f"{sid},{tissue},tissue,F,T,WGBS,F,non-disease,single,CEEHRC")
            (base / "EX18" / "samples_needed" / sid).mkdir()
            x = torch.zeros(2, 147)
            x[0, :21] = .8 if high else .2
            x[0, 63] = 1 if high else 0
            with torch.no_grad():
                pred = float(m(x)[0, 0])
            label = math.log2(tpm + 0.001)
            data = {'feat': x, 'exp_label': torch.tensor([label, -1.0]),
                    'test_id': torch.tensor([1, -1]), 'test_mask': torch.tensor([1, 0]),
                    'train_mask': torch.tensor([0, 0]), 'val_mask': torch.tensor([0, 0])}
            with (base / "EX18" / "graphs_tpmpc1e-3_hmpc1e-9_zscale_clip5_WGBS+HM" / "10sample_+thyroid" / f"{sid}.all.dgl").open("wb") as fd:
                pickle.dump(data, fd)
            pd.DataFrame([
                {'test_id': 1, 'logits': pred, 'labels': label, 'txid': txid, 'chr': 'chr18'},
                {'test_id': 2, 'logits': -9.9, 'labels': -9.9, 'txid': 'ENSTtrain.1', 'chr': 'chr1'},
            ]).to_csv(base / "EX51" / "logits_debug_pc1" / f"10sample_+thyroid_{sample_ids[0]}_{sid}_128_test_logits.csv", index=False)
            fd = features_root / sid
            fd.mkdir(parents=True)
            raw = {'transcript.id': txid, 'gene.id': 'ENSG00000012345.1', 'symbol': 'TESTCASE'}
            raw.update({f'beta_hg38{region}': (.8 if high else .2) for region in REGIONS})
            raw['beta_hg38gene_body'] = raw['beta_hg38intnintrons']
            raw['beta_hg38promoter'] = .1
            pd.DataFrame([raw]).to_csv(fd / f'{sid}.4.standard-WGBS.cpg.all.txt', index=False, sep='\t')
            hraw = {'transcript.id': txid, 'gene.id': 'ENSG00000012345.1', 'symbol': 'TESTCASE'}
            hraw.update({f'hg38H3K36me3{region}': (10 if high else .1) for region in REGIONS})
            hraw['hg38H3K36me3gene_body'] = 10 if high else .1
            pd.DataFrame([hraw]).to_csv(fd / f'{sid}.4.H3K36me3.pval.all.txt', index=False, sep='\t')
        (base / 'metadata' / 'SuppTable1_with_project.tsv').write_text('\n'.join(metadata) + '\n')
        out = temp / "run"
        env = {**os.environ, 'PYTHONPATH': str(mocks.parent) + os.pathsep + str(ROOT)}
        call([ROOT / 'screen_ihec_cases.py', 'screen', '--base', base,
              '--trained-sample', sample_ids[0], '--outdir', out,
              '--min-group-n', '3', '--top-per-type', '4'], env)
        call([ROOT / 'screen_ihec_cases.py', 'features', '--features-root', features_root,
              '--outdir', out], env)
        call([ROOT / 'make_attribution_manifest.py', '--screen-dir', out,
              '--logits-dir', base / 'EX51' / 'logits_debug_pc1',
              '--outdir', out, '--seed-cases', temp / 'no_seeds.csv'], env)
        call([ROOT / 'export_attribution_features.py', '--outdir', out,
              '--features-root', features_root], env)
        exported = pd.read_csv(out / 'attribution_target_raw_features.csv')
        assert exported.has_both_internal_introns.all(), exported
        tasks = pd.read_csv(out / 'attribution_tasks.csv')
        assert len(tasks) == 4, tasks  # 2 cases sides x 2 reps each
        for i in range(len(tasks)):
            call([ROOT / 'transcript_attribution.py', '--base', base,
                  '--features-root', features_root, '--outdir', out, '--task-index', str(i)], env)
        call([ROOT / 'summarize_attribution.py', '--outdir', out], env)
        call([ROOT / 'plot_candidate_panels.py', '--outdir', out, '--max-cases', '5'], env)
        assert list((out / 'candidate_panel_png').glob('*.png'))
        merged = pd.read_csv(out / 'candidate_per_sample_assay_attribution.csv')
        assert merged.assay.nunique() == 7
        dna = merged.query("assay=='DNAm'")
        assert np.all(dna.self_gradient_abs_share > 0.8), dna
        h3 = merged.query("assay=='H3K36me3'")
        assert np.all(h3.self_gradient_abs_share > 0.05), h3
        # Negative audit: swapping one WGBS input must be detected, not relabelled.
        bad = torch.zeros(2, 147)
        bad[0, :21] = 0.3
        try:
            verify_wgbs_first(features_root, sample_ids[0], {txid: 0}, bad)
        except ValueError as exc:
            assert 'FEATURE ORDER MISMATCH' in str(exc)
        else:
            raise AssertionError('Mismatched WGBS tensor was not rejected')
        print('PASS synthetic EX130 full workflow: fixed-model discovery, raw feature joins, manifests, scalar target gradient, DGL tensor/WGBS order audit, provenance, merge and mismatch rejection')


if __name__ == '__main__':
    run()
