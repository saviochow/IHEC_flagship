#!/usr/bin/env python3
"""Submission harness with a fully mocked sbatch/sinfo; never contacts Slurm."""
import os
import csv
import re
import subprocess
import tempfile
from pathlib import Path
from shutil import copyfile

SRC=Path(__file__).resolve().parent
WORKER_FILES=['ex130_attribution.sbatch','ex130_summarize.sbatch','ex130_tracks_prepare.sbatch','ex130_tracks.sbatch','ex130_tracks_collect.sbatch']
ORIG=(SRC.parent / 'EX130') if (SRC.parent / 'EX130').is_dir() else SRC

with tempfile.TemporaryDirectory() as tmp:
    root=Path(tmp)
    d=root/'project'; d.mkdir()
    (d/'results/fixed_source').mkdir(parents=True)
    (d/'results/fixed_source/attribution_targets.csv').write_text('source_model,eval_sample,txid\nx,x,ENST00000000001\n')
    with (d/'results/fixed_source/attribution_tasks.csv').open('w') as f:
        f.write('task_index,source_model,eval_sample,n_target_rows,n_unique_transcripts\n')
        for i in range(69): f.write(f'{i},x,x,1,1\n')
    with (d/'seed_cases.csv').open('w') as f:
        f.write('case_id,txid\n')
        for i in range(13): f.write(f'PILOT_{i:02},ENST00000000001\n')
    (d/'ex130_raw_tracks.py').write_text('#stub\n')
    (d/'transcript_attribution.py').write_text('#patched for python3.8\n')
    for fn in WORKER_FILES: copyfile(ORIG/fn,d/fn)
    bind=root/'bin';bind.mkdir(); log=root/'sbatch.log'; counter=root/'counter'; counter.write_text('200000')
    (bind/'sinfo').write_text('#!/bin/sh\ncase "$*" in *high*) printf "n003\\n";; *normal*) printf "n001\\nn002\\n";; esac\n')
    (bind/'sbatch').write_text('''#!/usr/bin/env bash
set -euo pipefail
printf '%s\\n' "$*" >> "$MOCK_SBATCH_LOG"
cur=$(cat "$MOCK_SBATCH_COUNTER")
next=$((cur+1))
printf '%s' "$next" > "$MOCK_SBATCH_COUNTER"
printf '%s\\n' "$next"
''')
    for p in bind.iterdir():p.chmod(0o755)
    env=os.environ.copy();env.update(EX130_DIR=str(d),PYTHON_BIN='/usr/bin/python3',TRACK_PYTHON_BIN='/usr/bin/python3',MOCK_SBATCH_LOG=str(log),MOCK_SBATCH_COUNTER=str(counter),PATH=str(bind)+':'+env['PATH'])
    def call(*args, expect=0, extra=None):
        e=dict(env, **(extra or {}))
        proc=subprocess.run(['bash',str(SRC/'submit_ex130_balanced.sh'),*args],env=e,capture_output=True,text=True)
        assert proc.returncode==expect,(proc.returncode,proc.stdout,proc.stderr)
        return proc

    # Dry run never attempts to schedule a job.
    p=call('all','--dry-run')
    assert 'ATTR_DISTRIBUTION total=69 high:n003=30 normal:n001=20 normal:n002=19 in_flight=2,1,1' in p.stdout,p.stdout
    assert 'TRACK_DISTRIBUTION total=13 high:n003=6 normal:n001=4 normal:n002=3 in_flight=2,1,1' in p.stdout,p.stdout
    assert not log.exists()

    # Actual submission only hits synthetic sbatch. All three arrays and their dependencies.
    p=call('all','--submit')
    lines=log.read_text().splitlines()
    assert len(lines)==9,len(lines)
    expect_nodes=['--partition=high --nodelist=n003','--partition=normal --nodelist=n001','--partition=normal --nodelist=n002']
    for i in range(3):
        assert expect_nodes[i] in lines[i],lines[i]
        assert '--array=' in lines[i],lines[i]
    assert '--dependency=afterok:200001:200002:200003' in lines[3],lines[3]
    assert '--partition=normal --nodelist=n001' in lines[4],lines[4]
    for i in range(3):
        assert expect_nodes[i] in lines[5+i],lines[5+i]
        assert '--dependency=afterok:200005' in lines[5+i],lines[5+i]
    assert '--dependency=afterany:200006:200007:200008' in lines[8],lines[8]
    assert len(list((d/'results/balanced_submissions').glob('*_all_jobs.txt')))==1
    print('PASS weighted arrays: attribution 30/20/19, tracks 6/4/3')
    print('PASS default throttles: 2/1/1 per workflow')
    print('PASS two independent three-way dependency fan-ins')
    print('PASS job partition/node pinning and receipts')

    # Test adjustable concurrency and negative controls.
    q=call('attribution','--dry-run',extra={'MAX_CONCURRENT':'8'})
    assert 'in_flight=4,2,2' in q.stdout,q.stdout
    q=call('attribution','--submit',extra={'MAX_CONCURRENT':'2'},expect=2)
    assert '>=3' in q.stderr
    (d/'transcript_attribution.py').write_text('ap.add_argument("--self-loop", action=argparse.BooleanOptionalAction, default=True)\n')
    q=call('attribution','--submit',expect=3)
    assert 'Python 3.8 incompatibility' in q.stderr
    print('PASS Py3.8 incompatibility and malformed resource limits refused before submission')

    # Idempotent patch helper can safely update the exact known source pattern.
    copyfile(SRC/'patch_ex130_python38.py', d/'patch_ex130_python38.py')
    proc=subprocess.run(['/usr/bin/python3',str(d/'patch_ex130_python38.py')],capture_output=True,text=True)
    assert proc.returncode==0,(proc.stdout,proc.stderr)
    assert 'ALREADY_PRESENT' in subprocess.run(['/usr/bin/python3',str(d/'patch_ex130_python38.py')],capture_output=True,text=True).stdout
    assert 'BooleanOptionalAction' not in (d/'transcript_attribution.py').read_text()
    print('PASS idempotent narrow Python 3.8 patch')
