#!/usr/bin/env python3
"""Test EX130 Slurm orchestration with a fake sbatch, never submit a job."""
import csv,os,subprocess,sys,tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parent

def run():
    with tempfile.TemporaryDirectory() as d:
        root=Path(d);bin_dir=root/'bin';bin_dir.mkdir();log=root/'slurm.log';counter=root/'count'
        mock=bin_dir/'sbatch';mock.write_text('''#!/bin/sh
printf '%s\\n' "$*" >> "$EX130_TEST_SUBMIT_LOG"
n=7001
[ -e "$EX130_TEST_COUNT" ] && n=$(($(cat "$EX130_TEST_COUNT")+1))
printf '%s' "$n" > "$EX130_TEST_COUNT"
printf '%s\\n' "$n"
''');mock.chmod(0o700)
        env=dict(os.environ,PATH=str(bin_dir)+os.pathsep+os.environ.get('PATH',''),
            EX130_DIR=str(ROOT),TRACK_PYTHON_BIN=sys.executable,
            EX130_TEST_SUBMIT_LOG=str(log),EX130_TEST_COUNT=str(counter),MAX_CONCURRENT='4')
        p=subprocess.run(['bash',str(ROOT/'submit_ex130_tracks.sh')],env=env,text=True,
                         capture_output=True,check=True)
        calls=log.read_text().splitlines();assert len(calls)==3,calls
        assert '--dependency=afterok:7001' in calls[1],calls[1]
        assert '--array=0-12%4' in calls[1],calls[1]
        assert '--dependency=afterany:7002' in calls[2],calls[2]
        assert 'CASE_TASKS=13' in p.stdout,p.stdout
        print('PASS three-job preparation, 13-case array throttle %4, afterany collector')
        # Create fake isolated case receipts for collector, then test both outcomes.
        for case in csv.DictReader((ROOT/'seed_cases.csv').open()):
            folder=root/'results/raw_genomic_tracks/cases'/case['case_id']
            folder.mkdir(parents=True)
            (folder/'case_plot_status.csv').write_text('case_id,status,issues\n'+case['case_id']+',PLOTTED,\n')
            (folder/'track_measurement_status.csv').write_text('case_id,sample_id,role,tissue,track,status,source,points,annotation\n')
        (root/'seed_cases.csv').write_text((ROOT/'seed_cases.csv').read_text())
        job=[sys.executable,str(ROOT/'collect_ex130_tracks.py'),'--ex130-dir',str(root)]
        p=subprocess.run(job,text=True,capture_output=True)
        assert p.returncode==0,p.stderr
        folder=root/'results/raw_genomic_tracks/cases/PILOT_13'
        (folder/'case_plot_status.csv').unlink()
        p=subprocess.run(job,text=True,capture_output=True)
        assert p.returncode==2,p.stderr
        print('PASS isolated array receipt collection and missing-case failure detection')

if __name__=='__main__':run()
