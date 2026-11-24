"""Managed adapter tests, using no live cluster or human research data."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import managed_payload as payload
import managed_artifacts as artifacts

class ManagedAdapterTest(unittest.TestCase):
    def config(self):
        return {"cluster":"yipcluster","ex130_root":"/tmp/ex130","logits_dir":"/tmp/logits",
                "gff":"/tmp/annotation.gff3","features_root":"/tmp/features",
                "python":"/usr/bin/python3","source_model":"FIXED","workers":4,
                "manifest_sha256":payload.MANIFEST_SHA}

    def test_original_bytes_and_bash_syntax(self):
        script=payload.render(ROOT,self.config())
        self.assertIn(b'RESEARCHSTACK_TASK_DIR/ex130-v2',script)
        self.assertIn(b'V2_SOURCE_MODEL=FIXED',script)
        self.assertNotIn(b'\n#SBATCH',script)
        proc=subprocess.run(["bash","-n"],input=script,capture_output=True)
        self.assertEqual(proc.returncode,0,proc.stderr)
        # All private data path and worker settings are explicit; no caller
        # variables are needed beyond the canonical runtime envelope.
        self.assertEqual(script,payload.render(ROOT,self.config()))

    def test_invalid_config_and_modified_source(self):
        for change in ({"features_root":"/tmp/../other"},{"workers":9},
                       {"source_model":"$(echo unsafe)"},{"manifest_sha256":"0"*64},
                       {"cluster":"unknown"}):
            with self.assertRaises(ValueError):
                payload.render(ROOT,{**self.config(),**change})
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            (p/"isoform_switch.py").write_text("changed\n")
            with self.assertRaises(ValueError):
                payload.render(p,self.config())

    def test_artifact_completeness_and_private_data_exclusion(self):
        with tempfile.TemporaryDirectory() as d:
            task=Path(d).resolve()
            a=task/"ex130-v2"/"analysis"
            case=a/"case_runs"/"000"
            case.mkdir(parents=True)
            (a/"discovery_receipt.json").write_text('{"status":"DISCOVERY_COMPLETE"}')
            (a/"collection_receipt.json").write_text(
                '{"status":"COMPLETE","expected_cases":1,"present_receipts":1}')
            (a/"cases.csv").write_text("case_id\nSWITCH_001\n")
            (a/"candidate_validation.csv").write_text("case_id\nSWITCH_001\n")
            (case/"receipt.json").write_text(
                '{"index":0,"case_id":"SWITCH_001","status":"VALIDATED"}')
            (case/"isoform_pair.png").write_bytes(b"synthetic-figure")
            (a/"expression_cache.sqlite").write_bytes(b"private")
            (case/"paired_features.csv").write_text("private\n")
            self.assertEqual(artifacts.collect(task,0)["status"],"INCOMPLETE")
            (case/"isoform_pair.svg").write_text("<svg/>")
            r=artifacts.collect(task,0)
            self.assertEqual(r["status"],"COMPLETE")
            paths=[x["path"] for x in r["files"]]
            self.assertFalse(any("sqlite" in x or "paired_features" in x for x in paths))
            self.assertEqual(artifacts.collect(task,2)["status"],"INCOMPLETE")
            f=next(x for x in r["files"] if x["kind"]=="figure_svg")
            self.assertEqual(f["sha256"],hashlib.sha256(b"<svg/>").hexdigest())

    def test_managed_payload_full_synthetic_plot_execution(self):
        # Mock the runtime envelope ONLY in disposable synthetic CI. This is
        # not a Slurm canary, admission receipt or authorized live allocation.
        import os
        import shutil
        from test_isoform_switch import create_fixture
        with tempfile.TemporaryDirectory() as d:
            root=Path(d).resolve()
            manifest,logits,gff,features=create_fixture(root)
            ex130=root/"ex130"
            target=ex130/"results/fixed_source/sample_manifest.csv"
            target.parent.mkdir(parents=True)
            shutil.copyfile(manifest,target)
            run=root/"managed-run"
            task=run/"task-0"
            task.mkdir(parents=True)
            context=run/"context.json"
            context.write_text(json.dumps({"cluster":"yipcluster","run_dir":str(run),
                "nodes":1,"trials":1}))
            context_sha=hashlib.sha256(context.read_bytes()).hexdigest()
            manifest_sha=hashlib.sha256(target.read_bytes()).hexdigest()
            cfg={**self.config(),"ex130_root":str(ex130),"logits_dir":str(logits),
                 "gff":str(gff),"features_root":str(features),
                 "python":sys.executable,"workers":2,"manifest_sha256":manifest_sha}
            # The production digest remains pinned. This patch is confined to
            # the test renderer and binds only the synthetic fixture bytes.
            with patch.object(payload,"MANIFEST_SHA",manifest_sha):
                script=payload.render(ROOT,cfg)
            home=root/"home"
            home.mkdir()
            env={"PATH":"/usr/bin:/bin","HOME":str(home),"LANG":"C.UTF-8",
                 "MPLCONFIGDIR":str(home/"mpl"),"RESEARCHSTACK_CONTEXT":str(context),
                 "RESEARCHSTACK_CONTEXT_SHA256":context_sha,
                 "RESEARCHSTACK_RUN_DIR":str(run),"RESEARCHSTACK_TASK_DIR":str(task),
                 "SLURM_JOB_ID":"1","SLURM_CPUS_PER_TASK":"2"}
            result=subprocess.run(["/bin/bash","-s","--"],input=script,cwd=task,
                env=env,capture_output=True,timeout=180)
            self.assertEqual(result.returncode,0,
                (result.stdout+result.stderr).decode(errors="replace")[-8000:])
            receipt=json.loads((task/"ex130_artifact_manifest.json").read_text())
            self.assertEqual(receipt["status"],"COMPLETE")
            self.assertEqual(receipt["expected_cases"],2)
            for kind,magic in (("figure_png",b"\x89PNG\r\n\x1a\n"),("figure_svg",b"<?xml")):
                figures=[f for f in receipt["files"] if f["kind"]==kind]
                self.assertEqual(len(figures),2)
                for f in figures:
                    data=(task/f["path"]).read_bytes()
                    self.assertTrue(data.startswith(magic))
                    self.assertGreater(len(data),1000)
                    self.assertEqual(hashlib.sha256(data).hexdigest(),f["sha256"])
                    self.assertEqual(len(data),f["bytes"])
            summary=json.loads((task/"ex130-v2/analysis/collection_receipt.json").read_text())
            self.assertEqual(summary["positive_expression_and_both_mark_flips"],1)
            self.assertEqual(summary["stable_expression_and_no_mark_change"],1)
            repeat=subprocess.run(["/bin/bash","-s","--"],input=script,cwd=task,
                env=env,capture_output=True,timeout=30)
            self.assertNotEqual(repeat.returncode,0)
            self.assertEqual(json.loads((task/"ex130_artifact_manifest.json").read_text()),receipt)

if __name__=="__main__":
    unittest.main()
