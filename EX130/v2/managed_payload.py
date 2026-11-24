#!/usr/bin/env python3
"""Render a self-contained Researchstack payload from reviewed EX130 v2 source.

No submission, approval, SSH, package installation or dataset modification.
The original source bundle remains the immutable scientific implementation.
"""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import shlex

EXPECTED = {
    "isoform_switch.py": "58840cbe5c177a37dd83a607f47d4381e607007b4d5dfae9dc92eafc78de4120",
    "slurm/run_single_allocation.sh": "fd54322db7a9db32e88667ef35be8ea31dd16dd0d8be2c8b1e9f4483126bfbd2",
}
MANIFEST_SHA = "411a795141301a87888175430ffd4248e497bc3133023434ed36921a74b004f7"

def sha(data):
    return hashlib.sha256(data).hexdigest()

def absolute(value):
    if (not isinstance(value, str) or not re.fullmatch(r"/[A-Za-z0-9_./+-]+", value)
            or str(PurePosixPath(value)) != value or ".." in PurePosixPath(value).parts
            or value == "/"):
        raise ValueError("require explicit canonical absolute path")
    return value

def checked_config(cfg):
    required = {"cluster", "ex130_root", "logits_dir", "gff", "features_root",
                "python", "source_model", "workers", "manifest_sha256"}
    if not isinstance(cfg, dict) or set(cfg) != required:
        raise ValueError("incorrect payload config fields")
    if cfg["cluster"] not in ("pines", "yipcluster"):
        raise ValueError("unknown site")
    for k in ("ex130_root", "logits_dir", "gff", "features_root", "python"):
        absolute(cfg[k])
    if not isinstance(cfg["source_model"], str) or not re.fullmatch(
            r"[A-Za-z0-9_.+-]{1,120}", cfg["source_model"]):
        raise ValueError("invalid source model")
    if type(cfg["workers"]) is not int or not 1 <= cfg["workers"] <= 8:
        raise ValueError("workers must be 1..8")
    if cfg["manifest_sha256"] != MANIFEST_SHA:
        raise ValueError("manifest identity differs from reviewed 358-sample manifest")
    return cfg

def embed(path, data):
    text = data.decode("utf-8")
    if not text.endswith("\n") or "\x00" in text:
        raise ValueError("embedded source must be newline-terminated UTF-8")
    delimiter = "EX130_SOURCE_" + sha(data)
    if delimiter in text.splitlines():
        raise ValueError("unexpected delimiter collision")
    return ("cat > " + path + " <<'" + delimiter + "'\n" + text + delimiter +
            "\n" + "printf '%s  %s\\n' " + sha(data) + " " + path +
            " | /usr/bin/sha256sum -c -\n")

def render(source_dir, cfg, artifact_source=None):
    cfg = checked_config(cfg)
    sources = {}
    for name, expected in EXPECTED.items():
        p = source_dir / name
        if p.is_symlink() or not p.is_file():
            raise ValueError("missing original source")
        data = p.read_bytes()
        if sha(data) != expected:
            raise ValueError("original source identity changed: " + name)
        sources[name] = data
    if artifact_source is None:
        artifact_source = Path(__file__).with_name("managed_artifacts.py").read_bytes()
    q = shlex.quote
    # The verified Researchstack runtime passes these values after its storage
    # identity, payload hash and allocation preflight. They are not inherited
    # from an operator's interactive shell.
    script = """#!/bin/bash
set -Eeuo pipefail
umask 077
: "$RESEARCHSTACK_CONTEXT" "$RESEARCHSTACK_CONTEXT_SHA256"
: "$RESEARCHSTACK_RUN_DIR" "$RESEARCHSTACK_TASK_DIR" "$SLURM_JOB_ID" "$SLURM_CPUS_PER_TASK"
export PATH=/usr/bin:/bin
export PYTHONDONTWRITEBYTECODE=1
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
[[ "$SLURM_JOB_ID" =~ ^[1-9][0-9]*$ ]]
[[ "$SLURM_CPUS_PER_TASK" =~ ^[1-9][0-9]*$ ]]
"""
    script += "export V2_PYTHON=" + q(cfg["python"]) + "\n"
    script += "[[ -x \"$V2_PYTHON\" ]]\n"
    script += "(( SLURM_CPUS_PER_TASK >= " + str(cfg["workers"]) + " ))\n"
    script += '"$V2_PYTHON" -I -B - ' + q(cfg["cluster"]) + " <<'EX130_CONTEXT_CHECK'\n"
    script += """import hashlib,json,os,sys
from pathlib import Path
p=Path(os.environ["RESEARCHSTACK_CONTEXT"])
assert p.is_file() and not p.is_symlink()
raw=p.read_bytes()
assert hashlib.sha256(raw).hexdigest()==os.environ["RESEARCHSTACK_CONTEXT_SHA256"]
c=json.loads(raw)
run=Path(os.environ["RESEARCHSTACK_RUN_DIR"])
task=Path(os.environ["RESEARCHSTACK_TASK_DIR"])
assert run.is_absolute() and run.resolve()==run and task.resolve()==task
assert c["cluster"]==sys.argv[1] and c["trials"]==1 and c.get("nodes",1)==1
assert Path(c["run_dir"])==run and task==run/"task-0"
assert p.parent==run and task.is_dir()
EX130_CONTEXT_CHECK
"""
    for var,key in (("V2_EX130_ROOT","ex130_root"),("V2_EX51_ROOT","logits_dir"),
                    ("V2_GFF","gff"),("V2_FEATURE_ROOT","features_root"),
                    ("V2_SOURCE_MODEL","source_model")):
        script += "export " + var + "=" + q(cfg[key]) + "\n"
    script += "export V2_WORKERS=" + str(cfg["workers"]) + "\nexport V2_RESERVED_CASES=35\n"
    script += "printf '%s  %s\\n' " + cfg["manifest_sha256"] + ' "$V2_EX130_ROOT/results/fixed_source/sample_manifest.csv" | /usr/bin/sha256sum -c -\n'
    script += """export V2_SOURCE_ROOT="$RESEARCHSTACK_TASK_DIR/ex130-source"
export V2_RUN_ROOT="$RESEARCHSTACK_TASK_DIR/ex130-v2"
[[ ! -e "$V2_SOURCE_ROOT" && ! -e "$V2_RUN_ROOT" ]]
mkdir -- "$V2_SOURCE_ROOT"
mkdir -p -- "$V2_SOURCE_ROOT/EX130/v2/slurm"
"""
    script += embed('"$V2_SOURCE_ROOT/EX130/v2/isoform_switch.py"', sources["isoform_switch.py"])
    script += embed('"$V2_SOURCE_ROOT/EX130/v2/slurm/run_single_allocation.sh"', sources["slurm/run_single_allocation.sh"])
    script += embed('"$V2_SOURCE_ROOT/EX130/v2/managed_artifacts.py"', artifact_source)
    script += """rc=0
/bin/bash "$V2_SOURCE_ROOT/EX130/v2/slurm/run_single_allocation.sh" || rc=$?
artifact_rc=0
"$V2_PYTHON" -I -B "$V2_SOURCE_ROOT/EX130/v2/managed_artifacts.py" \\
    --task-dir "$RESEARCHSTACK_TASK_DIR" --analysis-exit "$rc" || artifact_rc=$?
(( rc == 0 )) || exit "$rc"
exit "$artifact_rc"
"""
    payload = script.encode("utf-8")
    if len(payload) > 1024 * 1024:
        raise ValueError("payload exceeds Researchstack limit")
    return payload

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-dir", required=True, type=Path,
                   help="original extracted EX130/v2 directory")
    p.add_argument("--config", required=True, type=Path,
                   help="private reviewed site/input configuration JSON")
    p.add_argument("--output", required=True, type=Path)
    a=p.parse_args()
    raw=a.config.read_bytes()
    cfg=checked_config(json.loads(raw))
    payload=render(a.source_dir,cfg)
    a.output.mkdir(mode=0o700,parents=False,exist_ok=False)
    target=a.output/"payload.sh"
    with target.open("xb") as f: f.write(payload)
    receipt={"schema":"ihec.ex130-managed-payload/v1","payload_sha256":sha(payload),
             "config_sha256":sha(raw),"original_source_sha256":EXPECTED,
             "artifact_source_sha256":sha(Path(__file__).with_name("managed_artifacts.py").read_bytes()),
             "cluster":cfg["cluster"],"workers":cfg["workers"],
             "manifest_sha256":cfg["manifest_sha256"],
             "plan_created":False,"approval_authority":False,"submitted":False}
    with (a.output/"payload_receipt.json").open("x") as f:
        json.dump(receipt,f,indent=2,sort_keys=True)
        f.write("\n")
    print(json.dumps(receipt,sort_keys=True))

if __name__=="__main__":
    main()
