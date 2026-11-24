#!/usr/bin/env python3
"""Hash only explicit EX130 review outputs; never publish measurements to GitHub."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import stat
import os

def hash_regular(path):
    if path.is_symlink():
        raise ValueError("symlink artifact")
    fd=os.open(path,os.O_RDONLY|os.O_NOFOLLOW|os.O_NONBLOCK)
    with os.fdopen(fd,"rb") as f:
        before=os.fstat(f.fileno())
        if not stat.S_ISREG(before.st_mode) or before.st_size<=0:
            raise ValueError("invalid artifact")
        h=hashlib.sha256()
        for chunk in iter(lambda:f.read(1024*1024),b""):
            h.update(chunk)
        after=os.fstat(f.fileno())
        if (before.st_dev,before.st_ino,before.st_size,before.st_mtime_ns,before.st_ctime_ns)!=(
                after.st_dev,after.st_ino,after.st_size,after.st_mtime_ns,after.st_ctime_ns):
            raise ValueError("artifact changed")
    return before.st_size,h.hexdigest()

def read_json(p):
    if p.is_symlink() or p.stat().st_size > 16*1024*1024:
        raise ValueError("unsafe receipt")
    return json.loads(p.read_text())

def collect(task,analysis_exit):
    if not task.is_absolute() or task.resolve()!=task or not task.is_dir():
        raise ValueError("invalid task directory")
    out=task/"ex130-v2"/"analysis"
    files=[]
    missing=[]
    case_states=[]
    def add(rel,kind,required=True):
        path=task/rel
        if path.resolve()!=path or not path.is_file():
            if required: missing.append(rel)
            return False
        n,h=hash_regular(path)
        files.append({"path":rel,"kind":kind,"bytes":n,"sha256":h})
        return True
    prefix="ex130-v2/analysis/"
    for name in ("discovery_receipt.json","collection_receipt.json","cases.csv","candidate_validation.csv"):
        add(prefix+name,"review_metadata")
    cases=[]
    if (out/"cases.csv").is_file() and not (out/"cases.csv").is_symlink():
        with (out/"cases.csv").open(newline="") as f:
            cases=list(csv.DictReader(f))
    if len(cases)>35 or len({r.get("case_id") for r in cases})!=len(cases):
        raise ValueError("invalid candidate inventory")
    for i,c in enumerate(cases):
        cid=c.get("case_id","")
        if not re.fullmatch(r"SWITCH_[0-9]{3}",cid):
            raise ValueError("invalid case id")
        base=prefix+"case_runs/"+f"{i:03d}"+"/"
        if not add(base+"receipt.json","case_receipt"):
            case_states.append({"case_id":cid,"status":"MISSING_RECEIPT"})
            continue
        receipt=read_json(task/(base+"receipt.json"))
        if receipt.get("case_id")!=cid or receipt.get("index")!=i:
            raise ValueError("case receipt identity mismatch")
        status=receipt.get("status")
        if status not in ("VALIDATED","INSUFFICIENT_FEATURE_PAIRS"):
            raise ValueError("unknown case status")
        for suffix in ("png","svg"):
            add(base+"isoform_pair."+suffix,"figure_"+suffix,status=="VALIDATED")
        case_states.append({"case_id":cid,"status":status})
    completed=False
    r=out/"collection_receipt.json"
    if r.is_file():
        summary=read_json(r)
        completed=(summary.get("status")=="COMPLETE" and summary.get("expected_cases")==len(cases)
                   and summary.get("present_receipts")==len(cases))
    status=("COMPLETE" if analysis_exit==0 and completed and not missing and cases
            and all(c["status"]=="VALIDATED" for c in case_states) else
            "EMPTY_DISCOVERY" if analysis_exit==0 and completed and not missing and not cases
            else "INCOMPLETE")
    return {"schema":"ihec.ex130-artifacts/v1","status":status,
            "analysis_exit_code":analysis_exit,"expected_cases":len(cases),
            "cases":case_states,"files":files,"missing":missing,
            "scientific_support_claimed":False,"approval_authority":False,
            "raw_measurements_exported":False}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--task-dir",required=True,type=Path)
    p.add_argument("--analysis-exit",required=True,type=int)
    a=p.parse_args()
    report=collect(a.task_dir,a.analysis_exit)
    # This manifest is part of the private managed result directory, not a
    # GitHub upload or a replacement for authoritative scheduler accounting.
    with (a.task_dir/"ex130_artifact_manifest.json").open("x") as f:
        json.dump(report,f,indent=2,sort_keys=True)
        f.write("\n")
    print(json.dumps({"status":report["status"],"artifacts":len(report["files"]),
                      "missing":len(report["missing"])}))
    return 0 if report["status"] in ("COMPLETE","EMPTY_DISCOVERY") else 2

if __name__=="__main__":
    raise SystemExit(main())
