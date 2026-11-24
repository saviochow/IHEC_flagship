#!/usr/bin/env python3
"""Read-only login-side EX130 and shared-storage inventory; no Slurm submission."""
import csv
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys

MANIFEST_SHA="411a795141301a87888175430ffd4248e497bc3133023434ed36921a74b004f7"
MARKER_SHA="0dc9b935aeb1e6114022cfb3ea98c4710dd4702dd499d40c1c746772517afd6d"
SOURCE_MODEL="IHECRE00001900_mRNA-Seq"
PINEROOT=PurePosixPath("/data/labs/yip_lab/schow/Agent_projects/clustermcp-runs")
BASE=Path("/autofs/login001_ssd/hchow/pines_mount")
STAGED=Path("/autofs/login001_ssd/hchow/scratch/clustermcp-runs/IHEC_flagship_source_176a7aeb49ba0f7bf0050ff0340451baf5c47dd0/EX130/v2")
EXPECTED={"isoform_switch.py":"58840cbe5c177a37dd83a607f47d4381e607007b4d5dfae9dc92eafc78de4120",
          "slurm/run_single_allocation.sh":"fd54322db7a9db32e88667ef35be8ea31dd16dd0d8be2c8b1e9f4483126bfbd2"}

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def bounded(path,n):
    with path.open("rb") as f:
        v=f.read(n+1)
    if len(v)>n:
        raise ValueError("oversized bounded input")
    return v

def readable(path):
    return path.is_file() and os.access(path,os.R_OK) and path.stat().st_size>0

def header(path,delimiter,required):
    try:
        with path.open("r",encoding="utf-8",newline="") as f:
            line=f.readline(65537)
        return len(line)<=65536 and required.issubset(next(csv.reader([line],delimiter=delimiter)))
    except (OSError,UnicodeError,ValueError,StopIteration):
        return False

def mount_row(path):
    p=subprocess.run(["findmnt","--json","--target",str(path),
        "--output","TARGET,SOURCE,FSTYPE"],capture_output=True,timeout=10)
    if p.returncode:
        return None
    rows=json.loads(p.stdout).get("filesystems",[])
    return rows[0] if len(rows)==1 else None

def storage_view():
    candidates={BASE,BASE/"clustermcp-runs",BASE/"Agent_projects/clustermcp-runs"}
    matches=[]
    try:
        m=mount_row(BASE)
        if m:
            match=re.fullmatch(r"(?:schow@)?(?:pines|pines\.burnham\.org):(/[A-Za-z0-9_./-]+)",m.get("source",""))
            if match:
                root=PurePosixPath(match.group(1))
                if PINEROOT==root or root in PINEROOT.parents:
                    p=Path(m["target"])/PINEROOT.relative_to(root)
                    if p==BASE or BASE in p.parents:
                        candidates.add(p)
        for p in sorted(candidates):
            try:
                if not p.is_dir() or p.resolve()!=p:
                    continue
                m=mount_row(p)
                if not m or m.get("fstype") not in ("nfs","nfs4","fuse.sshfs"):
                    continue
                match=re.fullmatch(r"(?:schow@)?(?:pines|pines\.burnham\.org):(/[A-Za-z0-9_./-]+)",m.get("source",""))
                if not match:
                    continue
                rel=p.relative_to(Path(m["target"]))
                if PurePosixPath(match.group(1))/PurePosixPath(rel.as_posix())!=PINEROOT:
                    continue
                marker=p/".researchstack-storage-id"
                if marker.is_symlink() or not marker.is_file():
                    continue
                if sha(bounded(marker,256).strip())==MARKER_SHA:
                    matches.append(str(p))
            except (OSError,ValueError,subprocess.SubprocessError):
                continue
    except (OSError,ValueError,subprocess.SubprocessError):
        pass
    return {"matching_views":len(matches),"verified_yip_mount":matches[0] if len(matches)==1 else None,
            "compute_storage_verified":False}

def main():
    home=Path.home()
    manifest=home/"IHEC_Sept_2025/EX130/results/fixed_source/sample_manifest.csv"
    logits=home/"IHEC_Mar_2024/EX51/logits_debug_pc1"
    gff=home/"IHEC_Mar_2023/networks/g_b/gencode.v29.annotation.gff3"
    features=Path("/home/shared/IHEC/baseline_expt/features")
    py=home/"miniconda3/envs/DGL/bin/python3"
    r={"schema":"ihec.ex130-yip-login-preflight/v1","status":"PREFLIGHT_INCOMPLETE",
       "jobs_submitted":0,"compute_storage_verified":False,"approval_authority":False,
       "read_only":True,"source_model":SOURCE_MODEL}
    try:
        r["staged_sources_match"]=all(sha(bounded(STAGED/k,65536))==v for k,v in EXPECTED.items())
        raw=bounded(manifest,2*1024*1024)
        r["manifest_matches"]=sha(raw)==MANIFEST_SHA
        if not r["manifest_matches"]:
            raise ValueError("manifest")
        rows=list(csv.DictReader(raw.decode("utf-8").splitlines()))
        r["manifest_records"]=len(rows)
        seen=set(); selected=[]
        for row in rows:
            if row["rna_type"]!="mRNA" or row["disease"]!="non-disease" or not row["tissue"] or not row["base_id"]:
                continue
            if not all(re.fullmatch(r"[A-Za-z0-9_.-]+",row[k]) and row[k] not in (".","..") for k in ("sample_id","base_id")):
                raise ValueError("invalid identifier")
            key=(row["base_id"],row["rna_type"])
            if key in seen:
                continue
            seen.add(key); selected.append(row)
        r["eligible_distinct_samples"]=len(selected)
        r["logits_present"]=r["logits_header_verified"]=0
        r["h3_feature_header_verified"]=r["wgbs_feature_header_verified"]=0
        r["both_feature_headers_verified"]=0
        for row in selected:
            sid,base=row["sample_id"],row["base_id"]
            p=logits/("10sample_+thyroid_"+SOURCE_MODEL+"_"+sid+"_128_test_logits.csv")
            if readable(p):
                r["logits_present"]+=1
                r["logits_header_verified"]+=int(header(p,",",{"txid","chr","labels","logits"}))
            folder=features/base
            h3=sorted(folder.glob(base+".*.H3K36me3.pval.all.txt"))
            wgbs=sorted(folder.glob(base+".*.*WGBS.cpg.all.txt"))
            standard=[p for p in wgbs if ".standard-WGBS." in p.name]
            if len(standard)==1: wgbs=standard
            a=len(h3)==1 and header(h3[0],"\t",{"transcript.id","hg38H3K36me3intnintrons"})
            b=len(wgbs)==1 and header(wgbs[0],"\t",{"transcript.id","beta_hg38intnintrons"})
            r["h3_feature_header_verified"]+=int(a)
            r["wgbs_feature_header_verified"]+=int(b)
            r["both_feature_headers_verified"]+=int(a and b)
        r["gff_readable"]=readable(gff)
        r["gff_bytes"]=gff.stat().st_size if r["gff_readable"] else 0
        r["python_executable_present"]=py.is_file() and os.access(py,os.X_OK)
        if r["python_executable_present"]:
            code="import sys,json,importlib.metadata as m; print(json.dumps({'python':sys.version.split()[0], 'numpy':m.version('numpy'), 'matplotlib':m.version('matplotlib')}))"
            p=subprocess.run([str(py),"-I","-B","-c",code],capture_output=True,timeout=20)
            if p.returncode==0:
                v=json.loads(p.stdout)
                if all(isinstance(x,str) and re.fullmatch(r"[A-Za-z0-9.+_-]{1,80}",x) for x in v.values()):
                    r["runtime_package_metadata"]=v
        r["shared_storage"]=storage_view()
        r["status"]="LOGIN_PREFLIGHT_RECORDED"
    except Exception as e:
        r["error_class"]=type(e).__name__
    print(json.dumps(r,sort_keys=True))
    return 0 if r["status"]=="LOGIN_PREFLIGHT_RECORDED" else 2

if __name__=="__main__":
    sys.exit(main())
