#!/usr/bin/env bash
# Execute fixed read-only preflight using existing operator SSH. No root or Slurm.
set -Eeuo pipefail
umask 077
[[ "$(hostname -s)" == saviocluster && "$(id -u)" != 0 ]]
SOURCE="$(dirname -- "$0")/ex130_yip_login_preflight.py"
[[ -f "$SOURCE" && ! -L "$SOURCE" ]]
WORK="$(mktemp -d "$HOME/ex130-yip-login-preflight.XXXXXXXX")"
rc=0
timeout -k 5s 180s ssh -T -o BatchMode=yes -o StrictHostKeyChecking=yes \
  -o ForwardAgent=no -o ConnectTimeout=10 -o ConnectionAttempts=1 \
  -o ControlMaster=no -o ControlPath=none hchow@yipcluster \
  '/usr/bin/python3 -I -B -' <"$SOURCE" >"$WORK/private.json" 2>"$WORK/private.stderr" || rc=$?
/usr/bin/python3 -I -B - "$WORK/private.json" "$rc" >"$WORK/receipt.json" <<'PY'
import datetime,json,re,sys
r={"schema":"ihec.ex130-yip-login-preflight/v1","status":"INVALID_REMOTE_RESPONSE",
   "read_only":True,"jobs_submitted":0,"compute_storage_verified":False,"approval_authority":False}
try:
    x=json.load(open(sys.argv[1]))
    assert x.get("schema")==r["schema"] and x.get("jobs_submitted")==0 and x.get("read_only") is True
    for k in ("status","error_class","source_model"):
        if isinstance(x.get(k),str) and re.fullmatch(r"[A-Za-z0-9_.+-]{1,100}",x[k]): r[k]=x[k]
    for k in ("staged_sources_match","manifest_matches","gff_readable","python_executable_present"):
        if type(x.get(k)) is bool: r[k]=x[k]
    for k in ("manifest_records","eligible_distinct_samples","logits_present","logits_header_verified",
              "h3_feature_header_verified","wgbs_feature_header_verified","both_feature_headers_verified","gff_bytes"):
        if type(x.get(k)) is int and 0<=x[k]<=10**12: r[k]=x[k]
    v=x.get("runtime_package_metadata",{})
    r["runtime_package_metadata"]={k:y for k,y in v.items() if k in ("python","numpy","matplotlib") and isinstance(y,str) and re.fullmatch(r"[A-Za-z0-9.+_-]{1,80}",y)}
    v=x.get("shared_storage",{})
    n=v.get("matching_views")
    if type(n) is int and 0<=n<=20:
        r["matching_yip_storage_views"]=n
        p=v.get("verified_yip_mount")
        if n==1 and isinstance(p,str) and re.fullmatch(r"/autofs/login001_ssd/hchow/pines_mount(?:/[A-Za-z0-9_.-]+)*",p) and ".." not in p.split("/"):
            r["verified_yip_mount"]=p
except Exception:
    r["status"]="INVALID_REMOTE_RESPONSE"
r["exit_code"]=int(sys.argv[2])
r["utc"]=datetime.datetime.now(datetime.timezone.utc).isoformat()
print(json.dumps(r,sort_keys=True))
PY
NAME="ex130-yip-login-preflight-$(date -u +%Y%m%dT%H%M%SZ)-$(cat /proc/sys/kernel/random/uuid).json"
if timeout -k 3s 45s scp -q -o BatchMode=yes -o StrictHostKeyChecking=yes \
    -o ForwardAgent=no -o ConnectTimeout=10 -- "$WORK/receipt.json" \
    "hchow@yipcluster:/autofs/login001_ssd/hchow/scratch/clustermcp-runs/$NAME"; then
    printf 'CLUSTERMCP_READABLE_ROOT=yip_ssd\nCLUSTERMCP_READABLE_PATH=scratch/clustermcp-runs/%s\n' "$NAME"
else
    echo CLUSTERMCP_READABLE_REPORT=UPLOAD_FAILED_LOCAL_COPY_RETAINED
fi
cat "$WORK/receipt.json"
echo "LOCAL_PRIVATE_RUN=$WORK"
exit "$rc"
