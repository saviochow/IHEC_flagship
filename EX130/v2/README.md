# EX130 isoform-switch v2

Versioned methods and runnable Python source are packaged into 11 base64 text chunks under `bundle/`. To reconstruct:

```bash
cat EX130/v2/bundle/[0-9][0-9].b64 | base64 -d > EX130/v2/isoform_switch_v2.tar.gz
printf '%s  %s\n' '0aa9ed114e6f0c89077d8f2490f776c8f10b03152158cc93a6825edaa59bcdea' 'EX130/v2/isoform_switch_v2.tar.gz' | sha256sum -c -
tar -xzf EX130/v2/isoform_switch_v2.tar.gz
python3 -m unittest discover -s EX130/v2/tests -v
```

Prefer `scripts/bootstrap_savio_ex130_v2.sh` for a new exact-source SavioCluster checkout. The archive includes `isoform_switch.py`, `mcp_dispatch.py`, a full `README.md`, synthetic integration tests, and two Slurm scripts. No research data or results are bundled. The source treats transcription prediction and epigenomic correlation as descriptive rather than causal.
