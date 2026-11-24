#!/usr/bin/env python3
"""Minimal, idempotent Python 3.8 argparse compatibility patch for EX130.

Run with /usr/bin/python3 to avoid importing DGL or any scientific packages.
Does not alter model files, manifests, results, or graph tensors.
"""
from pathlib import Path

source = Path(__file__).resolve().parent / 'transcript_attribution.py'
if not source.exists():
    raise SystemExit('Expected transcript_attribution.py in the same directory as this helper')
old = 'ap.add_argument("--self-loop", action=argparse.BooleanOptionalAction, default=True)'
new = '''ap.add_argument("--self-loop", dest="self_loop", action="store_true", default=True)
    ap.add_argument("--no-self-loop", dest="self_loop", action="store_false")'''
content = source.read_text()
if old in content:
    if content.count(old) != 1:
        raise SystemExit('Multiple occurrences: refusing to alter ambiguous source')
    backup = source.with_suffix('.py.before_py38_patch')
    if backup.exists():
        raise SystemExit('Backup already exists; refusing to overwrite')
    backup.write_text(content)
    source.write_text(content.replace(old, new))
    print('EX130_PY38_ARGPARSE_PATCH=APPLIED')
elif new in content:
    print('EX130_PY38_ARGPARSE_PATCH=ALREADY_PRESENT')
else:
    raise SystemExit('Unexpected source version; no changes made')
