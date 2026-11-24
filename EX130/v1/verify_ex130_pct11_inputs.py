#!/usr/bin/env python3
"""Validate raw IHEC methylation percentage column and optionally cache missing BED indexes.

Never modifies source .bed.gz / BigWig data. Python 3.8-compatible, stdlib-only.
"""
import argparse
import csv
import gzip
import json
import shutil
import subprocess
from pathlib import Path

from ex130_raw_tracks import locate_tracks


def verify_header(path):
    """Return a schema/first-record report after reading the beginning of one gzip file."""
    with gzip.open(str(path), 'rt', encoding='utf-8', errors='replace') as handle:
        header = handle.readline().strip().split('\t')
        if len(header) < 11 or header[:3] != ['seqnames', 'start', 'end'] or header[10] != 'percentage':
            raise ValueError('WGBS header is not seqnames/start/end + percentage in one-based column 11: {}'.format(header[:15]))
        if header[9] != 'cov':
            raise ValueError('WGBS column 10 is not cov: {}'.format(header[:11]))
        first = None
        for _ in range(16):
            line = handle.readline()
            if not line:
                break
            if line.strip() and not line.startswith('#'):
                first = line.rstrip('\n').split('\t')
                break
        if first is None:
            raise ValueError('No first WGBS CpG record')
        if len(first) < len(header):
            raise ValueError('First record shorter than schema')
        try:
            start, end, coverage, percentage = int(first[1]), int(first[2]), float(first[9]), float(first[10])
        except (TypeError, ValueError) as exc:
            raise ValueError('Invalid BED coordinates or CpG percentage: {}'.format(exc))
        if not (0 <= start < end and coverage >= 0 and 0 <= percentage <= 100):
            raise ValueError('Invalid first record ranges: start={}, end={}, coverage={}, percentage={}'.format(start, end, coverage, percentage))
        return {'column_one_based': 11, 'column_name': 'percentage', 'scale': 'percent',
                'beta_formula': 'percentage / 100', 'first_percentage': percentage,
                'first_coverage': coverage, 'header_fields': len(header)}


def validate_and_cache(root, outdir, main_raw_root, partial_raw_root, build_indexes, tabix_bin):
    root, outdir = Path(root), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sample_csv = root / 'track_sample_manifest.csv'
    with sample_csv.open(newline='') as f:
        rows = list(csv.DictReader(f))
    samples = sorted(set(r['sample_id'] for r in rows))
    if not samples:
        raise RuntimeError('No track samples found')
    cache = outdir / 'wgbs_index_cache'
    checked = []
    failures = []
    for sample in samples:
        record = {'sample_id': sample}
        try:
            found = locate_tracks(sample, main_raw_root, ['H3K36me3'],
                                  secondary_root=partial_raw_root, index_cache_dir=cache)
            raw = Path(found['wgbs_path']) if found['wgbs_path'] else None
            if raw is None or not raw.is_file():
                raise ValueError('WGBS whole-region BED missing')
            record.update({'path': str(raw), 'schema': verify_header(raw),
                           'original_track_status': found['status'],
                           'index_present_before': bool(found['wgbs_index'])})
            if not found['wgbs_index'] and build_indexes:
                index_path = cache / found['base_id'] / raw.name
                index_path.parent.mkdir(parents=True, exist_ok=True)
                if index_path.exists() and not index_path.is_symlink():
                    raise RuntimeError('Refusing to overwrite existing non-symlink cache file')
                if index_path.is_symlink() and index_path.resolve() != raw.resolve():
                    raise RuntimeError('Refusing to alter mismatched index symlink')
                if not index_path.is_symlink():
                    index_path.symlink_to(raw.resolve())
                tabix = shutil.which(tabix_bin)
                if not tabix:
                    record['index_cache_status'] = 'TABIX_NOT_INSTALLED'
                else:
                    # The non-comment first line is a BED header: -S 1 is required.
                    # Index is created ONLY beside the symlink in the EX130 output tree.
                    proc = subprocess.run([tabix, '-S', '1', '-p', 'bed', str(index_path)],
                                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                          universal_newlines=True)
                    idx = Path(str(index_path) + '.tbi')
                    csi = Path(str(index_path) + '.csi')
                    record['index_cache_status'] = ('CREATED' if proc.returncode == 0 and
                                                    (idx.exists() or csi.exists()) else 'FAILED')
                    record['index_stderr'] = proc.stderr.strip()[:500]
            fresh = locate_tracks(sample, main_raw_root, ['H3K36me3'],
                                  secondary_root=partial_raw_root, index_cache_dir=cache)
            record['index_available_after'] = bool(fresh['wgbs_index'])
            record['index_status_after'] = fresh['status']
        except (OSError, ValueError, RuntimeError) as exc:
            record['schema_error'] = '{}: {}'.format(type(exc).__name__, exc)
            failures.append(sample)
        checked.append(record)
    proof = {'n_samples': len(samples), 'valid_percentage_column_11': len(samples)-len(failures),
             'schema_errors': failures, 'all_headers_verified': not failures,
             'methylation_beta_formula': 'percentage / 100', 'samples': checked}
    (outdir / 'percentage_column_audit.json').write_text(json.dumps(proof, indent=2) + '\n')
    print('Verified WGBS column 11 (percentage) in {}/{} samples'.format(len(samples)-len(failures),len(samples)), flush=True)
    print('Schema audit: {}'.format(outdir / 'percentage_column_audit.json'), flush=True)
    if failures:
        raise RuntimeError('Refusing to plot due to unsupported WGBS header: {}'.format(failures))
    return proof


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ex130-dir', type=Path, default=Path.home() / 'IHEC_Sept_2025/EX130')
    p.add_argument('--outdir', type=Path, required=True)
    p.add_argument('--build-index-cache', action='store_true')
    p.add_argument('--tabix-bin', default='tabix')
    a = p.parse_args()
    validate_and_cache(a.ex130_dir, a.outdir,
        Path.home() / 'data_store/IHEC_data_processing_Aug_2023',
        Path.home() / 'data_store/IHEC_data_processing_Aug_2023_partial_epigenomes',
        a.build_index_cache, a.tabix_bin)


if __name__ == '__main__':
    main()
