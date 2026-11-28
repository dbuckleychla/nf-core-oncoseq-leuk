#!/usr/bin/env python3
"""Verify a real QC acceptance run, not just output declarations."""
import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('results')
a = parser.parse_args()
root = Path(a.results) / 'qc_extended'
allowed = {'metrics.json', 'summary.tsv', 'read_lengths.tsv', 'coverage.tsv', 'targets.tsv',
           'methylation.tsv', 'classifier_cpgs.tsv', 'adequacy.tsv', 'preprocessing.tsv'}
for name in ('SYNTHETIC_A', 'SYNTHETIC_B'):
    folder = root / name
    assert {p.name for p in folder.iterdir()} == allowed
    data = json.loads((folder / 'metrics.json').read_text())
    assert data['sample_id'] == name
    assert data['components']['methylation']['status'] == 'disabled'
    assert next(r['value'] for r in data['metrics'] if r['metric_id'] == 'alignment.rate') == .75
assert {p.name for p in root.iterdir()} == {'SYNTHETIC_A', 'SYNTHETIC_B', 'run_summary.tsv', 'run_preprocessing.tsv'}
print('Real sample publication and private-intermediate exclusion: PASS')
