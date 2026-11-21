#!/usr/bin/env python3
"""Deterministic VEP boundary fixture; scientific adapters remain real."""
import argparse
import shutil
from pathlib import Path

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--help', action='store_true')
parser.add_argument('--dir_cache')
parser.add_argument('--output_file')
args, _ = parser.parse_known_args()
if args.help:
    print('ensembl-vep : 115')
else:
    shutil.copyfile(Path(args.dir_cache) / 'fixture.vep.vcf.gz', args.output_file)
