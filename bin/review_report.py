#!/usr/bin/env python3
"""Prepare, annotate and render the optional offline review report."""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / 'review'))

from engine import prepare, complete, import_results
from resources import prepare_vcf, validate, make_manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version='oncoseq-review 1.0.0')
    commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('prepare')
    for name in ('sample', 'inputs', 'references', 'panel', 'output', 'payload'):
        p.add_argument('--' + name, required=True)
    p = commands.add_parser('complete')
    for name in ('sample', 'catalogue', 'payload', 'vep', 'references', 'panel', 'known-fusions', 'output'):
        p.add_argument('--' + name, required=True)
    p.add_argument('--distance', type=int, default=100)
    p.add_argument('--image-lock')
    p.add_argument('--code-digest')
    p = commands.add_parser('import')
    for name in ('results', 'sample', 'output'):
        p.add_argument('--' + name, required=True)
    p = commands.add_parser('prepare-vcf')
    for name in ('source', 'fasta', 'output'):
        p.add_argument('--' + name, required=True)
    p.add_argument('--skip-masked-reference', action='store_true',
                   help='Archive and exclude records over FASTA N bases; still reject other REF mismatches')
    p = commands.add_parser('validate')
    p.add_argument('--references', required=True)
    p.add_argument('--analytical-fasta')
    p.add_argument('--output', default='validated.json')
    p = commands.add_parser('manifest')
    p.add_argument('--resources', required=True)
    p.add_argument('--output', required=True)
    p = commands.add_parser('patch-alignment-qc')
    p.add_argument('--input', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--alignment-qc', required=True)
    p = commands.add_parser('render')
    p.add_argument('--input', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--template', default=str(Path(__file__).resolve().parents[1] / 'assets/review/report.html'))
    p.add_argument('--figures')
    p.add_argument('--results-root', help='Optional existing local output directory for presentation-only QC and NASVAR plots')
    args = vars(parser.parse_args())
    command = args.pop('command')
    if command == 'import':
        import_results(**args)
    elif command == 'prepare-vcf':
        prepare_vcf(**args)
    elif command == 'validate':
        _, manifest = validate(args['references'])
        if args['analytical_fasta'] and manifest['resources']['fasta'].get('sha256'):
            from resources import checksum
            if checksum(args['analytical_fasta']) != manifest['resources']['fasta']['sha256']:
                raise ValueError('Analytical FASTA checksum mismatch')
        Path(args['output']).write_text(json.dumps(manifest, indent=2) + '\n')
    elif command == 'manifest':
        make_manifest(**args)
    elif command == 'patch-alignment-qc':
        from alignment_qc import patch_html
        patch_html(**args)
    elif command == 'render':
        from render import render
        render(**args)
    else:
        globals()[command](**args)


if __name__ == '__main__':
    main()
