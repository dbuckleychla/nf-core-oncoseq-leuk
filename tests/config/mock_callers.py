#!/usr/bin/env python3
"""Capture actual module CLI arguments without running inference."""
from pathlib import Path
import sys

args = sys.argv[1:]
if args == ['-v']:
    print('ClairS-TO mock')
elif Path(sys.argv[0]).name == 'run_clairs_to':
    output = Path(args[args.index('-o') + 1])
    output.mkdir()
    for name in ['snv.calls.vcf.gz', 'indel.calls.vcf.gz']:
        (output / name).write_text('\n'.join(args) + '\n')
else:
    for prefix in ['--output_vcf=', '--output_gvcf=']:
        output = Path(next(arg[len(prefix):] for arg in args if arg.startswith(prefix)))
        output.write_text('\n'.join(args) + '\n')
        Path(str(output) + '.tbi').touch()
