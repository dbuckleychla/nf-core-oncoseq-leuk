#!/usr/bin/env python3
"""Emit small valid unaligned BAMs to test orchestration without GPU inference."""
import os
from pathlib import Path
import subprocess
import sys

if sys.argv[1:] == ['--version']:
    print('dorado mock')
    raise SystemExit(0)

sample = os.environ['MOCK_SAMPLE']
shard = int(os.environ['MOCK_SHARD'])
if os.environ.get('MOCK_FAIL_SAMPLE') == sample and os.environ.get('MOCK_FAIL_SHARD') == str(shard):
    raise SystemExit(1)

reads = [read for path in sorted(Path('pod5_inputs').glob('*.pod5')) for read in path.read_text().splitlines()]
header = ['@HD\tVN:1.6\tSO:unknown', f'@RG\tID:rg_{sample}\tSM:{sample}', '@PG\tID:dorado\tPN:dorado\tVN:mock']
records = [f'{read}\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tFFFF\tRG:Z:rg_{sample}\tMM:Z:C+m,0;\tML:B:C,200' for read in reads]
subprocess.run(['samtools', 'view', '-b', '-'], input='\n'.join(header + records) + '\n', text=True, stdout=sys.stdout.buffer, check=True)
