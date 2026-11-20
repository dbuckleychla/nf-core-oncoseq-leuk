#!/usr/bin/env python3
"""Check that every process has a digest pin and that the config matches the lock."""

import json
from pathlib import Path
import re
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


def check_release():
    lock = json.loads((ROOT / 'images.lock.json').read_text())
    if lock['platform'] != 'linux/amd64':
        raise ValueError('The clinical candidate targets linux/amd64')

    definitions = set()
    for source in (ROOT / 'modules').glob('**/main.nf'):
        if 'tests' in source.parts:
            continue
        definitions.update(re.findall(r'^process\s+(\w+)\s*\{', source.read_text(), re.MULTILINE))

    configured = set()
    expected = ['// Resolved from images.lock.json for linux/amd64.', 'process {']
    for source, entry in sorted(lock['images'].items()):
        if not re.fullmatch(r'[^\s@]+@sha256:[a-f0-9]{64}', entry['image']):
            raise ValueError(f'Image requires a SHA-256 digest: {source}')
        names = entry['processes']
        duplicates = configured.intersection(names)
        if duplicates:
            raise ValueError(f'Duplicate process assignments: {sorted(duplicates)}')
        configured.update(names)
        expected.extend([
            f'    // {source}',
            "    withName: '" + '|'.join(names) + "' {",
            "        container = '" + entry['image'] + "'",
            '    }',
        ])
    expected.extend(['}', ''])
    if configured != definitions:
        raise ValueError(f'Process coverage differs. Missing pins {sorted(definitions - configured)}. Missing definitions {sorted(configured - definitions)}')
    if (ROOT / 'conf/containers_clinical_amd64.config').read_text() != '\n'.join(expected):
        raise ValueError('Container config differs from images.lock.json')

    status = subprocess.check_output(['git', 'submodule', 'status', '--recursive'], cwd=ROOT, text=True)
    if any(line[0] != ' ' for line in status.splitlines()):
        raise ValueError('Initialise submodules at their recorded commits with git submodule update --init --recursive')
    print(f'Checked {len(configured)} processes and {len(lock["images"])} image digests for {lock["platform"]}')


if __name__ == '__main__':
    try:
        check_release()
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        print(f'Release check failed: {error}', file=sys.stderr)
        sys.exit(1)
