#!/usr/bin/env python3
"""Reject recognizable credentials and clinical accessions without printing them.

This is a conservative guard, not a general PHI detector. Gitleaks and manual
review supplement it. Reference paths and synthetic fixture IDs are allowed.
"""
import argparse
import gzip
import re
import subprocess
import sys

PATTERNS = {
    'clinical accession': re.compile(rb'\bM\d{2}(?:[-_]|%2[Dd]|\\u002[dD])\d{1,7}\b'),
    'AWS access key': re.compile(rb'\b(?:AKIA|ASIA)[A-Z0-9]{16}\b'),
    'GitHub token': re.compile(rb'\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,})'),
    'GitLab token': re.compile(rb'\bglpat-[A-Za-z0-9_-]{20,}'),
    'private key': re.compile(rb'-----BEGIN (?:RSA |EC |DSA |OPENSSH )?PRIVATE KEY-----'),
}


def findings(content):
    if content.startswith(b'\x1f\x8b'):
        try:
            content = gzip.decompress(content)
        except (OSError, EOFError):
            raise ValueError('Compressed tracked content could not be inspected')
    return [name for name, pattern in PATTERNS.items() if pattern.search(content)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--staged', action='store_true', help='Inspect the Git index (all staged tracked files)')
    parser.parse_args()
    entries = subprocess.check_output(['git', 'ls-files', '--stage', '-z']).split(b'\0')
    failures = []
    inspected = 0
    proc = subprocess.Popen(['git', 'cat-file', '--batch'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    try:
        for entry in entries:
            if not entry:
                continue
            details, _path = entry.split(b'\t', 1)
            mode, oid, stage = details.split()
            if mode == b'160000':
                continue
            if stage != b'0':
                raise ValueError('Resolve merge conflicts before public-content review')
            proc.stdin.write(oid + b'\n')
            proc.stdin.flush()
            header = proc.stdout.readline().split()
            if len(header) != 3 or header[1] != b'blob':
                raise ValueError('Cannot inspect a tracked Git object')
            content = proc.stdout.read(int(header[2]))
            proc.stdout.read(1)
            inspected += 1
            hits = sorted(set(findings(content) + findings(_path)))
            if hits:
                # Do not echo filenames: they may themselves contain clinical IDs.
                failures.append((oid.decode()[:12], hits))
    finally:
        proc.stdin.close()
        proc.wait()
    for oid, hits in failures:
        print('Blocked tracked blob ' + oid + ': ' + ', '.join(hits), file=sys.stderr)
    if failures:
        print('Remove sensitive content; do not bypass the check.', file=sys.stderr)
        return 1
    print(f'Public-content guard: {inspected} tracked blobs checked; no recognized sensitive patterns.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
