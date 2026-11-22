#!/usr/bin/env python3
"""Scan every outgoing commit patch; reserve full history for audits/unknown baselines."""
import argparse
import subprocess
import sys


def scan_range(head, upstream=None, full_history=False):
    return '--all' if full_history or not upstream else upstream + '..' + head


def check_outgoing_content(selection):
    from check_public_content import findings
    objects = subprocess.check_output(['git', 'rev-list', '--objects', selection], text=True)
    proc = subprocess.Popen(['git', 'cat-file', '--batch'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    blocked = False
    try:
        for line in objects.splitlines():
            oid, _, path = line.partition(' ')
            proc.stdin.write((oid + '\n').encode())
            proc.stdin.flush()
            header = proc.stdout.readline().split()
            if len(header) != 3:
                raise RuntimeError('Cannot inspect outgoing Git object')
            content = proc.stdout.read(int(header[2]))
            proc.stdout.read(1)
            hits = findings(path.encode())
            if header[1] in (b'blob', b'commit', b'tag'):
                hits.extend(findings(content))
            if hits:
                print('Blocked outgoing object ' + oid[:12] + ': ' + ', '.join(sorted(set(hits))), file=sys.stderr)
                blocked = True
    finally:
        proc.stdin.close()
        proc.wait()
    return not blocked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full-history', action='store_true', help='Audit every reachable Git reference')
    args = parser.parse_args()
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    result = subprocess.run(['git', 'rev-parse', '--verify', '--quiet', '@{upstream}^{commit}'],
                            text=True, capture_output=True)
    upstream = result.stdout.strip() if result.returncode == 0 else None
    selection = scan_range(head, upstream, args.full_history)
    count = int(subprocess.check_output(['git', 'rev-list', '--count', selection], text=True))
    if not count:
        print('Secret scan: no outgoing commits.')
        return 0
    scope = 'full history (audit or unknown baseline)' if selection == '--all' else 'outgoing commits'
    print(f'Secret scan: {count} {scope}.', flush=True)
    if not check_outgoing_content(selection):
        return 1
    # Scan all intermediate patches, not just the final net diff: removed secrets
    # in an earlier outgoing commit must still block publication.
    return subprocess.call(['gitleaks', 'git', '.', '--redact', '--no-banner',
                            '--ignore-gitleaks-allow', '--log-opts=' + selection])


if __name__ == '__main__':
    sys.exit(main())
