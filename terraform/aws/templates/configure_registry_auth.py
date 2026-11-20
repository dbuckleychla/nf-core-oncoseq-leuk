#!/usr/bin/env python3
"""Worker-only GHCR auth: retrieve a secret without exposing it to logs/state."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile


def configure(aws_cli, region, secret_arn, config):
    # Capture both streams: AWS error responses must not leak credentials either.
    result = subprocess.run(
        [aws_cli, 'secretsmanager', 'get-secret-value', '--region', region,
         '--secret-id', secret_arn, '--query', 'SecretString', '--output', 'text',
         '--no-cli-pager'], capture_output=True, text=True, timeout=60)
    if result.returncode:
        raise ValueError('Cannot retrieve GHCR secret; check worker IAM, KMS and connectivity')
    try:
        secret = json.loads(result.stdout)
    except (ValueError, TypeError):
        raise ValueError('GHCR secret must contain a JSON object') from None
    if not isinstance(secret, dict):
        raise ValueError('GHCR secret must contain a JSON object')
    username = secret.get('username')
    # Accept MSQ/ECR pull-through credentials and the ECS username/password shape.
    password = secret.get('accessToken', secret.get('password'))
    if ('accessToken' in secret and 'password' in secret
            and secret['accessToken'] != secret['password']):
        raise ValueError('GHCR secret contains conflicting token fields')
    for value in (username, password):
        if not isinstance(value, str) or not value or any(c.isspace() or ord(c) < 32 for c in value):
            raise ValueError('GHCR secret requires nonempty username and accessToken (or password)')
    auth = json.dumps({'ghcr.io': {'username': username, 'password': password}},
                      separators=(',', ':'), ensure_ascii=True)
    # ECS init reads ecs.config when starting the agent. Preserve the cluster and
    # other Batch-managed settings; replace only our two registry-auth settings.
    current = config.read_text() if config.exists() else ''
    lines = [line for line in current.splitlines()
             if not re.match(r'^\s*(?:export\s+)?ECS_ENGINE_AUTH_(?:TYPE|DATA)\s*=', line)]
    content = '\n'.join(lines + ['ECS_ENGINE_AUTH_TYPE=docker',
                                'ECS_ENGINE_AUTH_DATA=' + auth]) + '\n'
    descriptor, temporary = tempfile.mkstemp(prefix='.registry-auth-', dir=config.parent)
    try:
        with os.fdopen(descriptor, 'w') as handle:
            os.fchmod(handle.fileno(), 0o600)
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, config)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--aws-cli', required=True)
    parser.add_argument('--region', required=True)
    parser.add_argument('--secret-arn', required=True)
    parser.add_argument('--config', type=Path, default=Path('/etc/ecs/ecs.config'))
    args = parser.parse_args()
    try:
        configure(args.aws_cli, args.region, args.secret_arn, args.config)
    except (OSError, ValueError, subprocess.SubprocessError):
        # Do not print exception bodies, subprocess output, or the secret.
        print('GHCR authentication setup failed; check secret format, worker IAM/KMS and network access. ECS startup blocked.', file=sys.stderr)
        return 1
    print('GHCR authentication configured for ECS')
    return 0


if __name__ == '__main__':
    sys.exit(main())
