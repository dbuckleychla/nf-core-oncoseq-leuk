import importlib.util
import base64
import json
from pathlib import Path
import stat
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'terraform/aws/templates/configure_registry_auth.py'
spec = importlib.util.spec_from_file_location('configure_registry_auth', SCRIPT)
auth = importlib.util.module_from_spec(spec)
spec.loader.exec_module(auth)


class RegistryAuthTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which('terraform'), 'Terraform is needed for template rendering')
    def test_actual_bootstrap_renders_without_provider_or_cloud_access(self):
        for existing in (False, True):
            with self.subTest(existing=existing), tempfile.TemporaryDirectory() as tmp:
                variables = {
                    'aws_cli_version': '2.31.0', 'use_existing_aws_cli': existing,
                    'existing_aws_cli_path': '/custom/aws/bin/aws',
                    'host_agent_config': '{}', 'region': 'us-west-2',
                    'ghcr_secret_arn': 'arn:aws:secretsmanager:us-west-2:123456789012:secret:test-Ab1234',
                    'registry_auth_script': SCRIPT.read_text(),
                }
                expression = 'base64encode(templatefile(%s, jsondecode(%s)))' % (
                    json.dumps(str(SCRIPT.parent / 'user_data.mime')),
                    json.dumps(json.dumps(variables)))
                rendered = subprocess.run(['terraform', 'console'], cwd=tmp,
                                          input=expression + '\n', capture_output=True, text=True)
                self.assertEqual(rendered.returncode, 0, rendered.stderr)
                content = base64.b64decode(json.loads(rendered.stdout)).decode()
                self.assertLess(len(content.encode()), 16384)
                shell = content.split('#!/bin/bash', 1)[1].rsplit('--ONCOSEQ_LEUK--', 1)[0]
                check = subprocess.run(['/bin/bash', '-n'], input=shell, capture_output=True, text=True)
                self.assertEqual(check.returncode, 0, check.stderr)
                expected_cli = '/custom/aws/bin/aws' if existing else '/opt/aws-cli/v2/current/bin/aws'
                self.assertIn('--aws-cli ' + expected_cli, content)
                self.assertLess(content.index('20-ghcr-auth.conf'), content.index('systemctl start ecs'))
                self.assertIn(SCRIPT.read_text(), content)

    def test_msq_and_ecs_formats_preserve_config_and_rotate_credentials(self):
        for field in ('accessToken', 'password'):
            with self.subTest(field=field), tempfile.TemporaryDirectory() as tmp:
                config = Path(tmp) / 'ecs.config'
                config.write_text('ECS_CLUSTER=batch-cluster\nECS_IMAGE_PULL_BEHAVIOR=once\n'
                                  'ECS_ENGINE_AUTH_TYPE=dockercfg\n'
                                  'export ECS_ENGINE_AUTH_DATA=old-token\n')
                with patch.object(auth.subprocess, 'run') as run:
                    for token in ('first-token', 'rotated-token'):
                        run.return_value = subprocess.CompletedProcess([], 0, json.dumps(
                            {'username': 'test-user', field: token}), '')
                        auth.configure('/opt/aws', 'us-west-2', 'test-secret', config)
                        result = config.read_text()
                        self.assertIn('ECS_CLUSTER=batch-cluster\n', result)
                        self.assertIn('ECS_IMAGE_PULL_BEHAVIOR=once\n', result)
                        self.assertEqual(result.count('ECS_ENGINE_AUTH_DATA='), 1)
                        self.assertIn('ECS_ENGINE_AUTH_TYPE=docker\n', result)
                        data = result.split('ECS_ENGINE_AUTH_DATA=', 1)[1].strip()
                        self.assertEqual(json.loads(data), {'ghcr.io': {
                            'username': 'test-user', 'password': token}})
                        self.assertEqual(stat.S_IMODE(config.stat().st_mode), 0o600)
                    self.assertNotIn('first-token', result)
                    self.assertEqual(run.call_args.args[0], [
                        '/opt/aws', 'secretsmanager', 'get-secret-value', '--region',
                        'us-west-2', '--secret-id', 'test-secret', '--query',
                        'SecretString', '--output', 'text', '--no-cli-pager'])
                    self.assertEqual(list(Path(tmp).iterdir()), [config])

    def test_bad_secret_and_retrieval_failures_leave_existing_config_intact(self):
        cases = ['not JSON', 'null', '[]', '{}',
                 json.dumps({'username': 'u', 'accessToken': 'bad\ntoken'}),
                 json.dumps({'username': 'u', 'accessToken': 123}),
                 json.dumps({'username': 'u', 'accessToken': 'a', 'password': 'b'})]
        with tempfile.TemporaryDirectory() as tmp:
            config = Path(tmp) / 'ecs.config'
            config.write_text('ECS_CLUSTER=original\n')
            with patch.object(auth.subprocess, 'run') as run:
                for result in [subprocess.CompletedProcess([], 0, x, '') for x in cases] + [
                        subprocess.CompletedProcess([], 1, 'secret', 'secret')]:
                    run.return_value = result
                    with self.assertRaises(ValueError):
                        auth.configure('/opt/aws', 'us-west-2', 'secret-arn', config)
                    self.assertEqual(config.read_text(), 'ECS_CLUSTER=original\n')

    def test_cli_errors_do_not_print_secret_or_aws_response(self):
        with tempfile.TemporaryDirectory() as tmp:
            fake_aws = Path(tmp) / 'aws'
            fake_aws.write_text('#!/bin/sh\necho sentinel-secret\necho sentinel-secret >&2\nexit 1\n')
            fake_aws.chmod(0o700)
            result = subprocess.run([sys.executable, str(SCRIPT), '--aws-cli', str(fake_aws),
                                     '--region', 'us-west-2', '--secret-arn', 'test-secret',
                                     '--config', str(Path(tmp) / 'ecs.config')],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 1)
            self.assertNotIn('sentinel-secret', result.stdout + result.stderr)
            self.assertIn('ECS startup blocked', result.stderr)


if __name__ == '__main__':
    unittest.main()
