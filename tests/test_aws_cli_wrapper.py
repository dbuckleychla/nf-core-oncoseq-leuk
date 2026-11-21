"""Check the Batch-mounted CLI wrapper independently of AWS credentials."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class AwsCliWrapperTests(unittest.TestCase):
    def test_wrapper_preserves_arguments_and_library_path(self):
        template = (ROOT / 'terraform/aws/templates/user_data.mime').read_text()
        wrapper = template.split("<<'AWS_WRAPPER'\n", 1)[1].split('\nAWS_WRAPPER', 1)[0]
        wrapper = wrapper.replace('$${', '${')
        with tempfile.TemporaryDirectory(prefix='cli test ') as temporary:
            root = Path(temporary)
            (root / 'bin').mkdir()
            (root / 'dist').mkdir()
            (root / 'runtime-libs').mkdir()
            executable = root / 'bin/aws'
            executable.write_text(wrapper)
            executable.chmod(0o755)
            original = root / 'dist/aws'
            original.write_text('#!/bin/sh\nprintf "%s\\n" "$LD_LIBRARY_PATH" "$@"\n')
            original.chmod(0o755)
            for previous in ('', '/existing/library/path'):
                result = subprocess.run([str(executable), 's3', 'cp', 'path with spaces'],
                                        env=dict(os.environ, LD_LIBRARY_PATH=previous),
                                        capture_output=True, text=True, check=True)
                expected = str(root / 'runtime-libs') + (':' + previous if previous else '')
                self.assertEqual(result.stdout.splitlines(),
                                 [expected, 's3', 'cp', 'path with spaces'])


if __name__ == '__main__':
    unittest.main()
