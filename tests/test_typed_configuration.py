"""Exercise the real typed parameter block without launching analysis tasks."""
import os
import re
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which('nextflow')


@unittest.skipUnless(NEXTFLOW, 'Nextflow is required for typed configuration checks')
class TypedConfigurationTests(unittest.TestCase):
    def check_initialisation(self, profile, high='', low=''):
        environment = dict(os.environ, NXF_OFFLINE='true',
                           ONCOSEQ_AWS_BAM_HIGH_VCPUS=high,
                           ONCOSEQ_AWS_BAM_LOW_VCPUS=low)
        path_arguments = []
        for name in re.findall(r'^\s+(\w+): Path\??', (ROOT / 'main.nf').read_text(), re.M):
            path_arguments.extend(['--' + name, str(ROOT / 'assets/NOFILE')])
        configuration = ['-profile', profile] if profile == 'aws,clinical' else ['-c', str(ROOT / 'conf/profiles' / (profile + '.config'))]
        with tempfile.TemporaryDirectory() as temporary:
            result = subprocess.run(
                [NEXTFLOW, '-log', str(Path(temporary) / 'nextflow.log'), 'run', str(ROOT),
                 *configuration, '--version', *path_arguments,
                 '--outdir', temporary, '--basecall_model', 'dna_r10.4.1_e8.2_400bps_sup@v5.2.0'],
                cwd=temporary, env=environment, capture_output=True, text=True, timeout=120)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('nf-core/oncoseq v', result.stdout)

    def test_aws_exported_integer_strings(self):
        self.check_initialisation('aws,clinical', '32', '8')

    def test_aws_empty_exports_use_defaults(self):
        self.check_initialisation('aws,clinical')

    def test_profile_resource_types(self):
        for profile in ('test_fastq', 'test_drac', 'test_cfdna', 'test_full_wgs', 'test_full_time_series'):
            with self.subTest(profile=profile):
                self.check_initialisation(profile)


if __name__ == '__main__':
    unittest.main()
