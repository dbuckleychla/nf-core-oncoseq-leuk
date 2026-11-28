"""Shell contract checks; these do not validate the actual modkit image."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]


class ModkitWrapperTests(unittest.TestCase):
    def run_wrapper(self, exit_code):
        source = (ROOT / 'modules/local/extended_qc/main.nf').read_text()
        process = source.split('process QC_EXT_MODKIT_PILEUP {', 1)[1].split('process QC_EXT_METHYLATION {', 1)[0]
        script = process.split('"""', 2)[1]
        script = script.replace('${qcQuote(bam)}', "'synthetic input.bam'")
        script = script.replace('${qcQuote(fasta)}', "'synthetic reference.fa'")
        script = script.replace('${Math.min(8, task.cpus as int)}', '8').replace('\\\\\n', '\\\n')
        with tempfile.TemporaryDirectory(prefix='synthetic-modkit-wrapper-') as directory:
            root = Path(directory)
            tool = root / 'modkit'
            tool.write_text("#!/bin/sh\n"
                            "if [ \"$1\" = --version ]; then printf 'modkit synthetic-test\\n'; exit 0; fi\n"
                            "printf '%s\\n' \"$@\" > arguments.txt\n"
                            f"exit {exit_code}\n")
            tool.chmod(0o755)
            # PATH contains only the stand-in modkit; no Python is available.
            result = subprocess.run(['/bin/bash', '-euo', 'pipefail', '-c', script], cwd=root,
                                    env={**os.environ, 'PATH': str(root)}, capture_output=True, text=True)
            args = (root / 'arguments.txt').read_text().splitlines()
            version = (root / 'modkit.version.txt').read_text() if (root / 'modkit.version.txt').exists() else None
            return result, args, version

    def test_runs_without_python_and_preserves_arguments(self):
        result, args, version = self.run_wrapper(0)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(args[:5], ['pileup', 'synthetic input.bam', 'pileup.bed.gz', '--ref', 'synthetic reference.fa'])
        self.assertEqual(args[args.index('--modified-bases') + 1:args.index('--modified-bases') + 3], ['5mC', '5hmC'])
        self.assertIn('--bgzf', args)
        self.assertEqual(version, 'modkit synthetic-test\n')

    def test_pileup_failure_stops_before_version_collection(self):
        result, _args, version = self.run_wrapper(42)
        self.assertEqual(result.returncode, 42)
        self.assertIsNone(version)


if __name__ == '__main__':
    unittest.main()
