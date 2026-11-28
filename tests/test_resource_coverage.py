"""Prevent new BAM/read modules and aliases from bypassing isolated workers."""
from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[1]


class ResourceCoverageTests(unittest.TestCase):
    def test_every_bam_or_read_input_and_import_alias_has_isolated_routing(self):
        selectors = re.findall(r"withName:\s*'([^']+)'", (ROOT / 'conf/profiles/aws_bam.config').read_text())
        needed = {'SAMTOOLS_INDEX', 'SAMTOOLS_MERGE', 'SAMTOOLS_VIEW', 'SNIFFLES', 'ENSEMBLVEP_VEP'}
        for path in (ROOT / 'modules').glob('**/main.nf'):
            for name, body in re.findall(r'^process\s+(\w+)\s*\{([\s\S]*?)(?=^process\s|\Z)', path.read_text(), re.M):
                inputs = body.split('input:', 1)[-1].split('output:', 1)[0] if 'input:' in body else ''
                if re.search(r'\bpath\(?\s*(?:bam|bams|ubam|in_bam|reads|fastq)\b', inputs):
                    needed.add(name)
        for path in list((ROOT / 'subworkflows/local').glob('**/*.nf')) + list((ROOT / 'workflows').glob('*.nf')):
            for original, alias in re.findall(r'include\s*\{\s*(\w+)\s+as\s+(\w+)\s*\}', path.read_text()):
                if original in needed:
                    needed.add(alias)
        needed.remove('DORADO_BASECALL')
        missing = sorted(name for name in needed if not any(re.fullmatch(selector, name) for selector in selectors))
        self.assertEqual(missing, [])


if __name__ == '__main__':
    unittest.main()
