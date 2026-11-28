"""Check that public-content rejection does not print or store real secrets."""
import gzip
import importlib.util
from pathlib import Path
import unittest

source = Path(__file__).resolve().parents[2] / 'scripts/check_public_content.py'
spec = importlib.util.spec_from_file_location('public_content', source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class PublicContentTests(unittest.TestCase):
    def test_reference_paths_and_synthetic_names_are_allowed(self):
        self.assertEqual(module.findings(b's3://reference-bucket/hg38/vep-cache/ SAMPLE_001 SYNTHETIC'), [])

    def test_accessions_are_rejected_in_plain_and_compressed_content(self):
        accession = b'M' + b'99' + b'-' + b'123456'
        self.assertIn('clinical accession', module.findings(accession))
        self.assertIn('clinical accession', module.findings(gzip.compress(accession)))
        self.assertIn('clinical accession', module.findings(accession.replace(b'-', b'%2D')))

    def test_generated_fake_tokens_are_rejected(self):
        examples = [(b'A' + b'KIA' + b'Z' * 16, 'AWS access key'),
                    (b'gh' + b'p_' + b'z' * 36, 'GitHub token'),
                    (b'gl' + b'pat-' + b'z' * 24, 'GitLab token')]
        for content, category in examples:
            self.assertIn(category, module.findings(content))

    def test_unreadable_compressed_content_is_rejected(self):
        with self.assertRaises(ValueError):
            module.findings(b'\x1f\x8b' + b'incomplete')


if __name__ == '__main__':
    unittest.main()
