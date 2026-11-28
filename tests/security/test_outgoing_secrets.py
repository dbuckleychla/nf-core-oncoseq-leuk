"""Ensure faster routine scans retain safe full-history fallbacks."""
import importlib.util
from pathlib import Path
import unittest

path = Path(__file__).resolve().parents[2] / 'scripts/check_outgoing_secrets.py'
spec = importlib.util.spec_from_file_location('outgoing_secrets', path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class OutgoingSecretTests(unittest.TestCase):
    def test_known_upstream_scans_all_outgoing_commits(self):
        self.assertEqual(module.scan_range('b' * 40, 'a' * 40), 'a' * 40 + '..' + 'b' * 40)

    def test_unknown_upstream_audits_all_history(self):
        self.assertEqual(module.scan_range('b' * 40), '--all')

    def test_explicit_audit_includes_all_references(self):
        self.assertEqual(module.scan_range('b' * 40, 'a' * 40, True), '--all')


if __name__ == '__main__':
    unittest.main()
