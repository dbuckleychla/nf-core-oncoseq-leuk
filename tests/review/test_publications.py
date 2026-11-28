"""Offline citation catalog integrity and reviewed publication replacements."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'bin/review'))
from publications import catalog


class PublicationTests(unittest.TestCase):
    def test_unique_complete_entries_and_no_placeholder_dois(self):
        data = catalog()
        self.assertEqual(data['schema_version'], 1)
        self.assertEqual(len({r['id'] for r in data['entries']}), len(data['entries']))
        for row in data['entries']:
            self.assertTrue(row['citation'])
            self.assertTrue(row['group'])
            self.assertTrue(row['software_url'].startswith('https://'))
            self.assertNotIn('XXXX', row['citation'])
            if row['doi']:
                self.assertRegex(row['doi'], r'^10\.\d{4,9}/\S+$')
                self.assertNotIn('XXXX', row['doi'])
            if row['kind'] == 'publication':
                self.assertIsNotNone(row['doi'])
                self.assertFalse(row['doi'].startswith('10.1101/'))

    def test_reviewed_workflow_and_tool_publications(self):
        rows = {r['name']: r for r in catalog()['entries']}
        self.assertEqual(rows['nf-core/oncoseq']['doi'], '10.1186/s13073-026-01815-7')
        self.assertEqual(rows['DeepSomatic']['doi'], '10.1038/s41587-025-02839-x')
        self.assertEqual(rows['SnpEff']['doi'], '10.4161/fly.19695')
        self.assertEqual(rows['ClairS-TO']['kind'], 'preprint')
        self.assertEqual(rows['Dorado']['kind'], 'software')
        self.assertIsNone(rows['Dorado']['doi'])
        self.assertNotIn('FrankONTstein', rows)
        self.assertEqual(catalog()['entries'][0]['name'], 'nf-core/oncoseq')

    def test_catalog_is_independent_between_reports(self):
        a = catalog()
        a['entries'][0]['citation'] = 'Changed for a synthetic test'
        a['sources'].clear()
        self.assertNotEqual(catalog()['entries'][0]['citation'], a['entries'][0]['citation'])
        self.assertTrue(catalog()['sources'])


if __name__ == '__main__':
    unittest.main()
