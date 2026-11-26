"""Presentation arithmetic must preserve the aggregate QC evidence."""
import copy
import math
from pathlib import Path
import sys
import unittest
import xml.etree.ElementTree as ET

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'bin'))
from qc_extended.report import bars, content, percent, rebin


class ReportTests(unittest.TestCase):
    def test_sparse_length_bins_preserve_counts_and_overflow(self):
        h = dict(denominator=20, bins=[
            dict(lower=0, upper=100, count=4),
            dict(lower=400, upper=500, count=5),
            dict(lower=500, upper=600, count=3),
            dict(lower=10000, upper=10100, count=7),
            dict(lower=144000, upper=144100, count=1),
        ])
        counts, overflow = rebin(h, list(range(0, 10001, 500)))
        self.assertEqual(counts[:3], [9, 3, 0])
        self.assertEqual(overflow, 8)
        self.assertEqual(sum(counts) + overflow, h['denominator'])
        all_counts, overflow = rebin(h, [0, 500, 1000, 2000, 5000, 10000, 20000, 50000, math.inf])
        self.assertEqual(all_counts, [9, 3, 0, 0, 0, 7, 0, 1])
        self.assertEqual(overflow, 0)

    def test_beta_upper_endpoint_and_floating_bin_edges(self):
        h = dict(bins=[dict(lower=.29, upper=.30, count=2),
                       dict(lower=.99, upper=1., upper_inclusive=True, count=3)])
        counts, overflow = rebin(h, [i/20 for i in range(21)])
        self.assertEqual(counts[5], 2)
        self.assertEqual(counts[-1], 3)
        self.assertEqual(sum(counts), 5)
        self.assertEqual(overflow, 0)

    def test_cannot_split_original_bins(self):
        with self.assertRaisesRegex(ValueError, 'splits'):
            rebin(dict(bins=[dict(lower=200, upper=300, count=10)]), [0, 250, 500])

    def test_zoom_keeps_whole_group_denominator(self):
        svg = bars('Synthetic group', [('0–500 bp', 2)], 10, [(0, '0'), (1, '500')], '#238b9a', .5, 'Length')
        tree = ET.fromstring(svg)
        self.assertIn('20.00% of group', tree.find('rect/title').text)
        self.assertEqual(float(tree.find('rect').attrib['height']), 60)
        self.assertNotIn('polyline', svg)
        small = bars('Rare tail', [('Tail <end>', 1)], 100000, [], '#238b9a', .5, 'Length')
        self.assertIn('<0.01%', ET.fromstring(small).find('rect/title').text)
        self.assertEqual(percent(0), '0.00%')
        self.assertEqual(percent(None), 'Unavailable')

    def test_empty_report_has_explicit_missing_states_and_escapes_input(self):
        data = dict(schema_name='oncoseq_extended_qc', schema_version='1.0.0',
                    sample_id='<script>synthetic</script>', metrics=[], histograms=[], warnings=[],
                    components={'methylation': dict(status='unavailable', reason='no_usable_modifications')})
        original = copy.deepcopy(data)
        html = content(data)
        self.assertEqual(data, original)
        self.assertNotIn('<script>', html)
        self.assertIn('&lt;script&gt;synthetic&lt;/script&gt;', html)
        self.assertIn('No observations available', html)
        self.assertIn('Unavailable', html)
        self.assertIn('no usable modifications', html)
        self.assertNotIn('polyline', html)


if __name__ == '__main__':
    unittest.main()
