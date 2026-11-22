import sys
import tempfile
import unittest
from pathlib import Path
import xml.etree.ElementTree as ET
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'bin/review'))
from adapters import classifiers
from presentation import enrich, genome_plot


class PresentationTests(unittest.TestCase):
    def test_native_score_shapes_and_cutoff_roles(self):
        data = {'sturgeon': {'models': [{'model_root': '/models/general', 'top_class': {'class': 'Synthetic class', 'score': .9}, 'scores': [{'class': 'Synthetic class', 'score': .9}], 'threshold': .5, 'cutoff': .8}]},
                'tucan': {'top_class': {'name': 'Synthetic tissue', 'score': .6}, 'averaged_scores': [{'name': 'Synthetic tissue', 'score': .6}], 'cutoff': .7},
                'nanomix': {'scores': [{'tissue': 'Synthetic tissue', 'proportion': .4}], 'cutoff': .2}}
        result = classifiers(data)
        sturgeon = next(r for r in result if r['name'].startswith('sturgeon'))
        self.assertEqual((sturgeon['label'], sturgeon['score'], sturgeon['cutoff']), ('Synthetic class', .9, .8))
        self.assertEqual(next(r for r in result if r['name'] == 'nanomix')['cutoff'], None)
        self.assertEqual(next(r for r in result if r['name'] == 'lamprey')['status'], 'unavailable')

    def test_weighted_coverage_and_native_plot(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data = {'coverage': [dict(gene='DEMO', start=0, end=10, primary=10), dict(gene='DEMO', start=10, end=40, primary=30)]}
            enrich(data, root)
            self.assertEqual(data['coverage_genes'][0]['depth'], 25)
            cov = root/'coverage.tsv'; maf = root/'maf.tsv'; plot=root/'plot.svg'
            cov.write_text('chromosome\tbin_start\tbin_end\tn_reads\tn_reads_gc_adj\nchr2\t0\t1000000\t10\t12\nchr1\t0\t1000000\t8\t9\n')
            maf.write_text('chrom\tpos\tref\talt\nchr1\t100\t15\t5\nchr2\t200\t1\t1\n')
            genome_plot(cov, maf, plot)
            ET.fromstring(plot.read_text())
            self.assertIn('depth ≥20', plot.read_text())
            self.assertIn('not base depth', plot.read_text())
