import sys
import tempfile
import unittest
from pathlib import Path
import xml.etree.ElementTree as ET
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'bin/review'))
from adapters import classifiers
from presentation import add_classy_summary, enrich, genome_plot


class PresentationTests(unittest.TestCase):
    def test_classy_summary_from_copied_pipeline_artifact(self):
        data = dict(sample='SYNTHETIC', figures=[dict(category='classification', path='figures/0012_SYNTHETIC_combined_classification_combined_task_pies.svg')])
        add_classy_summary(data)
        add_classy_summary(data)
        self.assertEqual(len(data['figures']), 1)
        self.assertTrue(data['figures'][0]['classy_summary'])

    def test_classy_summary_discovery_is_sample_specific_and_optional(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            folder = root / 'reports/tumor_classifiers'
            folder.mkdir(parents=True)
            (folder / 'SYNTHETIC_OTHER_combined_classification_combined_task_pies.svg').write_text('<svg/>')
            data = dict(sample='SYNTHETIC')
            add_classy_summary(data, root)
            self.assertEqual(data['figures'], [])
            path = folder / 'SYNTHETIC_combined_classification_combined_task_pies.svg'
            path.write_text('<svg/>')
            add_classy_summary(data, root)
            add_classy_summary(data, root)
            self.assertEqual(len(data['figures']), 1)
            self.assertTrue(data['figures'][0]['classy_summary'])
            self.assertEqual(Path(data['figures'][0]['path']), path.resolve())

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
            maf.write_text('chrom\tpos\tref\talt\n1\t100\t15\t5\n1\t700000\t5\t15\n2\t200\t1\t1\n')
            baf = genome_plot(cov, maf, plot)
            ET.fromstring(plot.read_text())
            self.assertIn('depth ≥20', baf.read_text())
            self.assertIn('not log ratio', plot.read_text())
            self.assertIn('1.00', baf.read_text())
            cells = [e for e in ET.fromstring(baf.read_text()).iter() if e.tag.endswith('rect') and e.get('fill') == '#e9d5ff']
            self.assertEqual(len(cells), 2)
            self.assertLess(min(float(e.get('y')) for e in cells), 205)
            self.assertGreater(max(float(e.get('y')) for e in cells), 205)
            self.assertTrue(data['qc_plots'])
