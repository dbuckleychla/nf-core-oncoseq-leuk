"""Synthetic regression contracts for the additive blood classifier reports."""
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'bin/review'))
from adapters import classifiers
from presentation import add_classy_summary, enrich, classy_publication
from provenance import enrich as provenance

spec = importlib.util.spec_from_file_location('blood_tables', ROOT / 'bin/classy_blood_tables.py')
tables = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tables)


def synthetic():
    data = dict(sample_label='SYNTHETIC', oncoseq_classifier=dict(workflow='classy_blood', enable_alma3=True))
    data['combined'] = dict(classification_scores=[])
    for key in ('marlin', 'lamprey', 'alma3', 'optic'):
        data[key] = dict(model_version='synthetic-v1', model_version_source='Synthetic fixture',
                         target_coverage=dict(observed_target_count=6000))
        data['combined']['classification_scores'].append(dict(classifier_run_id=key + ':synthetic',
            native_path=[dict(name='Synthetic class')], score=0.9, cutoff=0.8,
            accepted_call=key != 'optic'))
    return data


class BloodReportTests(unittest.TestCase):
    def test_alma_uses_native_resolved_path_instead_of_comparing_head_scores(self):
        raw = synthetic()
        raw['alma3']['result'] = dict(result_summary='Synthetic subtype accepted', path=[
            dict(level='presence', status='resolved', classification='Synthetic present'),
            dict(level='subtype', status='resolved', classification='Synthetic subtype')])
        raw['combined']['classification_scores'] += [
            dict(classifier_run_id='alma3:presence', native_label='Synthetic present', score=0.99, accepted_call=True),
            dict(classifier_run_id='alma3:subtype', native_label='Synthetic subtype', score=0.85, accepted_call=True)]
        alma = next(c for c in classifiers(raw) if c['name'] == 'alma3')
        self.assertEqual(alma['label'], 'Synthetic subtype accepted')
        self.assertEqual(alma['score'], 0.85)

    def test_rollback_selects_legacy_even_when_new_outputs_are_stale(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            legacy = root / 'reports/classy_legacy/SYNTHETIC'
            legacy.mkdir(parents=True)
            (legacy / 'SYNTHETIC_combined_classification.json').write_text('{}')
            self.assertEqual(classy_publication(root, 'SYNTHETIC'), legacy)
            (root / 'pipeline_info').mkdir()
            (root / 'pipeline_info/params_synthetic.json').write_text('{"classy_blood":false}')
            new = root / 'reports/tumor_classifiers'
            new.mkdir()
            (new / 'SYNTHETIC_combined_classification.json').write_text('{}')
            self.assertEqual(classy_publication(root, 'SYNTHETIC'), legacy)
    def test_native_decision_and_approximation_are_preserved(self):
        cards = {c['name']: c for c in classifiers(synthetic())}
        self.assertFalse(cards['optic']['accepted_call'])
        self.assertGreater(cards['optic']['score'], cards['optic']['cutoff'])
        self.assertIn('native approximation', cards['alma3']['display_name'])
        self.assertTrue(cards['alma3']['approximation'])
        data = dict(classification=list(cards.values()))
        provenance(data)
        self.assertEqual(len(data['run_provenance']['models']), 4)

    def test_disabled_missing_and_zero_observations(self):
        raw = synthetic()
        raw['oncoseq_classifier']['enable_alma3'] = False
        raw['lamprey']['target_coverage']['observed_target_count'] = 0
        raw['combined']['classification_scores'] = [r for r in raw['combined']['classification_scores'] if not r['classifier_run_id'].startswith('optic')]
        cards = {c['name']: c for c in classifiers(raw)}
        self.assertEqual(cards['alma3']['status'], 'disabled')
        self.assertEqual(cards['optic']['status'], 'unavailable')
        self.assertIsNone(cards['lamprey']['score'])
        self.assertIsNone(cards['lamprey']['accepted_call'])

    def test_native_json_unchanged_and_reports_use_new_plots(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            published = root / 'reports/tumor_classifiers'
            published.mkdir(parents=True)
            native = published / 'SYNTHETIC_combined_classification.json'
            raw = synthetic()
            del raw['oncoseq_classifier']
            original = json.dumps(raw)
            native.write_text(original)
            plot = 'SYNTHETIC_combined_classification_combined_blood.svg'
            (published / plot).write_text('<svg xmlns="http://www.w3.org/2000/svg"/>')
            output = root / 'metrics'
            tables.prepare(native, 'SYNTHETIC', published, output, True, '5mc', 'synthetic@sha256:synthetic')
            self.assertEqual(native.read_text(), original)
            self.assertEqual(len(list((output / 'tables').glob('*.tsv'))), 4)
            self.assertIn('Not accepted', (output / 'tables/optic.tsv').read_text())
            qmd = (output / 'native_report_inputs/SYNTHETIC-Methylation.qmd').read_text()
            self.assertIn('Analysis performed using ALMA3.', qmd)
            self.assertIn('numerical parity', qmd)
            self.assertIn(plot, qmd)
            metadata = (output / 'classification.metadata.json').read_text()
            (published / 'SYNTHETIC_classification.metadata.json').write_text(metadata)
            # Stale legacy outputs from a resumed outdir must not enter the report.
            (published / 'SYNTHETIC_combined_classification_combined_brain.svg').write_text('<svg/>')
            data = dict(sample='SYNTHETIC', figures=[dict(category='classification', path='old_combined_brain.svg')])
            enrich(data, root / 'private-review', root)
            self.assertEqual(data['classifier_workflow'], 'classy_blood')
            self.assertEqual({f['classy_plot'] for f in data['figures']}, {'blood'})
            add_classy_summary(data, root)
            self.assertEqual(len(data['figures']), 1)
            self.assertEqual({c['name'] for c in data['classification']}, {'marlin', 'lamprey', 'alma3', 'optic'})
            with self.assertRaisesRegex(ValueError, 'sample mismatch'):
                tables.prepare(native, 'SYNTHETIC_OTHER', published, output, True, '5mc', 'synthetic')
