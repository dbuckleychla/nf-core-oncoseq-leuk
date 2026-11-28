import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'bin/review'))
from provenance import attach_report, completed_run, enrich, recorded_context, version_rows


class ProvenanceTests(unittest.TestCase):
    def test_attachment_preserves_html_source_and_existing_trace_observations(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data = dict(runtime_provenance=dict(parameters=dict(basecall_model='dna_synthetic@v1')),
                        run_provenance=dict(containers=[dict(process='CALL', image='synthetic-image', source='Execution trace')]))
            html = '<html><script type="application/json" id="review-data">' + json.dumps(data) + '</script><p>Preserved markup</p></html>'
            source = root / 'source.html'
            source.write_text(html)
            versions = root / 'versions.json'
            versions.write_text(json.dumps(['CALL:\n  synthetic: 1.20\n']))
            output = root / 'output.html'
            attach_report(source, output, versions)
            self.assertEqual(source.read_text(), html)
            self.assertIn('<p>Preserved markup</p>', output.read_text())
            self.assertIn('"version":"1.20"', output.read_text())
            self.assertIn('synthetic-image', output.read_text())

    def test_versions_preserve_literal_version_and_process(self):
        rows = version_rows('"WORKFLOW:CALL":\n  tool: 1.20\n  build: "2.2.0+synthetic"\n', 'Synthetic versions')
        self.assertEqual([(r['process'], r['tool'], r['version']) for r in rows], [('CALL', 'tool', '1.20'), ('CALL', 'build', '2.2.0+synthetic')])

    def test_context_only_retains_allowed_metadata(self):
        raw = dict(workflow=dict(pipeline='SYNTHETIC', unrelated='excluded'),
                   parameters=dict(basecall_model='dna_synthetic@v1', input='excluded', unrelated='excluded'),
                   software_version_texts=[])
        safe = recorded_context(raw)
        self.assertEqual(safe['parameters'], {'basecall_model': 'dna_synthetic@v1'})
        self.assertEqual(safe['workflow'], {'pipeline': 'SYNTHETIC'})

    def test_trace_containers_are_observed_not_inferred_from_lock(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            info = root / 'pipeline_info'
            info.mkdir()
            (info / 'params_2026-01-01.json').write_text(json.dumps(dict(basecall_model='dna_synthetic@v1', unrelated='excluded')))
            image = 'example.invalid/tool@sha256:' + 'a'*64
            (info / 'execution_trace_2026-01-01.txt').write_text('name\tstatus\tcontainer\nWF:CALL (SYNTHETIC)\tCOMPLETED\t' + image + '\nWF:FAILED\tFAILED\tother-image\n')
            run = completed_run(root, 'SYNTHETIC')
            self.assertEqual(run['containers'], [dict(process='CALL', image=image, source='Execution trace')])
            self.assertNotIn('unrelated', run['parameters'])
            data = dict(sample='SYNTHETIC', figures=[], references={},
                        images_provenance=dict(images={'unused': dict(image='configured-only', processes=['UNUSED'])}))
            enrich(data, root)
            self.assertEqual(data['run_provenance']['containers'], run['containers'])
            self.assertEqual(data['run_provenance']['configured_containers'][0]['image'], 'configured-only')

    def test_models_and_qc_versions_keep_sources_separate(self):
        data = dict(runtime_provenance=dict(parameters=dict(basecall_model='dna_synthetic@v1', m_bases='5mCG@v2'),
                                           software_version_texts=['CALL:\n  tool: 1.20\n']),
                    extended_qc=dict(provenance=dict(basecalling_models=['dna_observed@v3'], versions={'pysam':'0.23.3'}, modkit_version='0.6.1')))
        enrich(data)
        models = data['run_provenance']['models']
        self.assertEqual(models[0]['source'], 'BAM read-group metadata')
        self.assertEqual(models[-1]['specification'], 'dna_synthetic@v1_5mCG@v2')
        self.assertEqual(models[-1]['source'], 'Derived from configured model requests')
        self.assertIn('1.20', [r['version'] for r in data['run_provenance']['versions']])
        self.assertEqual(data['run_provenance']['containers'], [])


if __name__ == '__main__':
    unittest.main()
