import json
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'bin/review'))
from alignment_qc import load_alignment, primary_counts, apply_alignment, patch_html


class AlignmentQCTests(unittest.TestCase):
    def test_secondary_and_supplementary_do_not_enter_denominator(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'qc.json'
            path.write_text(json.dumps({'QC-passed reads':{'primary':8,'primary mapped':6,'secondary':20,'supplementary':10},'QC-failed reads':{'primary':2,'primary mapped':1}}))
            qc=load_alignment(path)
            self.assertEqual((qc['primary_reads'],qc['mapped_primary_reads'],qc['unmapped_primary_reads'],qc['alignment_rate_percent']),(10,7,3,70))

    def test_zero_total_is_unavailable_and_updates_are_idempotent(self):
        data={'qc_summary':[{'metric':'Other','value':1}], 'qc_plots':[]}
        qc=primary_counts(0,0,'Synthetic')
        apply_alignment(data,qc);apply_alignment(data,qc)
        self.assertEqual(data['qc_summary'][0]['value'],'Not available')
        self.assertEqual(len(data['qc_summary']),5)
        self.assertEqual(len(data['qc_plots']),1)

    def test_html_patch_preserves_images_and_non_payload_text(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); source=root/'report.html';target=root/'updated.html';qc=root/'qc.json'
            data={'sample':'SYNTHETIC','images':{'synthetic':'data:image/png;base64,AAAA'},'qc_summary':[]}
            source.write_text('<p>Before</p><script type="application/json" id="review-data">'+json.dumps(data)+'</script><p>After</p>')
            qc.write_text(json.dumps(primary_counts(9,1,'Synthetic')))
            patch_html(source,target,qc)
            self.assertTrue(target.read_text().startswith('<p>Before</p>'))
            self.assertTrue(target.read_text().endswith('<p>After</p>'))
            self.assertIn('data:image/png;base64,AAAA',target.read_text())
            self.assertIn('90.0',target.read_text())
