import gzip
import json
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest

import jsonschema
import pysam

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'bin'))
sys.path.insert(0, str(Path(__file__).parent))
from build_fixture import build
from qc_extended import alignment, methylation, regions
from qc_extended.assemble import assemble
from qc_extended.contracts import required_depth, state, write_tables
from qc_extended.report import content


class QCTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='synthetic-qc-')
        self.root = build(self.tmp.name)
        regions.reference(str(self.root / 'reference.fa'), self.root / 'reference.json', self.root / 'cpgs.gz')
        self.ref = json.loads((self.root / 'reference.json').read_text())

    def tearDown(self):
        self.tmp.cleanup()

    def collect(self, variant='plain'):
        return alignment.collect(str(self.root / (variant + '.bam')), self.ref,
                                 self.root / 'enrichment.bed', self.root / 'targets.bed', 2)

    def test_frankontstein_parity(self):
        a = self.collect()
        self.assertEqual(a['counts']['primary_records'], 4)
        self.assertEqual(a['counts']['mapped_primary'], 3)
        self.assertEqual(a['read_lengths']['all']['bases'], 45)
        self.assertEqual(a['read_lengths']['all']['n50'], 10)
        self.assertEqual(a['read_lengths']['all']['median'], 10)
        self.assertEqual(a['read_lengths']['on_target']['records'], 1)
        self.assertEqual(a['read_lengths']['on_enrichment']['records'], 2)
        self.assertEqual(a['coverage']['0']['enrichment']['bases'], 15)
        self.assertEqual(a['coverage']['0']['genome_primary']['aligned_bases'], 40)
        self.assertEqual(a['coverage']['20']['genome_primary']['aligned_bases'], 30)
        self.assertEqual(a['coverage']['20']['enrichment']['aligned_bases'], 10)

    def test_overlap_excludes_deletions_skips_and_clips(self):
        r = pysam.AlignedSegment()
        r.query_sequence, r.reference_id, r.reference_start = 'C' * 20, 0, 0
        for cigar in ('10M10D10M', '10M10N10M'):
            r.cigarstring = cigar
            self.assertIn('off_target', alignment.read_groups(r, regions.Regions([(0, 100)]), regions.Regions([(10, 20)])))
        r.flag = 256
        self.assertEqual(alignment.read_groups(r, regions.Regions([]), regions.Regions([])), [])

    def test_missing_length_and_empty_denominators(self):
        d = alignment.length_summary({}, 2)
        self.assertEqual(d['records'], 2)
        self.assertIsNone(d['n50'])
        self.assertEqual(alignment.depth_summary({}, 100)['histogram'], {0: 100})
        self.assertIsNone(alignment.depth_summary({}, 0)['mean'])
        self.assertEqual(alignment.length_summary({10: 1, 20: 1})['median'], 15)

    def test_empty_bam_has_explicit_zero_counts(self):
        path = self.root / 'empty.bam'
        with pysam.AlignmentFile(str(path), 'wb', header={'HD': {'SO': 'coordinate'}, 'SQ': [{'SN': 'chr1', 'LN': 100}]}):
            pass
        pysam.index(str(path))
        result = self.collect('empty')
        self.assertEqual(result['counts']['primary_records'], 0)
        self.assertEqual(result['coverage']['0']['genome_primary']['histogram'], {0: 100})
        self.assertIsNone(result['read_lengths']['all']['n50'])

    def test_modification_integrity(self):
        tagged = self.collect('tagged')
        self.assertEqual(tagged['modifications']['usable'], 3)
        self.assertEqual(tagged['modifications']['mapped_invalid'], 0)
        invalid = self.collect('invalid')
        self.assertEqual(invalid['modifications']['mapped_invalid'], 3)
        self.assertEqual(self.collect()['modifications']['usable'], 0)

    def test_methylation_shared_depth_and_availability(self):
        m = methylation.collect(self.root / 'mods.bed.gz', self.ref, self.root / 'cpgs.gz',
                                self.root / 'enrichment.bed', self.root / 'targets.bed', ['m', 'h'])
        d = m['scopes']['genome_primary']['combined']
        self.assertEqual((d['sites'], d['valid'], d['modified'], d['reference_cpgs']), (1, 10, 5, 50))
        self.assertEqual(m['scopes']['genome_primary']['5mC']['modified'], 3)
        only = methylation.collect(self.root / 'mods.bed.gz', self.ref, self.root / 'cpgs.gz',
                                   self.root / 'enrichment.bed', self.root / 'targets.bed', ['m'])
        self.assertFalse(only['scopes']['genome_primary']['5hmC']['available'])

    def test_serialization_reports_and_no_private_records(self):
        data = assemble('SYNTHETIC_A', json.loads(json.dumps(self.collect())), dict(status=state('unavailable', 'no_usable_modifications')), self.ref)
        schema = json.loads((ROOT / 'assets/qc_extended/metrics.schema.json').read_text())
        jsonschema.validate(data, schema)
        write_tables(self.root, data)
        text = content(data)
        self.assertIn('Read lengths on/off exact targets', text)
        self.assertIn('Theoretical low-VAF', text)
        self.assertIn('no_usable_modifications', text)
        self.assertNotIn('SYNTHETIC_SUPP', text)
        self.assertNotIn('SYNTHETIC_SUPP', json.dumps(data))
        self.assertTrue((self.root / 'adequacy.tsv').exists())
        zero = next(r for r in data['metrics'] if r['metric_id'] == 'coverage.zero_fraction' and r['scope'] == 'genome_primary' and r['mapq'] == 0)
        self.assertEqual(zero['value'], .6)

    def test_binomial_oracles(self):
        for vaf, expected in [(.01, (628, 913)), (.02, (313, 456)), (.05, (124, 181)), (.1, (61, 89))]:
            self.assertEqual(tuple(required_depth(vaf, k) for k in (3, 5)), expected)

    def test_overlapping_targets_match_bruteforce(self):
        (self.root / 'targets.bed').write_text('chr1\t0\t10\ta\nchr1\t5\t15\tb\nchr1\t10\t20\tc\nchr1\t10\t20\tduplicate\nchr1\t60\t70\tempty\n')
        a = self.collect()
        depths = [0] * 100
        for lo, hi in [(0, 10), (10, 20), (20, 30), (35, 45)]:
            for pos in range(lo, hi):
                depths[pos] += 1
        for row in a['target_intervals']:
            d = a['contigs'][0]['coverage']['0']['targets'][row['target_id']]
            self.assertEqual(d['aligned_bases'], sum(depths[row['start']:row['end']]))

    def test_reference_boundaries_and_validation(self):
        self.assertEqual(self.ref['cpg_counts']['chr1'], 50)
        (self.root / 'targets.bed').write_text('chr1\t99\t101\tbad\n')
        with self.assertRaises(ValueError):
            regions.read_bed(self.root / 'targets.bed', self.ref['lengths'])

    def test_malformed_methylation_denominator_is_rejected(self):
        with gzip.open(self.root / 'bad.bed.gz', 'wt') as out:
            out.write('chr1\t0\t1\tm\t10\t.\t0\t1\t0\t10\t30\t3\t9\t0\t0\t0\t0\t0\n')
        with self.assertRaises(ValueError):
            methylation.collect(self.root / 'bad.bed.gz', self.ref, self.root / 'cpgs.gz',
                                self.root / 'enrichment.bed', self.root / 'targets.bed', ['m'])

    def test_feature_dedup_and_reference_validation(self):
        (self.root / 'features.bed').write_text('chr1\t0\t1\tfirst\nchr1\t0\t1\tduplicate\n')
        (self.root / 'panels.tsv').write_text('panel_id\tassembly\tmodel_name\tmodel_version\tbed\nSYNTHETIC_PANEL\thg38\tSYNTHETIC_MODEL\t1\tfeatures.bed\n')
        p = methylation.panels(self.root / 'panels.tsv', str(self.root / 'reference.fa'), self.ref['lengths'], 'hg38')
        self.assertEqual(len(p[0]['features']), 1)
        self.assertEqual(p[0]['input_rows'], 2)
        with self.assertRaises(ValueError):
            methylation.panels(self.root / 'panels.tsv', str(self.root / 'reference.fa'), self.ref['lengths'], 'hs1')

    def test_portable_review_render_embeds_qc(self):
        sys.path.insert(0, str(ROOT / 'bin/review'))
        from render import render
        qc = assemble('SYNTHETIC_A', self.collect(), dict(status=state('disabled', 'user_disabled')), self.ref)
        (self.root / 'qc.extended.json').write_text(json.dumps(qc))
        data = dict(schema_version=1, sample='SYNTHETIC_A', assembly='GRCh38',
                    small_variants={'records': [], 'total_records': 0}, structural_events=[],
                    classification=[], coverage=[], cnv=[], nasvar={}, figures=[], references={},
                    panel={}, warnings=[], statuses={})
        (self.root / 'review.json').write_text(json.dumps(data))
        render(self.root / 'review.json', self.root / 'report.html', ROOT / 'assets/review/report.html')
        self.assertIn('extended_qc_html', (self.root / 'report.html').read_text())
        self.assertIn('Read lengths on/off enrichment', (self.root / 'report.html').read_text())

    def test_attach_native_report_preserves_source_and_navigation(self):
        source = self.root / "original's native.html"
        original = '<html><nav id="TOC"><ul><li>Coverage</li></ul></nav><main><section id="coverage">SYNTHETIC</section></main></html>'
        source.write_text(original)
        qc = assemble('SYNTHETIC_A', self.collect(), dict(status=state('disabled', 'user_disabled')), self.ref)
        (self.root / 'fragment.html').write_text(content(qc))
        result = self.root / 'final.html'
        subprocess.run([sys.executable, str(ROOT / 'bin/extended_qc.py'), 'attach-report',
                        '--input', str(source), '--fragment', str(self.root / 'fragment.html'), '--output', str(result)], check=True)
        self.assertEqual(source.read_text(), original)
        self.assertIn('ONCOSEQ_EXTENDED_QC', result.read_text())
        self.assertIn('href="#extended-qc"', result.read_text())
        self.assertIn('Read lengths on/off exact targets', result.read_text())


if __name__ == '__main__':
    unittest.main()
