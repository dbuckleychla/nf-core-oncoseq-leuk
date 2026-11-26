"""Scientific and report contract regressions using synthetic inputs."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'bin/review'))

from fixtures import make
from engine import prepare, complete
from report_prioritization import small_candidate, sample_evidence
from structural_consensus import consolidate, junction
from render import render
from annotate_structural_variants import endpoint
from resources import validate, prepare_vcf
from indexed_reference import IndexedReference
from structural_consensus import compatibility_edges, compatible
from adapters import nasvar_events
from structural_consensus import support
from annotate_structural_variants import annotate as annotate_sv
import pysam
import pysam.bcftools
import sqlite3


class ReviewTests(unittest.TestCase):
    def test_cached_lookups_match_direct_queries_and_respect_budget(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = make(Path(tmp) / 'refs')
            path = root / 'cosmic.vcf.gz'
            for budget in (16 * 1024 * 1024, 1):
                cache = IndexedReference(path, window=128, budget_bytes=budget, max_windows=1)
                with pysam.VariantFile(str(path)) as direct:
                    queries = [(149, 'A', 'C'), (149, 'A', 'C'), (199, 'A', 'T'),
                               (149, 'A', 'C'), (899, 'A', 'G'), (149, 'A', 'G')]
                    for start, ref, alt in queries:
                        query = direct.header.new_record(contig='chr1', start=start, alleles=(ref, alt))
                        expected = [str(r) for r in direct.fetch('chr1', start, start + 1)
                                    if r.pos == query.pos and r.ref == ref and r.alts == (alt,)]
                        self.assertEqual([str(r) for r in cache.exact(query)], expected)
                        self.assertLessEqual(cache.cached_bytes, budget)
                    self.assertEqual(cache.exact(direct.header.new_record(contig='chr2', start=9, alleles=('A', 'C'))), [])
                if budget == 1:
                    self.assertGreater(cache.oversized_windows, 0)
                cache.close()

    def test_indexed_junction_graph_matches_all_pairs(self):
        import random
        rng = random.Random(217)
        rows = []
        for i in range(120):
            rows.append(dict(id='SYNTHETIC_' + str(i), sample='SYNTHETIC',
                             caller=rng.choice(['sniffles', 'severus', 'NASVAR']),
                             chrom=rng.choice(['chr1', 'chr2']), pos=rng.randrange(100, 600),
                             partner=dict(chrom=rng.choice(['chr1', 'chr2']), pos=rng.randrange(100, 600),
                                          local_side=rng.choice(['left', 'right', None]),
                                          remote_side=rng.choice(['left', 'right', None])),
                             evidence=dict(info=dict(SVTYPE=rng.choice(['BND', 'TRA', 'DEL', 'DUP:TANDEM'])), format={})))
        for distance in (0, 50, 100):
            expected = {i: set() for i in range(len(rows))}
            for i, row in enumerate(rows):
                for j in range(i + 1, len(rows)):
                    if compatible(row, rows[j], distance):
                        expected[i].add(j)
                        expected[j].add(i)
            self.assertEqual(compatibility_edges(rows, distance), expected)

    def test_symbolic_sv_without_strands_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = make(Path(tmp) / 'refs')
            source = root / 'symbolic.vcf'
            header = pysam.VariantHeader()
            header.contigs.add('chr1', length=1000)
            header.info.add('END', 1, 'Integer', 'End position')
            header.info.add('SVTYPE', 1, 'String', 'Variant type')
            header.formats.add('GT', 1, 'String', 'Genotype')
            header.add_sample('SYNTHETIC')
            record = header.new_record(contig='chr1', start=149, stop=250,
                                       alleles=('A', '<DEL>'), id='SYNTHETIC_DEL')
            record.info['SVTYPE'] = 'DEL'
            record.samples['SYNTHETIC']['GT'] = (0, 1)
            result = endpoint(record)
            self.assertEqual(result['pos'], 250)
            self.assertIsNone(result['local_side'])
            self.assertIsNone(result['remote_side'])
            with pysam.VariantFile(str(source), 'w', header=header) as writer:
                writer.write(record)
            summary = annotate_sv('SYNTHETIC', 'sniffles', source, root / 'genes.gff',
                                  root / 'genome.fa.fai', root / 'sv-output')
            self.assertEqual(summary['status'], 'completed')
            output = root / 'sv-output/SYNTHETIC.sniffles.structural.annotated.vcf.gz'
            with pysam.VariantFile(str(output)) as handle:
                annotated = next(iter(handle))
                self.assertEqual(annotated.samples['SYNTHETIC']['GT'], (0, 1))
                self.assertEqual(annotated.stop, 250)
                self.assertNotIn('STRANDS', annotated.info)

    def test_development_checksums_optional_but_supplied_hashes_checked(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = make(Path(tmp) / 'refs')
            path = root / 'references.json'
            manifest = json.loads(path.read_text())
            manifest['pinning'] = 'development'
            for entry in manifest['resources'].values():
                entry.pop('sha256')
                if entry.get('index'):
                    entry['index'].pop('sha256')
            path.write_text(json.dumps(manifest))
            refs, _ = validate(path)
            self.assertTrue(Path(refs['fasta']).is_file())
            manifest['resources']['fasta']['sha256'] = '0' * 64
            path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
                validate(path)
            manifest['resources']['fasta'].pop('sha256')
            manifest.pop('pinning')
            path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
                validate(path)

    def test_masked_reference_preparation_remains_strict_elsewhere(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / 'genome.fa'
            fasta.write_text('>chrY\nANAAAAAAAA\n')
            pysam.faidx(str(fasta))
            source = root / 'source.vcf'
            header = ('##fileformat=VCFv4.2\n##contig=<ID=Y,length=10>\n'
                      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
            source.write_text(header + 'Y\t2\tmasked\tC\tT\t.\t.\t.\n'
                              'Y\t4\tvalid\tA\tT\t.\t.\t.\n')
            with self.assertRaises(pysam.utils.SamtoolsError):
                prepare_vcf(str(source), str(fasta), str(root / 'strict.vcf.gz'))
            output = root / 'prepared.vcf.gz'
            audit = prepare_vcf(str(source), str(fasta), str(output), True)
            self.assertEqual(audit['masked_reference']['records'], 1)
            self.assertEqual(audit['masked_reference']['contigs'], {'Y': 1})
            with pysam.VariantFile(str(output)) as handle:
                self.assertEqual([(r.contig, r.pos, r.id) for r in handle], [('chrY', 4, 'valid')])
            with open(audit['masked_reference']['path'], 'rb') as stream, pysam.VariantFile(stream) as handle:
                self.assertEqual([r.id for r in handle], ['masked'])
            source.write_text(header + 'Y\t5\twrong\tC\tT\t.\t.\t.\n')
            with self.assertRaises(pysam.utils.SamtoolsError):
                prepare_vcf(str(source), str(fasta), str(root / 'wrong.vcf.gz'), True)

    def test_native_annotation_and_cosmic(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = make(Path(tmp) / 'refs')
            entries = json.loads((root / 'inputs.json').read_text())
            for row in entries:
                row['path'] = str(root / row['path'])
            (root / 'inputs.json').write_text(json.dumps(entries))
            prepare('SYNTHETIC', root / 'inputs.json', root / 'references.json', root / 'panel.tsv',
                    Path(tmp) / 'catalogue', Path(tmp) / 'payload')
            data = complete('SYNTHETIC', Path(tmp) / 'catalogue', Path(tmp) / 'payload',
                            root / 'fixture.vep.vcf.gz', root / 'references.json', root / 'panel.tsv',
                            Path(__file__).resolve().parents[2] / 'assets/fusion-list.txt', 100, Path(tmp) / 'output')
            self.assertEqual(len(data['small_variants']['records']), 2)
            self.assertTrue(all(r['exact_pass_agreement'] for r in data['small_variants']['records']))
            snv = next(r for r in data['small_variants']['records'] if r['coordinates'][1] == 150)
            hit = snv['evidence'][0]['cosmic']['sources'][0]
            self.assertEqual(len(hit['observations']), 2)
            self.assertEqual(hit['recurrence'][0]['values'], [7, 9])
            self.assertTrue(hit['recurrence'][0]['inconsistent'])
            with pysam.VariantFile(str(Path(tmp) / 'output/clair3/SYNTHETIC.clair3.germline.targets.annotated.vcf.gz')) as vcf:
                tagged = next(iter(vcf))
                self.assertEqual(tagged.info['ONCOSEQ_REVIEW_COSMIC'], ('COSV1',))
                self.assertTrue(tagged.info['ONCOSEQ_REVIEW_COSMIC_BLOOD'])
                self.assertEqual(list(tagged.filter), ['PASS'])
            joined = snv['evidence'][0]['cosmic']['targeted_evidence'][0]
            self.assertEqual(joined['blood_sample_count'], 1)
            self.assertEqual(joined['sample_count'], 1)
            self.assertEqual(len(joined['evidence']), 2)
            self.assertTrue(snv['selections'][0]['cosmic_primary'])
            indel = next(r for r in data['small_variants']['records'] if r['coordinates'][1] == 160)
            self.assertFalse(indel['selections'][0]['cosmic_primary'])
            self.assertEqual(indel['evidence'][0]['cosmic']['targeted_evidence'][0]['phenotypes'][0]['classification']['PRIMARY_SITE'], 'lung')
            self.assertEqual(data['coverage'][0]['primary'], 80)
            self.assertEqual(data['classification'][0]['score'], 0.9)
            html = render(Path(tmp) / 'output/review.json', Path(tmp) / 'report.html',
                          Path(__file__).resolve().parents[2] / 'assets/review/report.html')
            self.assertIn('SYNTHETIC', html.read_text())
            self.assertIn('"classy_summary":true', html.read_text())
            self.assertNotIn('https://cdn', html.read_text())

    def row(self, caller, pos, other=500, orientation=('right', 'left')):
        return dict(id=caller + ':' + str(pos), caller=caller, sample='SYNTHETIC', chrom='chr1', pos=pos,
                    partner=dict(chrom='chr2', pos=other, local_side=orientation[0], remote_side=orientation[1]),
                    evidence=dict(info=dict(SVTYPE='BND', SUPPORT=5), format={}), filter=['PASS'],
                    genes1=[dict(name='RUNX1')], genes2=[dict(name='RUNX1T1')],
                    gene_pair=['RUNX1', 'RUNX1T1'], coordinate_convention='VCF 1-based', source_records=[])

    def test_consensus_rejects_chains_and_ambiguity(self):
        known = {('RUNX1', 'RUNX1T1')}
        self.assertEqual(len(consolidate([self.row('a', 500), self.row('b', 600)], known, 100)), 1)
        chain = consolidate([self.row('a', 500), self.row('b', 590), self.row('c', 680)], known, 100)
        self.assertEqual(len(chain), 3)
        self.assertTrue(all(r['ambiguous'] for r in chain))
        duplicate = consolidate([self.row('a', 500), self.row('a', 510), self.row('b', 505)], known, 100)
        self.assertEqual(len(duplicate), 3)
        opposite = consolidate([self.row('a', 500), self.row('b', 500, orientation=('left', 'right'))], known)
        self.assertEqual(len(opposite), 2)
        self.assertEqual(junction(self.row('a', 500))[2][1], 499)

    def test_failed_filter_and_depth_semantics(self):
        ann = dict(population=dict(status='retained_rare'), clinvar=[], consequences=[])
        self.assertIsNone(small_candidate(ann, ['.'], {}))
        row = self.row('a', 500)
        row['filter'] = ['LowQual']
        self.assertFalse(consolidate([row], {('RUNX1', 'RUNX1T1')})[0]['selected'])
        text = 'chr1\t150\tx\tA\tT\t.\tPASS\t.\tGT:DP:AD\t0/1:100:70,30'
        evidence = sample_evidence(text, 'clairsto')['clairsto']
        self.assertAlmostEqual(evidence['af'], 0.3)
        self.assertEqual(evidence['supporting_reads'], 30)

    def test_all_bnd_encodings_and_nasvar_offsets(self):
        header = pysam.VariantHeader()
        header.contigs.add('chr1', length=1000)
        header.contigs.add('chr2', length=1000)
        for alt, expected in [('A[chr2:500[', ('right', 'left')), ('A]chr2:500]', ('right', 'right')),
                              ('[chr2:500[A', ('left', 'left')), (']chr2:500]A', ('left', 'right'))]:
            row = header.new_record(contig='chr1', start=499, alleles=('A', alt))
            result = endpoint(row)
            self.assertEqual((result['local_side'], result['remote_side']), expected)
        event = dict(gene1=dict(name='RUNX1', chr='chr1', pos=501),
                     gene2=dict(name='RUNX1T1', chr='chr2', pos=500), supporting_reads=10,
                     breakpoints=[dict(gene0_chr='chr1', gene0_pos=501, gene0_dir='left',
                                       gene1_chr='chr2', gene1_pos=500, gene1_dir='right', n_reads=4)])
        row = nasvar_events(dict(fusions=dict(fusions=[event])), 'SYNTHETIC')[0]
        self.assertEqual(junction(row), ('BND', ('chr1', 500, 'right'), ('chr2', 499, 'left')))
        self.assertEqual(row['evidence']['supporting_reads'], 4)

    def test_reference_hash_and_html_escaping(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = make(Path(tmp) / 'refs')
            (root / 'targets.bed').write_text('changed\n')
            with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
                validate(root / 'references.json')
            data = dict(schema_version=1, sample='SAMPLE<script>bad()</script>', assembly='GRCh38', figures=[])
            (root / 'review.json').write_text(json.dumps(data))
            target = render(root / 'review.json', root / 'preview.html',
                            Path(__file__).resolve().parents[2] / 'assets/review/report.html')
            self.assertNotIn('<script>bad()</script>', target.read_text())

    def test_original_multiallelic_genotype_is_retained(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = make(Path(tmp) / 'refs')
            header = pysam.VariantFile(root / 'clair3.vcf.gz').header.copy()
            path = root / 'multi.vcf.gz'
            with pysam.VariantFile(path, 'wz', header=header) as writer:
                row = writer.new_record(contig='chr1', start=149, alleles=('A', 'C', 'T'), id='multi')
                row.filter.add('PASS')
                row.samples['SYNTHETIC']['GT'] = (0, 2)
                row.samples['SYNTHETIC']['DP'] = 100
                row.samples['SYNTHETIC']['AD'] = (60, 10, 30)
                writer.write(row)
            (root / 'inputs.json').write_text(json.dumps([dict(kind='small:clairsto', path=str(path))]))
            prepare('SYNTHETIC', root / 'inputs.json', root / 'references.json', root / 'panel.tsv',
                    Path(tmp) / 'catalogue', Path(tmp) / 'payload')
            with sqlite3.connect(Path(tmp) / 'catalogue/catalog.sqlite') as db:
                originals = [json.loads(r[0]) for r in db.execute('SELECT native_json FROM native_records')]
                self.assertEqual(len(originals), 2)
                self.assertTrue(all('0/2' in records[0] for records in originals))
                for (record,) in db.execute('SELECT record FROM records'):
                    self.assertIsNone(sample_evidence(record, 'clairsto')['clairsto']['af'])

    def test_guarded_native_sv_sample_alias(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = make(Path(tmp) / 'refs')
            mapping = root / 'samples.txt'
            mapping.write_text('SYNTHETIC SAMPLE\n')
            path = root / 'renamed.vcf.gz'
            pysam.bcftools.reheader('-s', str(mapping), '-o', str(path), str(root / 'sniffles.sv.vcf.gz'), catch_stdout=False)
            result = annotate_sv('SYNTHETIC', 'sniffles', path, root / 'genes.gff', root / 'genome.fa.fai', root / 'sv')
            self.assertEqual(result['records'][0]['vcf_sample'], 'SAMPLE')
            self.assertEqual(support(result['records'][0]), 5)
            with self.assertRaisesRegex(ValueError, 'sample mismatch'):
                annotate_sv('OTHER', 'severus', path, root / 'genes.gff', root / 'genome.fa.fai', root / 'bad')


if __name__ == '__main__':
    unittest.main()
