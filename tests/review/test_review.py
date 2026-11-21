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
from adapters import nasvar_events
from structural_consensus import support
from annotate_structural_variants import annotate as annotate_sv
import pysam
import pysam.bcftools
import sqlite3


class ReviewTests(unittest.TestCase):
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
            self.assertEqual(data['coverage'][0]['primary'], 80)
            self.assertEqual(data['classification'][0]['score'], 0.9)
            html = render(Path(tmp) / 'output/review.json', Path(tmp) / 'report.html',
                          Path(__file__).resolve().parents[2] / 'assets/review/report.html')
            self.assertIn('SYNTHETIC', html.read_text())
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
