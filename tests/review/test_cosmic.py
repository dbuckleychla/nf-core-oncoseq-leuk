"""Synthetic COSMIC joins, tissue selection and reference identity regressions."""
import copy
import csv
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'bin/review'))
import pysam
from cosmic import annotate
from indexed_reference import IndexedReference
from cosmic_evidence import CosmicEvidence, EVIDENCE_FIELDS, PHENOTYPE_FIELDS, prepare_evidence
from fixtures import make
from report_prioritization import small_candidate, consensus_candidate
from resources import validate
from variant_presentation import enrich_variants


class CosmicTests(unittest.TestCase):
    def test_positive_screens_distinct_samples_and_missing_joins(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = make(Path(tmp))
            mutation = root / 'cosmic.mutations.tsv.gz'
            with gzip.open(mutation, 'at') as stream:
                writer = csv.DictWriter(stream, fieldnames=[*EVIDENCE_FIELDS, 'POSITIVE_SCREEN'], delimiter='\t')
                base = dict(GENOMIC_MUTATION_ID='COSV1', COSMIC_SAMPLE_ID='SYNTHETIC_COSS2', COSMIC_PHENOTYPE_ID='SYNTHETIC_COSO1', PUBMED_PMID='SYNTHETIC_PMID2', MUTATION_SOMATIC_STATUS='Variant of unknown origin', GENE_SYMBOL='TP53', TRANSCRIPT_ACCESSION='TX1', POSITIVE_SCREEN='y')
                writer.writerow(base)
                writer.writerow(dict(base, COSMIC_SAMPLE_ID='SYNTHETIC_NEGATIVE', POSITIVE_SCREEN='n'))
                writer.writerow(dict(base, GENOMIC_MUTATION_ID='COSV4', COSMIC_PHENOTYPE_ID='SYNTHETIC_MISSING'))
                writer.writerow(dict(base, GENOMIC_MUTATION_ID='', POSITIVE_SCREEN='n'))
            path = root / 'joined.sqlite'
            audit = prepare_evidence(mutation, root / 'cosmic.classification.tsv.gz', root / 'cosmic.vcf.gz', 'synthetic', path)
            db = CosmicEvidence(path)
            try:
                hit = db.lookup('COSV1')
                self.assertEqual(hit['sample_count'], 2)
                self.assertEqual(hit['blood_sample_count'], 2)
                self.assertEqual(len(hit['evidence']), 3)
                self.assertTrue(hit['blood_match'])
                self.assertEqual(len(hit['phenotypes']), 1)
                self.assertEqual(len(hit['phenotypes'][0]['somatic_statuses']), 2)
                missing = db.lookup('COSV4')
                self.assertFalse(missing['blood_match'])
                self.assertEqual(missing['phenotypes'][0]['status'], 'missing_classification')
                self.assertEqual(db.lookup('COSV5')['status'], 'missing_targeted_evidence')
                self.assertEqual(audit['counts']['negative_rows'], 2)
                self.assertEqual(audit['counts']['unclassified_evidence_rows'], 1)
            finally:
                db.close()

    def test_mandatory_resources_and_bound_vcf(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = make(Path(tmp))
            path = root / 'references.json'
            original = json.loads(path.read_text())
            for name in ('cosmic_targeted', 'cosmic_evidence'):
                changed = copy.deepcopy(original)
                del changed['resources'][name]
                path.write_text(json.dumps(changed))
                with self.assertRaisesRegex(ValueError, 'Missing review resources'):
                    validate(path)
            changed = copy.deepcopy(original)
            changed['resources']['cosmic_targeted']['release'] = 'different-release'
            path.write_text(json.dumps(changed))
            with self.assertRaisesRegex(ValueError, 'releases do not match'):
                validate(path)
            changed = copy.deepcopy(original)
            changed['pinning'] = 'development'
            entry = changed['resources']['cosmic_targeted']
            entry['uri'] = 'clinvar.vcf.gz'
            entry.pop('sha256')
            entry['index'] = dict(uri='clinvar.vcf.gz.tbi')
            path.write_text(json.dumps(changed))
            with self.assertRaisesRegex(ValueError, 'different targeted VCF'):
                validate(path)

    def test_bad_headers_and_conflicting_classification_fail_atomically(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = make(Path(tmp))
            classification = root / 'cosmic.classification.tsv.gz'
            with gzip.open(classification, 'at') as stream:
                writer = csv.DictWriter(stream, fieldnames=PHENOTYPE_FIELDS, delimiter='\t')
                writer.writerow(dict(COSMIC_PHENOTYPE_ID='SYNTHETIC_COSO1', PRIMARY_SITE='lung'))
            path = root / 'bad.sqlite'
            with self.assertRaisesRegex(ValueError, 'Conflicting COSMIC phenotype'):
                prepare_evidence(root / 'cosmic.mutations.tsv.gz', classification, root / 'cosmic.vcf.gz', 'synthetic', path)
            self.assertFalse(path.exists())
            self.assertFalse(path.with_name('bad.sqlite.building').exists())
            classification.write_bytes(gzip.compress(b'WRONG_COLUMN\ninvalid\n'))
            with self.assertRaisesRegex(ValueError, 'Missing COSMIC TSV columns'):
                prepare_evidence(root / 'cosmic.mutations.tsv.gz', classification, root / 'cosmic.vcf.gz', 'synthetic', path)

    def test_exact_allele_join_never_matches_position_alone(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = make(Path(tmp))
            reference = IndexedReference(root / 'cosmic.vcf.gz')
            evidence = CosmicEvidence(root / 'cosmic.evidence.sqlite')
            try:
                header = pysam.VariantHeader()
                header.contigs.add('chr1', length=1000)
                for ref, alt, expected in [('A', 'T', 'matched'), ('A', 'C', 'no_match'), ('C', 'T', 'no_match')]:
                    record = header.new_record(contig='chr1', start=149, alleles=(ref, alt))
                    result = annotate({'cosmic_targeted': reference}, record, evidence)
                    self.assertEqual(result['status'], expected)
                    self.assertEqual(result['blood_match'], expected == 'matched')
            finally:
                reference.close()
                evidence.close()

    def test_primary_blood_expanded_other_and_filter_gates(self):
        annotation = dict(population=dict(status='retained_rare', combined_af=.001), clinvar=[], consequences=[],
                          cosmic=dict(status='matched', blood_match=True, targeted_evidence=[], sources=[dict(observations=[dict(info=dict(GENE='DEMO', TRANSCRIPT='TX1', HGVSC='TX1:c.1A>T'))])]))
        primary = small_candidate(annotation, ['PASS'], {})
        self.assertEqual(primary['review_group'], 'findings')
        self.assertEqual(primary['gene'], 'DEMO')
        self.assertIsNone(small_candidate(annotation, ['LowQual'], {}))
        self.assertIsNone(small_candidate(annotation, ['.'], {}))
        annotation['population']['status'] = 'excluded_common'
        self.assertIsNone(small_candidate(annotation, ['PASS'], {}))
        self.assertIsNone(consensus_candidate(annotation, ['PASS'], {}))
        annotation['population']['status'] = 'retained_rare'
        annotation['cosmic']['blood_match'] = False
        expanded = small_candidate(annotation, ['PASS'], {})
        self.assertEqual(expanded['review_group'], 'coding')
        for blood, expected in ((True, 'findings'), (False, 'coding')):
            annotation['cosmic']['blood_match'] = blood
            row = dict(selections=[small_candidate(annotation, ['PASS'], {})], evidence=[dict(cosmic=annotation['cosmic'], clinvar=[dict(CLNSIG='Benign')])])
            enrich_variants(dict(small_variants=dict(records=[row])))
            self.assertEqual(row['presentation']['group'], expected)
            self.assertTrue(any('benign' in r.lower() for r in row['presentation']['reasons']))


if __name__ == '__main__':
    unittest.main()
