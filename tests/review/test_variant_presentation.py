import copy
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'bin/review'))
from variant_presentation import enrich_variants


def selection(reason, impact='MODERATE', group='findings'):
    return dict(gene='DEMO', reason=reason, review_group=group, panel_gene=True, priority=0,
                consequence='missense_variant', gene_context=dict(transcript=dict(SYMBOL='DEMO', Feature='TX1', IMPACT=impact,
                    Consequence='missense_variant', BIOTYPE='protein_coding')), inclusion_reasons=[reason])


def record(selections, assertions=None, population=None):
    return dict(coordinates=['chr1', 100, 'A', 'T'], selections=selections,
                evidence=[dict(caller='clairsto', filter='PASS', clinvar=assertions or [], population=population or {})], exact_pass_agreement=True)


class VariantPresentationTests(unittest.TestCase):
    def test_reviewed_pathogenic_is_not_conflicting_due_to_no_conflicts_text(self):
        ordinary=selection('ClinVar pathogenic / likely pathogenic')
        broad=selection('ClinVar P/LP', group='consensus_review')
        assertions=[dict(CLNSIG=['Pathogenic'], CLNREVSTAT=['criteria_provided','_multiple_submitters','_no_conflicts'])]
        row=record([ordinary,broad], assertions, dict(combined_af=.00001,status='retained_rare'))
        original=copy.deepcopy(row)
        enrich_variants(dict(small_variants=dict(records=[row])))
        p=row['presentation']
        self.assertEqual((p['group'], p['rank'], p['review_stars']), ('findings',0,2))
        self.assertEqual(p['selected_impact'], 'MODERATE')
        self.assertTrue(any('Maximum reported population frequency' in r for r in p['reasons']))
        self.assertEqual(row['selections'], original['selections'])
        self.assertEqual(row['evidence'], original['evidence'])
        self.assertEqual(len(p['transcripts']),1)

    def test_benign_and_common_broad_calls_remain_lower_priority(self):
        benign=record([selection('VEP HIGH', 'HIGH', 'consensus_review')], [dict(CLNSIG=['Benign'])])
        common=record([selection('VEP HIGH', 'HIGH', 'consensus_review')], population=dict(combined_af=.1,status='excluded_common'))
        enrich_variants(dict(small_variants=dict(records=[benign,common])))
        self.assertEqual(benign['presentation']['group'],'supporting')
        self.assertEqual(common['presentation']['group'],'supporting')
        self.assertGreater(common['presentation']['rank'],benign['presentation']['rank'])
        self.assertTrue(any('Common-population' in r for r in common['presentation']['reasons']))

    def test_generic_high_is_not_promoted_to_curated_loss_of_function(self):
        broad=record([selection('VEP HIGH', 'HIGH', 'consensus_review')])
        selected=selection('Predicted loss of function','HIGH')
        selected['consequence']='stop_gained'
        selected['gene_context']['transcript']['Consequence']='stop_gained'
        lof=record([selected])
        enrich_variants(dict(small_variants=dict(records=[broad,lof])))
        self.assertEqual(broad['presentation']['group'],'coding')
        self.assertEqual(lof['presentation']['group'],'findings')
        self.assertLess(lof['presentation']['rank'],broad['presentation']['rank'])

    def test_selected_impact_is_not_replaced_by_high_on_another_transcript(self):
        selected=selection('VEP HIGH', 'MODERATE', 'consensus_review')
        selected['transcripts']=[dict(SYMBOL='DEMO', Feature='TX2', IMPACT='HIGH', Consequence='splice_donor_variant')]
        row=record([selected])
        enrich_variants(dict(small_variants=dict(records=[row])))
        self.assertEqual(row['presentation']['selected_impact'],'MODERATE')
        self.assertTrue(any(t['impact']=='HIGH' and 'Broader review' in t['roles'] for t in row['presentation']['transcripts']))
        self.assertEqual(row['presentation']['group'],'coding')
