import copy
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'bin/review'))
from variant_presentation import enrich_variants, cosmic_summary


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


class NonCosmicPriorityTests(unittest.TestCase):
    def test_previous_priority_criteria_get_a_separate_table(self):
        pathogenic=record([selection('ClinVar pathogenic / likely pathogenic')], [dict(CLNSIG='Pathogenic',CLNREVSTAT='reviewed_by_expert_panel')])
        lof=record([selection('Predicted loss of function','HIGH')])
        broad=record([selection('VEP HIGH','HIGH','consensus_review')])
        benign=record([selection('Predicted loss of function','HIGH')],[dict(CLNSIG='Benign')])
        common=record([selection('Predicted loss of function','HIGH')],population=dict(status='excluded_common',combined_af=.1))
        rows=[pathogenic,lof,broad,benign,common]
        for row in rows:row['evidence'][0]['cosmic']=dict(status='no_match',sources=[],targeted_evidence=[],blood_match=False)
        enrich_variants(dict(small_variants=dict(records=rows)))
        self.assertEqual([r['presentation']['group'] for r in rows],['noncosmic_findings','noncosmic_findings','coding','coding','coding'])
        self.assertTrue(all(not r['presentation']['cosmic']['matched'] for r in rows))


class CosmicSummaryTests(unittest.TestCase):
    def test_distinct_samples_and_publications_across_ids_and_callers(self):
        phenotype=dict(phenotype_id='SYNTHETIC_PHENOTYPE', classification=dict(PRIMARY_SITE='haematopoietic_and_lymphoid_tissue'))
        evidence=dict(COSMIC_SAMPLE_ID='SYNTHETIC_SAMPLE', COSMIC_PHENOTYPE_ID='SYNTHETIC_PHENOTYPE', PUBMED_PMID='SYNTHETIC_PUBLICATION', MUTATION_SOMATIC_STATUS='Confirmed somatic variant')
        target=dict(identifier='COSV1', status='matched', phenotypes=[phenotype], evidence=[evidence, dict(evidence, TRANSCRIPT_ACCESSION='SYNTHETIC_TX')])
        cosmic=dict(blood_match=True, sources=[dict(identifiers=['COSV1','COSV2'], observations=[dict(info=dict(GENE='DEMO',TIER='1')),dict(info=dict(GENE='OTHER',TIER='2'))])], targeted_evidence=[target,dict(target,identifier='COSV2')])
        row=dict(evidence=[dict(cosmic=cosmic),dict(cosmic=copy.deepcopy(cosmic))])
        summary=cosmic_summary(row)
        self.assertEqual((summary['sample_count'],summary['blood_sample_count'],summary['publication_count']),(1,1,1))
        target['evidence'].append(dict(evidence,COSMIC_SAMPLE_ID='SYNTHETIC_SAMPLE_2'))
        row['evidence'][1]['cosmic']=copy.deepcopy(cosmic)
        multiple_samples=cosmic_summary(row)
        self.assertEqual((multiple_samples['sample_count'],multiple_samples['publication_count']),(2,1))
        self.assertNotIn('tiers',summary)
        self.assertEqual(summary['somatic_categories'],['Confirmed somatic reports'])
        self.assertTrue(summary['tissue_known'])
        row['selections']=[selection('COSMIC blood/lymphoid match')]
        enrich_variants(dict(small_variants=dict(records=[row])))
        self.assertFalse(any('do not establish somatic origin' in r for r in row['presentation']['reasons']))

    def test_missing_join_has_unavailable_counts(self):
        row=dict(evidence=[dict(cosmic=dict(sources=[dict(identifiers=['COSV1'])], targeted_evidence=[dict(identifier='COSV1',status='missing_targeted_evidence',phenotypes=[],evidence=[])]))])
        summary=cosmic_summary(row)
        self.assertTrue(summary['matched'])
        self.assertIsNone(summary['sample_count'])
        self.assertFalse(summary['tissue_known'])
        row['evidence'][0]['cosmic']['sources'][0]['identifiers'].append('COSV2')
        row['evidence'][0]['cosmic']['targeted_evidence'][0]['status']='matched'
        self.assertIsNone(cosmic_summary(row)['sample_count'])
