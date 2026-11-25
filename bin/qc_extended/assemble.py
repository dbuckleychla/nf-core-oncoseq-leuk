"""Assemble public metric contracts from private calculation components."""
from collections import Counter

from .alignment import THRESHOLDS, depth_summary
from .contracts import FILES, Metrics, identifier, required_depth, state


def assemble(sample, alignment, methylation, reference, preprocessing=(), provenance=None):
    identifier(sample)
    m = Metrics()
    for key, value in alignment['counts'].items():
        m.add('alignment.' + key, value, scope='all')
    counts = alignment['counts']
    m.ratio('alignment.rate', counts.get('mapped_primary', 0), counts.get('primary_records', 0), scope='all')
    for key, value in alignment['modifications'].items():
        m.add('alignment.modification_' + key, value, scope='all')
    for group, data in alignment['read_lengths'].items():
        for key in ('records', 'known_length_records', 'missing_length_records', 'bases', 'mean', 'median', 'n50', 'minimum', 'maximum'):
            m.add('read_length.' + key, data[key], 'bp' if key in ('bases', 'mean', 'median', 'n50', 'minimum', 'maximum') else 'records', scope=group)
        m.histogram('read_length.distribution', data['histogram'], 'bp', data['known_length_records'], width=100, scope=group)
        for bin_data in m.histograms[-1]['bins']:
            key = int(bin_data['lower'])
            bin_data['base_sum'] = data['histogram_bases'].get(key, data['histogram_bases'].get(str(key)))
        for label, fractions in data.get('fractions', {}).items():
            m.ratio('read_length.' + label + '_read_fraction', fractions['reads'], data['known_length_records'], scope=group)
            m.ratio('read_length.' + label + '_base_fraction', fractions['bases'], data['bases'], scope=group)
    m.histogram('alignment.mapq', alignment['mapq'], 'mapq', counts.get('mapped_primary', 0), scope='mapped')
    m.add('alignment.mapq_unavailable_records', alignment['mapq'].get(255, alignment['mapq'].get('255', 0)), scope='mapped')
    mapped = alignment['read_lengths']['mapped']
    m.add('alignment.soft_clipped_bases', alignment['soft_clipped_bases'], 'bp', scope='mapped')
    m.ratio('alignment.soft_clipped_fraction', alignment['soft_clipped_bases'], mapped['bases'], scope='mapped')
    for key, value in alignment['quality'].items():
        m.add('read_quality.' + key, value, 'phred' if key in ('qs_mean', 'error_phred') else 'fraction' if key == 'error_probability' else 'bp' if key == 'error_bases' else 'records', scope='all')

    def coverage(data, scope, mapq, target_id=None, baseline=None):
        dims = dict(scope=scope, mapq=mapq, target_id=target_id)
        for key in ('bases', 'aligned_bases', 'mean', 'median', 'cv', 'fraction_ge_0_2_mean'):
            name = {'mean': 'mean_depth', 'median': 'median_depth'}.get(key, key)
            m.add('coverage.' + name, data[key], 'bp' if key in ('bases', 'aligned_bases') else 'depth' if key in ('mean', 'median') else 'fraction', **dims)
        for threshold in THRESHOLDS:
            m.add('coverage.breadth', data['breadth'][str(threshold)], 'fraction',
                  sum(v for k, v in data['histogram'].items() if int(k) >= threshold), data['bases'], threshold=threshold, **dims)
        if data['bases']:
            m.ratio('coverage.zero_fraction', data['histogram'].get(0, data['histogram'].get('0', 0)), data['bases'], **dims)
        if baseline is not None:
            m.add('coverage.mapq_retention', data['mean'] / baseline['mean'] if baseline['mean'] else None, 'fraction', **dims)
        if not target_id:
            m.histogram('coverage.depth_distribution', data['histogram'], 'depth', data['bases'], **dims)
        if scope == 'targets':
            for vaf in (.01, .02, .05, .1):
                for support in (3, 5):
                    depth = required_depth(vaf, support)
                    n = sum(v for k, v in data['histogram'].items() if int(k) >= depth)
                    m.add('adequacy.target_breadth', n / data['bases'] if data['bases'] else None, 'fraction',
                          n, data['bases'], threshold=f'{vaf}:{support}', **dims)

    for mapq, scopes in alignment['coverage'].items():
        for scope, data in scopes.items():
            coverage(data, scope, int(mapq), baseline=alignment['coverage']['0'][scope] if mapq == '20' else None)
        e, bg, t, genome = (scopes[k] for k in ('enrichment', 'off_enrichment', 'targets', 'genome_primary'))
        m.add('coverage.observed_enrichment', e['mean'] / bg['mean'] if bg['mean'] and e['mean'] is not None else None, 'ratio', scope='enrichment', mapq=int(mapq))
        m.ratio('allocation.target_aligned_base_fraction', t['aligned_bases'], genome['aligned_bases'], scope='genome_primary', mapq=int(mapq))
        m.ratio('allocation.padding_aligned_base_fraction', scopes['padding']['aligned_bases'], e['aligned_bases'], scope='enrichment', mapq=int(mapq))
        total = alignment['read_lengths']['all']
        m.add('allocation.target_yield_per_primary_query_base', t['aligned_bases'] / total['bases'] if total['bases'] and not total['missing_length_records'] else None,
              'fraction', t['aligned_bases'], total['bases'] if not total['missing_length_records'] else None,
              reason='incomplete_query_length_denominator' if total['missing_length_records'] else None, scope='targets', mapq=int(mapq))
    for contig in alignment['contigs']:
        for mapq, rows in contig['coverage'].items():
            coverage(rows['scopes']['genome_primary'], 'contig:' + contig['contig'], int(mapq))
            for key, data in rows['targets'].items():
                coverage(data, 'targets', int(mapq), key, contig['coverage']['0']['targets'][key] if mapq == '20' else None)
    m.ratio('allocation.targets_outside_enrichment_fraction', alignment['targets_outside_enrichment_bases'], alignment['target_union_bases'], scope='targets')

    for scope, mods in methylation.get('scopes', {}).items():
        for mod, data in mods.items():
            dims = dict(scope=scope, modification=mod)
            available = data['available']
            for key in ('reference_cpgs', 'sites', 'valid', 'modified', 'filter_fail', 'no_call', 'deletion', 'differing_base'):
                m.add('methylation.' + key, data[key] if available else None, 'sites' if key in ('reference_cpgs', 'sites') else 'calls', reason='modification_not_measured', **dims)
            for name, numerator, denominator in [('callable_fraction', data['sites'], data['reference_cpgs']),
                                                 ('mean_beta', data['beta_sum'], data['sites']),
                                                 ('weighted_beta', data['modified'], data['valid']),
                                                 ('filter_fail_fraction', data['filter_fail'], data['valid'] + data['filter_fail']),
                                                 ('no_call_fraction', data['no_call'], data['valid'] + data['filter_fail'] + data['no_call'])]:
                m.add('methylation.' + name, numerator / denominator if denominator and available else None, 'fraction',
                      numerator if available else None, denominator if available else None,
                      reason='modification_not_measured' if not available else None, **dims)
            if available:
                for threshold in (1, 5, 10, 20):
                    n = sum(v for k, v in data['depth'].items() if int(k) >= threshold)
                    m.add('methylation.sites_at_depth', n, 'sites', threshold=threshold, **dims)
                    m.ratio('methylation.reference_fraction_at_depth', n, data['reference_cpgs'], threshold=threshold, **dims)
                m.histogram('methylation.valid_depth', data['depth'], 'calls', data['sites'], **dims)
                for name, key in [('beta', 'beta'), ('beta_depth10', 'beta10')]:
                    m.histogram('methylation.' + name, {int(k) / 100: v for k, v in data[key].items()}, 'beta', sum(data[key].values()), width=.01, **dims)
    for panel in methylation.get('features', []):
        dims = dict(scope='classifier_panel', panel_id=panel['panel_id'])
        m.add('classifier_features.reference_sites', panel['size'], 'sites', **dims)
        m.add('classifier_features.duplicate_input_rows', panel['input_rows'] - panel['size'], 'sites', **dims)
        depths = {int(k): v for k, v in panel['depth'].items()} if panel['depth'] is not None else None
        m.add('classifier_features.median_valid_depth', depth_summary(depths, panel['size'])['median'] if depths is not None else None,
              'calls', reason=methylation['status']['reason'] if depths is None else None, **dims)
        for threshold in (1, 5, 10, 20):
            n = sum(v for k, v in panel['depth'].items() if int(k) >= threshold) if panel['depth'] is not None else None
            m.add('classifier_features.reference_fraction_at_depth', n / panel['size'] if n is not None and panel['size'] else None,
                  'fraction', n, panel['size'], reason=methylation['status']['reason'] if n is None else None, threshold=threshold, **dims)
    for row in preprocessing:
        for key in ('reads', 'bases'):
            m.add('preprocessing.' + key, row[key], 'records' if key == 'reads' else 'bp', scope=row['scope_type'] + ':' + row['scope_id'], stage=row['stage'])
    # The mapping stage is always measurable even when import history is absent.
    m.add('preprocessing.reads', counts.get('primary_records', 0), stage='mapping_output', scope='sample:' + sample)
    m.add('preprocessing.bases', alignment['read_lengths']['all']['bases'], 'bp', stage='mapping_output', scope='sample:' + sample)
    status = methylation.get('status', state('unavailable', 'no_usable_modifications'))
    partial = status['status'] in ('unavailable', 'partial') or alignment['read_lengths']['all']['missing_length_records'] > 0
    warnings = []
    for condition, code in [(not counts.get('mapped_primary'), 'no_mapped_primary_reads'),
                            (alignment['targets_outside_enrichment_bases'] > 0, 'targets_outside_enrichment'),
                            (len(alignment['models']) > 1, 'mixed_basecalling_models')]:
        if condition:
            warnings.append(dict(code=code, component='alignment', count=1, message=code.replace('_', ' ')))
    if status['reason']:
        warnings.append(dict(code=status['reason'], component='methylation', count=1, message=status['reason'].replace('_', ' ')))
    if methylation.get('scopes') and not methylation['scopes']['genome_primary']['combined']['sites']:
        warnings.append(dict(code='no_valid_cpg_calls', component='methylation', count=1, message='No valid primary-genome CpG calls observed'))
    maximum_depth = max((int(k) for k, v in alignment['coverage']['0']['genome_primary']['histogram'].items() if v), default=0)
    if maximum_depth >= 8000:
        warnings.append(dict(code='potential_methylation_depth_cap', component='methylation', count=1,
                             message='Source depth may exceed the 8000-record modkit cap; eligible-record policies differ'))
    provenance = dict(provenance or {}, reference_sha256=reference['reference_sha256'],
                      target_intervals=alignment['target_intervals'], basecalling_models=alignment['models'],
                      classifier_panels=[{k: v for k, v in p.items() if k != 'depth'} for p in methylation.get('features', [])],
                      modkit_version=methylation.get('modkit_version', 'not_run'),
                      filter_policies=dict(coverage='Primary; includes duplicate/QC-failed; MAPQ 0/20; BQ0; excludes D/N',
                                           lengths='Primary query lengths including soft clips, missing lengths excluded explicitly',
                                           methylation=methylation.get('filtering', {})))
    return dict(schema_name='oncoseq_extended_qc', schema_version='1.0.0', sample_id=sample,
                input_scope='full_mapping_stage_sample', status='partial' if partial else 'complete',
                components=dict(alignment=state(), coverage=state(), methylation=status,
                                classifier_features=(state() if status['status'] in ('complete', 'partial') else status) if methylation.get('features') else state('unavailable', 'no_panel_supplied'),
                                preprocessing=state('partial', 'prior_stages_unavailable') if not preprocessing else state()),
                provenance=provenance, definitions=m.definitions, metrics=m.metrics, histograms=m.histograms,
                warnings=warnings, tables=[dict(basename=f, schema_version='1.0.0') for f in FILES])
