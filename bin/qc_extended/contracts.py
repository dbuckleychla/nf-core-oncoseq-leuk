"""Versioned metric contracts and safe serialization."""
import csv
import hashlib
import json
import math
import re
from pathlib import Path

VERSION = '1.0.0'
FILES = ('summary.tsv', 'read_lengths.tsv', 'coverage.tsv', 'targets.tsv',
         'methylation.tsv', 'classifier_cpgs.tsv', 'adequacy.tsv', 'preprocessing.tsv')
DIMENSIONS = ('scope', 'mapq', 'target_id', 'panel_id', 'stage', 'modification', 'threshold')
FORMULAS = {
    'alignment.rate': 'Mapped primary records / all primary records, including unmapped and QC-failed',
    'coverage.mean_depth': 'Sum of eligible primary aligned-base depth / reference coordinate bases, including zeros and Ns',
    'coverage.breadth': 'Reference coordinate bases with eligible depth at least threshold / all coordinate bases in scope',
    'coverage.zero_fraction': 'Reference coordinate bases with zero eligible depth / all coordinate bases in scope',
    'coverage.mapq_retention': 'Mean depth at numerical MAPQ >=20 / mean depth at numerical MAPQ >=0',
    'coverage.observed_enrichment': 'Mean enrichment coordinate depth / mean off-enrichment coordinate depth',
    'read_length.mean': 'Known query bases including soft clips / primary records with known query lengths',
    'read_length.median': 'Exact record-weighted median of known primary query lengths',
    'read_length.n50': 'Descending known primary query length at cumulative half of total known query bases',
    'allocation.target_aligned_base_fraction': 'Exact-target aligned bases / all primary-genome aligned bases',
    'allocation.padding_aligned_base_fraction': 'Aligned bases in enrichment minus exact targets / enrichment aligned bases',
    'allocation.target_yield_per_primary_query_base': 'Primary-genome exact-target aligned bases / known query bases of all primary records; unavailable if any lengths are unknown',
    'allocation.targets_outside_enrichment_fraction': 'Exact-target union coordinate bases outside enrichment / exact-target union bases',
    'adequacy.target_breadth': 'Exact-target bases meeting 95% binomial alternate-support sampling probability / all target bases; assumed VAF and support in threshold',
    'methylation.callable_fraction': 'Dyads with valid modification calls / reference CpG dyads in scope',
    'methylation.mean_beta': 'Sum of site modification fractions / valid CpG sites; each dyad counted once',
    'methylation.weighted_beta': 'Modified calls / valid calls, with one shared denominator across m/h rows',
    'methylation.filter_fail_fraction': 'Filter-fail calls / (valid + filter-fail calls)',
    'methylation.no_call_fraction': 'No-call observations / (valid + filter-fail + no-call observations)',
    'methylation.reference_fraction_at_depth': 'CpG dyads reaching threshold valid depth / all reference CpG dyads in scope',
    'classifier_features.reference_fraction_at_depth': 'Unique panel CpG features reaching threshold valid depth / unique reference panel CpGs',
}


def identifier(value):
    if not value or value in ('.', '..') or re.search(r'[/\\\x00-\x1f]', value):
        raise ValueError('Unsafe QC identifier')
    return value


def checksum(path):
    h = hashlib.sha256()
    with open(path, 'rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, data):
    destination = Path(path)
    temporary = destination.with_suffix(destination.suffix + '.tmp')
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    temporary.replace(destination)


def state(status='complete', reason=None, exhaustive=True):
    return dict(status=status, reason=reason, exhaustive=exhaustive)


class Metrics:
    def __init__(self):
        self.metrics, self.histograms, self.definitions = [], [], {}

    def add(self, name, value, unit='records', numerator=None, denominator=None,
            reason=None, description=None, **dimensions):
        if isinstance(value, float) and not math.isfinite(value):
            raise ValueError('Nonfinite QC metric: ' + name)
        unknown = set(dimensions) - set(DIMENSIONS)
        if unknown:
            raise ValueError('Unknown metric dimensions: ' + str(sorted(unknown)))
        self.metrics.append(dict(metric_id=name, **{d: dimensions.get(d) for d in DIMENSIONS},
                                 unit=unit, value=value, numerator=numerator, denominator=denominator,
                                 status='complete' if value is not None else 'unavailable',
                                 reason=None if value is not None else reason or 'zero_denominator'))
        self.definitions[name] = dict(description=description or name.replace('_', ' ').replace('.', ': '),
                                      unit=unit, policy='See provenance.filter_policies and metric dimensions',
                                      formula=FORMULAS.get(name, 'Observed aggregate or exact descriptive statistic in the declared unit and scope'),
                                      null_policy='Unavailable observations or an undefined denominator; never substituted with zero')

    def ratio(self, name, numerator, denominator, **dimensions):
        self.add(name, numerator / denominator if denominator else None, 'fraction',
                 numerator, denominator, **dimensions)

    def histogram(self, name, values, unit, denominator, width=None, **dimensions):
        bins = []
        for lower, count in sorted(values.items(), key=lambda item: float(item[0])):
            lower = float(lower) if width else int(lower)
            upper = lower + width if width else lower
            bins.append(dict(lower=lower, upper=upper, upper_inclusive=width is None or upper >= 1 and unit == 'beta', count=count))
        self.histograms.append(dict(metric_id=name, **{d: dimensions.get(d) for d in DIMENSIONS},
                                    unit=unit, binning='fixed_width' if width else 'discrete',
                                    denominator=denominator, status='complete', reason=None, bins=bins))


def write_tables(directory, data):
    columns = ['sample_id', 'input_scope', *DIMENSIONS, 'contig', 'start', 'end', 'label',
               'assumed_vaf', 'min_support', 'probability_target', 'required_depth',
               'scope_type', 'scope_id', 'exhaustive', 'records_examined', 'model_name', 'model_version', 'assembly',
               'metric_id', 'unit', 'value', 'numerator', 'denominator', 'status', 'reason']
    targets = {r['target_id']: r for r in data['provenance'].get('target_intervals', [])}
    panels = {r['panel_id']: r for r in data['provenance'].get('classifier_panels', [])}
    for filename in FILES:
        with open(Path(directory) / filename, 'w') as out:
            writer = csv.DictWriter(out, fieldnames=columns, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            for metric in data['metrics']:
                name = metric['metric_id']
                selected = {
                    'summary.tsv': not metric['target_id'] and not metric['panel_id'],
                    'read_lengths.tsv': name.startswith(('read_length.', 'read_quality.')),
                    'coverage.tsv': name.startswith(('coverage.', 'allocation.')) and not metric['target_id'],
                    'targets.tsv': bool(metric['target_id']) and not name.startswith('adequacy.'),
                    'methylation.tsv': name.startswith('methylation.'),
                    'classifier_cpgs.tsv': name.startswith('classifier_features.'),
                    'adequacy.tsv': name.startswith('adequacy.'),
                    'preprocessing.tsv': name.startswith('preprocessing.'),
                }[filename]
                if not selected:
                    continue
                row = dict(sample_id=data['sample_id'], input_scope=data['input_scope'], **metric)
                if metric['target_id']:
                    row.update({k: targets[metric['target_id']][k] for k in ('contig', 'start', 'end', 'label')})
                if metric['panel_id']:
                    row.update({k: panels[metric['panel_id']][k] for k in ('model_name', 'model_version', 'assembly')})
                if name.startswith('adequacy.'):
                    v, k = metric['threshold'].split(':')
                    row.update(assumed_vaf=v, min_support=k, probability_target=.95,
                               required_depth=required_depth(float(v), int(k)))
                if name.startswith('preprocessing.'):
                    scope_type, scope_id = metric['scope'].split(':', 1)
                    row.update(scope_type=scope_type, scope_id=scope_id, exhaustive=True)
                writer.writerow({k: 'NA' if row.get(k) is None else row[k] for k in columns})


def required_depth(vaf, support):
    if not 0 < vaf < 1 or support < 1:
        raise ValueError('Invalid binomial sampling assumptions')
    n = support
    while True:
        probability = (1 - vaf) ** n
        cumulative = probability
        for j in range(1, support):
            probability *= (n - j + 1) / j * vaf / (1 - vaf)
            cumulative += probability
        if cumulative <= .05:
            return n
        n += 1
