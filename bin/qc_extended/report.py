"""Portable report content from aggregate QC only; no external assets."""
from html import escape
import math


def validate(data):
    if data.get('schema_name') != 'oncoseq_extended_qc' or data.get('schema_version') != '1.0.0':
        raise ValueError('Unsupported extended QC schema')
    for row in data['metrics']:
        if row['value'] is not None and (not isinstance(row['value'], (int, float)) or not math.isfinite(row['value'])):
            raise ValueError('Invalid metric value')


def table(headers, rows):
    return '<div class="qc-ext-table"><table><thead><tr>' + ''.join('<th>' + escape(str(h)) + '</th>' for h in headers) + '</tr></thead><tbody>' + ''.join('<tr>' + ''.join('<td>' + escape(str(v)) + '</td>' for v in row) + '</tr>' for row in rows) + '</tbody></table></div>'


def value(row):
    n = row['value']
    if n is None:
        return 'Unavailable (' + (row.get('reason') or 'unknown') + ')'
    if row['unit'] == 'fraction':
        return f'{100*n:.2f}%'
    return f'{n:,.3f}'.rstrip('0').rstrip('.') if isinstance(n, float) else f'{n:,}'


def chart(title, series, log=False):
    """Histogram lines with per-series denominators and an explicit bin axis."""
    available = [s for s in series if s.get('bins') and s.get('denominator')]
    if not available:
        return '<h4>' + escape(title) + '</h4><p>No observed distribution available.</p>'
    xs = [max(1, (b['lower'] + b['upper']) / 2) if log else (b['lower'] + b['upper']) / 2 for s in available for b in s['bins']]
    low, high = min(xs), max(xs)
    if log:
        low, high = math.log10(low), math.log10(high)
    high = max(high, low + .01)
    peak = max(b['count'] / s['denominator'] for s in available for b in s['bins']) or 1
    colors = ('#1e88e5', '#e67700', '#8e44ad')
    svg = '<svg viewBox="0 0 640 240" role="img" aria-label="' + escape(title) + '"><path d="M50 15V195H620" fill="none" stroke="currentColor"/>'
    for i, s in enumerate(available):
        points = []
        marks = ''
        for b in s['bins']:
            x = (b['lower'] + b['upper']) / 2
            x = math.log10(max(1, x)) if log else x
            px, py = 50+(x-low)/(high-low)*570, 195-b['count']/s['denominator']/peak*170
            points.append(f'{px:.2f},{py:.2f}')
            marks += f'<circle cx="{px:.2f}" cy="{py:.2f}" r="3" fill="{colors[i%3]}"><title>{b["lower"]:g}–{b["upper"]:g}: {b["count"]} observations</title></circle>'
        svg += f'<polyline points="{" ".join(points)}" fill="none" stroke="{colors[i%3]}" stroke-width="2"/>'
        svg += marks
        svg += f'<text x="{55+i*185}" y="235" fill="{colors[i%3]}" font-size="12">{escape(s.get("scope") or s.get("modification") or "distribution")}</text>'
    svg += f'<text x="5" y="25" font-size="11">{peak*100:.1f}%</text><text x="50" y="215" font-size="11">{10**low if log else low:g}</text><text x="575" y="215" font-size="11">{10**high if log else high:g}</text></svg>'
    return '<h4>' + escape(title) + '</h4>' + svg + '<p class="qc-ext-note">Fraction of observations in each group; ' + ('100-bp bin midpoints on a logarithmic length axis.' if log else 'original histogram bins.') + '</p>'


def content(data):
    validate(data)
    rows = data['metrics']
    chosen = {'alignment.rate', 'alignment.primary_records', 'coverage.mean_depth', 'coverage.observed_enrichment',
              'allocation.target_aligned_base_fraction', 'allocation.padding_aligned_base_fraction',
              'methylation.callable_fraction', 'methylation.mean_beta', 'methylation.weighted_beta'}
    summary = [[r['metric_id'].replace('.', ': ').replace('_', ' '), r['scope'], r['mapq'] if r['mapq'] is not None else '', r['modification'] or '', value(r)]
               for r in rows if r['metric_id'] in chosen and not r['target_id'] and r['scope'] in ('all', 'targets', 'enrichment', 'off_enrichment', 'genome_primary')]
    out = '<section class="qc-ext"><style>.qc-ext svg{width:100%;max-width:760px;color:inherit}.qc-ext-table{overflow:auto}.qc-ext table{border-collapse:collapse;width:100%;margin:1em 0}.qc-ext th,.qc-ext td{padding:.45em;border-bottom:1px solid #8885;text-align:left}.qc-ext-note{font-size:.85em;opacity:.8}</style>'
    out += '<h3>Extended sequencing QC</h3><p>Full mapping-stage sample: ' + escape(data['sample_id']) + '. These metrics describe retained input records before downstream filtering, phasing, or classifier subsampling. They are not window-specific or clinical pass/fail results.</p>'
    out += table(['Component', 'Status', 'Reason'], [[k.replace('_', ' '), s['status'], s['reason'] or ''] for k, s in data['components'].items()])
    out += table(['Metric', 'Region', 'MAPQ', 'Modification', 'Value'], summary)
    for group, title in [('target', 'Read lengths on/off exact targets'), ('enrichment', 'Read lengths on/off enrichment')]:
        out += chart(title, [h for h in data['histograms'] if h['metric_id'] == 'read_length.distribution' and h['scope'] in ('on_' + group, 'off_' + group)], True)
    out += chart('CpG 5mC beta: primary genome', [h for h in data['histograms'] if h['metric_id'] == 'methylation.beta' and h['scope'] == 'genome_primary' and h['modification'] == '5mC'])
    out += chart('Valid CpG depth: exact targets and off enrichment', [h for h in data['histograms'] if h['metric_id'] == 'methylation.valid_depth' and h['scope'] in ('targets', 'off_enrichment') and h['modification'] == 'combined'])
    lengths = [[r['scope'], r['metric_id'].split('.')[-1], value(r)] for r in rows if r['metric_id'] in ('read_length.records', 'read_length.mean', 'read_length.median', 'read_length.n50')]
    out += '<h4>Regional read statistics</h4>' + table(['Group', 'Metric', 'Value'], lengths)
    breadth = [[r['scope'], r['mapq'], r['threshold'], value(r)] for r in rows if r['metric_id'] == 'coverage.breadth' and not r['target_id'] and r['scope'] in ('targets', 'enrichment', 'off_enrichment')]
    out += '<h4>Coverage breadth</h4>' + table(['Region', 'MAPQ', 'Minimum depth', 'Fraction of reference bases'], breadth)
    adequacy = [[r['threshold'].split(':')[0], r['threshold'].split(':')[1], value(r)] for r in rows if r['metric_id'] == 'adequacy.target_breadth' and not r['target_id'] and r['mapq'] == 20]
    out += '<h4>Theoretical low-VAF sampling adequacy</h4><p>Fraction of exact-target bases meeting 95% binomial sampling probability. Sequencing errors, allele-specific mapping, purity, copy number, molecule dependence and caller filters are not modeled; this is not assay sensitivity or MRD detection.</p>' + table(['Assumed VAF', 'Alternate supports', 'Target breadth'], adequacy)
    features = [[r['panel_id'], r['threshold'], value(r)] for r in rows if r['metric_id'] == 'classifier_features.reference_fraction_at_depth']
    out += '<h4>Classifier CpG feature coverage</h4><p>Full-BAM feature availability does not describe a subsampled native classifier input.</p>' + (table(['Panel', 'Valid depth', 'Fraction of panel'], features) if features else '<p>No matching feature panel or methylation evidence available.</p>')
    out += '<p>Observed enrichment and regional short-read patterns do not establish adaptive-sampling acceptance or rejection rates. Fractions use the denominators recorded in the metric JSON.</p>'
    if data['warnings']:
        out += table(['QC observations'], [[w['message']] for w in data['warnings']])
    return out + '</section>'
