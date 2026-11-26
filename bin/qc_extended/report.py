"""Portable, presentation-only QC summaries from aggregate metrics."""
from html import escape
from bisect import bisect_right
import math


LABELS = {
    'all': 'All primary reads', 'mapped': 'Mapped reads', 'unmapped': 'Unmapped reads',
    'on_target': 'On exact targets', 'off_target': 'Off exact targets',
    'on_enrichment': 'On enrichment regions', 'off_enrichment': 'Outside enrichment regions',
    'targets': 'Exact targets', 'enrichment': 'Enrichment regions', 'genome_primary': 'Primary genome',
}
COLORS = ('#238b9a', '#d98c39')


def validate(data):
    if data.get('schema_name') != 'oncoseq_extended_qc' or data.get('schema_version') != '1.0.0':
        raise ValueError('Unsupported extended QC schema')
    for row in data['metrics']:
        if row['value'] is not None and (not isinstance(row['value'], (int, float)) or not math.isfinite(row['value'])):
            raise ValueError('Invalid metric value')


def table(headers, rows):
    return '<div class="qc-ext-table"><table><thead><tr>' + ''.join('<th>' + escape(str(h)) + '</th>' for h in headers) + '</tr></thead><tbody>' + ''.join('<tr>' + ''.join('<td>' + escape(str(v)) + '</td>' for v in row) + '</tr>' for row in rows) + '</tbody></table></div>'


def percent(n):
    if n is None:
        return 'Unavailable'
    return '<0.01%' if 0 < n < .0001 else f'{100*n:.2f}%'


def value(row):
    if not row or row['value'] is None:
        return 'Unavailable'
    n = row['value']
    if row['unit'] == 'fraction':
        return percent(n)
    return f'{n:,.1f}'.rstrip('0').rstrip('.') if isinstance(n, float) else f'{n:,}'


def details(title, body):
    return '<details class="qc-ext-details"><summary>' + escape(title) + '</summary>' + body + '</details>'


def box(title, body, note=''):
    return '<div class="qc-ext-box"><h4>' + escape(title) + '</h4>' + ('<p class="qc-ext-note">' + escape(note) + '</p>' if note else '') + body + '</div>'


def rebin(histogram, edges):
    """Merge whole source bins only; preserve counts and disclose overflow.

    Display edges align with the original 100-bp or 0.01-beta grid. No
    interpolation, smoothing, sampling, or renormalization is performed.
    """
    counts = [0] * (len(edges) - 1)
    overflow = 0
    for b in histogram.get('bins', []):
        lo, hi = b['lower'], b['upper']
        i = bisect_right(edges, lo + 1e-9) - 1
        if i == len(counts) and b.get('upper_inclusive') and math.isclose(lo, edges[-1]):
            i -= 1
        if i < 0:
            raise ValueError('Histogram starts below display range')
        if i >= len(counts):
            overflow += b['count']
        else:
            if hi > edges[i + 1] + 1e-8:
                raise ValueError('Display edge splits a source histogram bin')
            counts[i] += b['count']
    return counts, overflow


def bars(title, counts, denominator, ticks, color, peak, axis):
    """Bars use an explicit percent-of-group y axis; ticks are slot boundaries."""
    if not denominator:
        return '<p class="qc-ext-empty">No observations available for this group.</p>'
    left, width, top, height = 52, 440, 25, 150
    svg = '<svg viewBox="0 0 510 245" role="img" aria-label="' + escape(title + '; ' + axis) + '">'
    svg += '<text x="52" y="13" class="qc-ext-axis">% of group</text>'
    for fraction in (0, .5, 1):
        y = top + height * (1 - fraction)
        svg += f'<path d="M{left} {y}h{width}" class="qc-ext-grid"/><text x="44" y="{y+4}" text-anchor="end" class="qc-ext-axis">{peak*fraction*100:g}%</text>'
    slot = width / len(counts)
    for i, (label, count) in enumerate(counts):
        bar_height = count / denominator / peak * height
        svg += f'<rect x="{left+i*slot+1:.2f}" y="{top+height-bar_height:.2f}" width="{max(.5,slot-2):.2f}" height="{bar_height:.2f}" fill="{color}"><title>{escape(label)}: {count:,} observations ({escape(percent(count/denominator))} of group)</title></rect>'
    for position, label in ticks:
        svg += f'<text x="{left+position*slot:.2f}" y="198" text-anchor="middle" class="qc-ext-axis">{escape(label)}</text>'
    svg += f'<text x="272" y="231" text-anchor="middle" class="qc-ext-axis">{escape(axis)}</text></svg>'
    return svg


def peak_scale(groups):
    peak = max((max(counts, default=0) / h['denominator'] for h, counts, _overflow in groups if h.get('denominator')), default=0)
    return max(.05, math.ceil(peak * 20) / 20)


def length_comparison(data, lookup, group):
    scopes = ('on_' + group, 'off_' + group)
    series = [next((h for h in data['histograms'] if h['metric_id'] == 'read_length.distribution' and h['scope'] == scope),
                   dict(scope=scope, bins=[], denominator=0)) for scope in scopes]
    views = [
        ('0–10 kb', list(range(0, 10001, 500)), [(i * 4, str(i * 2)) for i in range(6)], 'Read length (kb); 500-bp bins'),
        ('0–2 kb', list(range(0, 2001, 100)), [(i * 5, f'{i*.5:g}') for i in range(5)], 'Read length (kb); 100-bp bins'),
        ('All lengths', [0, 500, 1000, 2000, 5000, 10000, 20000, 50000, math.inf],
         [(i+.5, label) for i, label in enumerate(['<0.5', '.5–1', '1–2', '2–5', '5–10', '10–20', '20–50', '≥50'])], 'Length bands (kb); unequal widths'),
    ]
    out = '<div class="qc-ext-lengths" role="group" aria-label="Read length display range">'
    for i, (label, _edges, _ticks, _axis) in enumerate(views):
        key = f'qc-ext-{group}-{i}'
        out += f'<input class="qc-ext-range qc-ext-range-{i}" type="radio" name="qc-ext-{group}" id="{key}"' + (' checked' if i == 0 else '') + f'><label for="{key}">{escape(label)}</label>'
    for index, (_label, edges, ticks, axis) in enumerate(views):
        groups = [(h, *rebin(h, edges)) for h in series]
        peak = peak_scale(groups)
        out += f'<div class="qc-ext-view qc-ext-view-{index}"><div class="qc-ext-two">'
        for i, (h, counts, overflow) in enumerate(groups):
            scope, denominator = h['scope'], h.get('denominator', 0)
            title = LABELS[scope]
            summary = f'{denominator:,} measured reads · Median {value(lookup("read_length.median", scope))} bp · N50 {value(lookup("read_length.n50", scope))} bp'
            labels = [f'{lo:g}–{hi:g} bp' if math.isfinite(hi) else f'≥{lo:g} bp' for lo, hi in zip(edges, edges[1:])]
            plot = bars(title, list(zip(labels, counts)), denominator, ticks, COLORS[i], peak, axis)
            if index < 2 and denominator:
                plot += f'<p class="qc-ext-tail">Beyond this view (≥{edges[-1]/1000:g} kb): <strong>{escape(percent(overflow/denominator))}</strong> · {overflow:,} reads</p>'
            out += '<div class="qc-ext-distribution"><h5><span style="background:' + COLORS[i] + '"></span>' + title + '</h5><p class="qc-ext-note">' + escape(summary) + '</p>' + plot + '</div>'
        out += '</div></div>'
    out += '</div><p class="qc-ext-note">Each group sums to 100% of reads with known length. The two plots share a y-axis scale within each view. Longer reads stay in the denominator when zoomed; use All lengths to include the full tail. Missing lengths are excluded. Hover over a bar for counts.</p>'
    return out


def heatmap(headers, labels, values, row_label='Region'):
    out = '<div class="qc-ext-table"><table class="qc-ext-heatmap"><thead><tr><th>' + escape(row_label) + '</th>' + ''.join('<th>' + escape(h) + '</th>' for h in headers) + '</tr></thead><tbody>'
    for label, cells in zip(labels, values):
        out += '<tr><th scope="row">' + escape(label) + '</th>'
        for n in cells:
            style = f' style="background:rgba(35,139,154,{.06+.6*n:.3f})"' if n is not None else ''
            out += '<td' + style + '>' + escape(percent(n)) + '</td>'
        out += '</tr>'
    return out + '</tbody></table></div>'


def content(data):
    validate(data)
    rows = data['metrics']
    index = {(r['metric_id'], r['scope'], r.get('mapq'), r.get('modification'), r.get('threshold')): r for r in rows if not r.get('target_id') and not r.get('panel_id')}

    def lookup(metric, scope='all', mapq=None, modification=None, threshold=None):
        return index.get((metric, scope, mapq, modification, threshold))

    def number(*args, **kwargs):
        row = lookup(*args, **kwargs)
        return row['value'] if row else None

    out = '<section class="qc-ext">' + STYLE + '<h3>Extended sequencing QC</h3><p class="qc-ext-note">Full mapping-stage sample: ' + escape(data['sample_id']) + '. Primary records before downstream filtering or classifier subsampling. These summaries are not clinical pass/fail assessments.</p>'
    cards = [
        ('Mapped reads', lookup('alignment.rate'), value(lookup('alignment.primary_records')) + ' primary records'),
        ('On-target median length', lookup('read_length.median', 'on_target'), 'bp · exact-target overlap'),
        ('Mean exact-target depth', lookup('coverage.mean_depth', 'targets', 20), '× · primary alignments, MAPQ ≥20'),
        ('Target CpGs at ≥10 calls', lookup('methylation.reference_fraction_at_depth', 'targets', modification='combined', threshold=10), 'of all reference CpGs in exact targets'),
    ]
    out += '<div class="qc-ext-cards">' + ''.join('<div class="qc-ext-card"><span>' + escape(title) + '</span><strong>' + escape(value(row)) + '</strong><small>' + escape(note) + '</small></div>' for title, row, note in cards) + '</div>'
    statuses = []
    for key, state in data['components'].items():
        statuses.append('<span class="qc-ext-status"><strong>' + escape(key.replace('_', ' ').capitalize()) + '</strong> · ' + escape(state['status']) + (': ' + escape(state['reason'].replace('_', ' ')) if state.get('reason') else '') + '</span>')
    out += '<div class="qc-ext-statuses">' + ''.join(statuses) + '</div>'
    if data['warnings']:
        out += box('QC observations', '<ul>' + ''.join('<li>' + escape(w['message']) + '</li>' for w in data['warnings']) + '</ul>')

    lengths = length_comparison(data, lookup, 'target')
    lengths += details('Compare read lengths on/off enrichment regions', '<h5>Read lengths on/off enrichment</h5>' + length_comparison(data, lookup, 'enrichment'))
    statistics = [[LABELS.get(scope, scope)] + [value(lookup('read_length.' + key, scope)) for key in ('records', 'missing_length_records', 'mean', 'median', 'n50')] for scope in LABELS if scope in ('all', 'mapped', 'unmapped', 'on_target', 'off_target', 'on_enrichment', 'off_enrichment')]
    lengths += details('All read-length statistics', table(['Read group', 'Reads', 'Missing length', 'Mean (bp)', 'Median (bp)', 'N50 (bp)'], statistics))
    out += box('Read lengths on/off exact targets', lengths, 'An on-target read has an aligned match overlapping the exact-target BED. Off-target reads are mapped reads without that overlap. Unmapped reads are excluded from this comparison.')

    scopes = ('targets', 'enrichment', 'off_enrichment')
    thresholds = (1, 5, 10, 20, 30, 50, 100)
    coverage = ''
    for mq in (20, 0):
        grid = heatmap([f'≥{t}×' for t in thresholds], [LABELS[s] for s in scopes], [[number('coverage.breadth', s, mq, threshold=t) for t in thresholds] for s in scopes])
        grid += '<p class="qc-ext-note">Darker cells mean more reference bases meet the depth threshold; color does not indicate pass/fail. Zero-coverage bases are included.</p>'
        coverage += grid if mq == 20 else details('Coverage including MAPQ 0 alignments', grid)
    out += box('Coverage breadth', coverage, 'Percent of reference bases with at least the indicated depth. Default: primary alignments with MAPQ ≥20. Exact targets and enrichment regions can overlap; these rows are not additive.')

    meth_status = data['components'].get('methylation', {})
    if meth_status.get('status') in ('complete', 'partial'):
        methylation = heatmap([f'≥{t} calls' for t in (1, 5, 10, 20)], [LABELS[s] for s in scopes], [[number('methylation.reference_fraction_at_depth', s, modification='combined', threshold=t) for t in (1, 5, 10, 20)] for s in scopes])
        methylation += '<p class="qc-ext-note">Percent of all reference CpG dyads with sufficient valid calls. 5mC and 5hmC share the combined valid-call denominator.</p>'
        beta = []
        for mod in ('5mC', '5hmC'):
            h = next((h for h in data['histograms'] if h['metric_id'] == 'methylation.beta_depth10' and h['scope'] == 'genome_primary' and h['modification'] == mod), None)
            if h and h.get('denominator'):
                counts, overflow = rebin(h, [i/20 for i in range(21)])
                beta.append((h, counts, overflow))
        if beta:
            plots = '<div class="qc-ext-two">'
            peak = peak_scale(beta)
            for i, (h, counts, _overflow) in enumerate(beta):
                labels = [(f'{j/20:.2f}–{(j+1)/20:.2f} beta', n) for j, n in enumerate(counts)]
                plots += '<div><h5>' + escape(h['modification']) + '</h5>' + bars(h['modification'] + ' beta at ≥10 valid calls', labels, h['denominator'], [(j*5, f'{j/4:g}') for j in range(5)], COLORS[i % 2], peak, 'Beta (modified / valid calls)') + f'<p class="qc-ext-note">{h["denominator"]:,} observed CpGs · primary genome · ≥10 valid calls</p></div>'
            methylation += details('Methylation distributions: CpGs with ≥10 valid calls', plots + '</div><p class="qc-ext-note">Low-depth CpGs are excluded from this display. Beta distributions describe observed sites, not the reference CpG coverage shown above.</p>')
        out += box('Methylation QC', methylation)
    else:
        out += box('Methylation QC', '<p>' + escape((meth_status.get('status') or 'unavailable').capitalize()) + ' · ' + escape((meth_status.get('reason') or 'No methylation evidence').replace('_', ' ')) + '</p>')

    adequacy = [r for r in rows if r['metric_id'] == 'adequacy.target_breadth' and not r.get('target_id') and r['mapq'] == 20]
    vafs = sorted({float(r['threshold'].split(':')[0]) for r in adequacy})
    adequacy_grid = heatmap(['≥3 alternate reads', '≥5 alternate reads'], [f'{v*100:g}% VAF' for v in vafs], [[next((r['value'] for r in adequacy if r['threshold'] == f'{v}:{support}'), None) for support in (3, 5)] for v in vafs], 'Assumed VAF')
    out += details('Theoretical low-VAF sampling adequacy', '<p>Percent of exact-target bases meeting 95% binomial sampling probability at MAPQ ≥20. Sequencing errors, allele-specific mapping, purity, copy number, molecule dependence and caller filters are not modeled. This is not assay sensitivity or MRD detection.</p>' + adequacy_grid)
    features = [[r['panel_id'], r['threshold'], value(r)] for r in rows if r['metric_id'] == 'classifier_features.reference_fraction_at_depth']
    out += details('Classifier CpG feature coverage', '<p>Full-BAM feature availability does not describe a subsampled native classifier input.</p>' + (table(['Panel', 'Minimum valid calls', 'Fraction of reference panel'], features) if features else '<p>No matching feature panel or methylation evidence available.</p>'))
    chosen = {'alignment.rate', 'alignment.primary_records', 'coverage.mean_depth', 'coverage.observed_enrichment', 'allocation.target_aligned_base_fraction', 'allocation.padding_aligned_base_fraction', 'methylation.callable_fraction', 'methylation.mean_beta', 'methylation.weighted_beta'}
    summary = [[r['metric_id'].replace('.', ': ').replace('_', ' '), LABELS.get(r['scope'], r['scope']), r['mapq'] if r['mapq'] is not None else '—', r['modification'] or '—', value(r)] for r in rows if r['metric_id'] in chosen and not r['target_id'] and r['scope'] in ('all', 'targets', 'enrichment', 'off_enrichment', 'genome_primary')]
    out += details('Detailed regional metrics and component status', table(['Metric', 'Region', 'Minimum MAPQ', 'Modification', 'Value'], summary) + table(['Component', 'Status', 'Reason'], [[k, s['status'], s.get('reason') or '—'] for k, s in data['components'].items()]))
    out += '<p class="qc-ext-note">Observed enrichment and regional short-read patterns do not establish adaptive-sampling acceptance or rejection rates. Exact values and denominators remain in the QC metric JSON and TSV files.</p>'
    return out + '</section>'


STYLE = '''<style>
.qc-ext{--qc-line:var(--line,#cbd5df);--qc-surface:var(--surface,#f7f9fb);font:inherit;color:inherit}
.qc-ext h4{font-size:1.13em;margin:0 0 12px}.qc-ext h5{font-size:1em;margin:16px 0 4px}
.qc-ext p{line-height:1.5}.qc-ext-note,.qc-ext-tail{font-size:.88em;line-height:1.5;opacity:.85}
.qc-ext-box{border:1px solid var(--qc-line);border-radius:10px;padding:22px;margin:20px 0}
.qc-ext-cards{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px;margin:20px 0}
.qc-ext-card{border:1px solid var(--qc-line);background:var(--qc-surface);border-radius:9px;padding:18px}
.qc-ext-card span,.qc-ext-card small{display:block;font-size:.84em}.qc-ext-card strong{display:block;font-size:1.85em;line-height:1.4;margin:6px 0;font-variant-numeric:tabular-nums}.qc-ext-card small{opacity:.8}
.qc-ext-statuses{display:flex;flex-wrap:wrap;gap:8px}.qc-ext-status{font-size:.8em;border:1px solid var(--qc-line);padding:5px 9px;border-radius:6px}
.qc-ext-two{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:24px}
.qc-ext svg{display:block;width:100%;height:auto;color:inherit;overflow:visible}.qc-ext svg text{fill:currentColor;font-family:inherit}.qc-ext-axis{font-size:12px}.qc-ext-grid{fill:none;stroke:currentColor;opacity:.17}
.qc-ext-distribution h5 span{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:7px}.qc-ext-tail{margin:0}
.qc-ext-table{overflow:auto}.qc-ext table{border-collapse:collapse;width:100%;margin:12px 0;font-size:.9em}.qc-ext th,.qc-ext td{padding:10px;border-bottom:1px solid var(--qc-line);text-align:left;font-variant-numeric:tabular-nums}
.qc-ext .qc-ext-heatmap td{text-align:center;white-space:nowrap;border:3px solid var(--qc-surface);min-width:72px}.qc-ext .qc-ext-heatmap thead th{text-align:center}.qc-ext .qc-ext-heatmap th:first-child{text-align:left}
.qc-ext-details{border-top:1px solid var(--qc-line);margin-top:14px;padding:12px 0}.qc-ext-details>summary{cursor:pointer;font-weight:600;white-space:normal}.qc-ext-details[open]>summary{margin-bottom:16px}
.qc-ext-lengths{margin-top:18px}.qc-ext-range{position:absolute;opacity:0;width:1px;height:1px}.qc-ext-range+label{display:inline-block;padding:7px 14px;border:1px solid var(--qc-line);cursor:pointer;font-size:.88em;margin:0 6px 8px 0;border-radius:6px}
.qc-ext-range:checked+label{background:#238b9a;color:#fff;border-color:#238b9a}.qc-ext-range:focus-visible+label{outline:3px solid #d98c39;outline-offset:2px}.qc-ext-view{display:none}
.qc-ext-range-0:checked~.qc-ext-view-0,.qc-ext-range-1:checked~.qc-ext-view-1,.qc-ext-range-2:checked~.qc-ext-view-2{display:block}
@media(max-width:800px){.qc-ext-cards{grid-template-columns:repeat(2,minmax(0,1fr))}.qc-ext-two{grid-template-columns:1fr}.qc-ext-box{padding:14px}.qc-ext-axis{font-size:16px}}
@media print{.qc-ext-box{break-inside:avoid}.qc-ext-range+label{display:none}}
</style>'''
