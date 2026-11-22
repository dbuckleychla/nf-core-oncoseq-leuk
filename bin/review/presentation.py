"""Presentation-only summaries; never alter analytical calls or selection."""
import csv
import json
from collections import defaultdict
from html import escape
from pathlib import Path
from adapters import classifiers, chromosome_key


def enrich(data, root, results_root=None):
    root = Path(root)
    native = root / 'classification.native.json'
    if native.exists():
        data['classification'] = classifiers(json.loads(native.read_text()))
    summaries = []
    genes = defaultdict(lambda: [0, 0])
    for key in ('primary', 'mapq60', 'nofilter'):
        rows = [r for r in data.get('coverage', []) if isinstance(r.get(key), (int, float)) and r['end'] > r['start']]
        total = sum(r['end'] - r['start'] for r in rows)
        if total:
            summaries.append(dict(metric=key + ' interval-weighted mean depth', value=round(sum(r[key] * (r['end'] - r['start']) for r in rows) / total, 2)))
    for r in data.get('coverage', []):
        if r.get('gene') and isinstance(r.get('primary'), (int, float)):
            length = r['end'] - r['start']
            genes[r['gene']][0] += length * r['primary']
            genes[r['gene']][1] += length
    data['coverage_genes'] = sorted([dict(gene=g, depth=round(v[0] / v[1], 2)) for g, v in genes.items() if v[1] > 0], key=lambda r: r['depth'])[:20]
    qc = data.get('nasvar', {}).get('qc', {})
    for key in ('reads_aligned', 'reads_on_target', 'target_coverage', 'median_read_quality_enriched'):
        if isinstance(qc.get(key), (int, float)):
            summaries.append(dict(metric='NASVAR ' + key.replace('_', ' '), value=qc[key]))
    cramino = root / 'qc.cramino.tsv'
    seqkit = [('input', root / 'qc.seqkit.tsv')]
    if results_root:
        run = Path(results_root)
        cramino = run / 'program_outputs' / 'cramino' / (data['sample'] + '.txt')
        seqkit = [(scope, run / 'reports' / 'seqkit' / (data['sample'] + '_' + scope + '.tsv')) for scope in ('pass', 'fail')]
    if cramino.exists():
        for line in cramino.read_text().splitlines():
            fields = line.split('\t')
            if len(fields) >= 2 and fields[0] in ('Number of reads', 'Yield [Gb]', 'Mean coverage', 'N50', 'N75', 'Median length', 'Mean length', 'Median identity', 'Mean identity'):
                summaries.append(dict(metric='Cramino ' + fields[0], value=fields[1]))
    for scope, path in seqkit:
        if path.exists():
            with path.open() as stream:
                for row in csv.DictReader(stream, delimiter='\t'):
                    for key in ('num_seqs', 'sum_len', 'N50', 'AvgQual', 'Q20(%)', 'Q30(%)'):
                        if key in row:
                            summaries.append(dict(metric='Basecalling ' + scope + ' ' + key, value=row[key]))
    data['qc_summary'] = summaries
    arms = data.get('nasvar', {}).get('karyotype', {}).get('karyotype', {})
    if arms:
        path = root / 'presentation' / 'nasvar.copy_estimates.svg'
        path.parent.mkdir(exist_ok=True)
        arm_plot(arms, path)
        data['figures'] = [f for f in data.get('figures', []) if f.get('name') != 'NASVAR chromosome-arm copy estimates']
        data['figures'].append(dict(category='copy_number', name='NASVAR chromosome-arm copy estimates', path=str(path.resolve()), reformatted=True))
    coverage = root / 'nasvar.coverage.tsv'
    maf = root / 'nasvar.maf.tsv'
    if results_root:
        coverage = Path(results_root) / 'nasvar' / (data['sample'] + '.coverage.tsv')
        maf = Path(results_root) / 'nasvar' / (data['sample'] + '.maf')
    if coverage.exists() and maf.exists():
        plot = root / 'presentation' / 'nasvar.overview.svg'
        plot.parent.mkdir(exist_ok=True)
        genome_plot(coverage, maf, plot)
        data['figures'] = [f for f in data.get('figures', []) if f.get('name') != 'NASVAR genome coverage and folded allele fraction']
        data['figures'].append(dict(category='copy_number', name='NASVAR genome coverage and folded allele fraction', path=str(plot.resolve()), reformatted=True))


def genome_plot(coverage, maf, output):
    """Display native binned counts and folded allele fractions, not new CN calls."""
    cov, alleles, lengths = [], defaultdict(lambda: [0.0, 0]), {}
    with open(coverage) as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            chrom = row['chromosome']
            if chromosome_key(chrom)[0] > 24:
                continue
            start, end = int(row['bin_start']), int(row['bin_end'])
            value = float(row.get('n_reads_gc_adj') or row['n_reads'])
            cov.append((chrom, (start + end) / 2, value))
            lengths[chrom] = max(lengths.get(chrom, 0), end)
    with open(maf) as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            chrom = row['chrom']; pos = int(row['pos'])
            ref, alt = int(row['ref']), int(row['alt'])
            if chrom not in lengths or ref + alt < 20:
                continue
            bucket = (chrom, pos // 1000000)
            alleles[bucket][0] += min(ref, alt) / (ref + alt)
            alleles[bucket][1] += 1
    chromosomes = sorted(lengths, key=chromosome_key)
    offsets = {}; total = 0
    for chrom in chromosomes:
        offsets[chrom] = total; total += lengths[chrom]
    if not total:
        raise ValueError('No canonical NASVAR coverage bins')
    x = lambda chrom, pos: 80 + 1060 * (offsets[chrom] + pos) / total
    values = sorted(v for _, _, v in cov)
    lower = min(0, values[int((len(values)-1)*.01)])
    maximum = max(lower + 1, values[int((len(values)-1)*.99)])
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 560"><rect width="1200" height="560" fill="#111827"/><g fill="#dbeafe" font-family="sans-serif" font-size="14">']
    def text(xp, yp, value):
        parts.append(f'<text x="{xp}" y="{yp}">{escape(str(value))}</text>')
    text(30, 30, 'NASVAR genome overview — native measurements')
    text(30, 60, 'Native GC-adjusted binned coverage signal (not base depth or inferred copy number)')
    text(30, 295, 'Folded allele fraction: min(ref, alt) / total; depth ≥20; 1 Mb mean for display')
    text(12, 90, round(maximum)); text(12, 245, round(lower)); text(35, 330, '0.5'); text(45, 470, 0)
    for chrom in chromosomes:
        xp = x(chrom, 0)
        parts.append(f'<path d="M{xp:.1f} 75V480" stroke="#334155"/>')
        text(x(chrom, lengths[chrom] / 2) - 5, 505, chrom.removeprefix('chr'))
    for chrom, pos, value in cov:
        shown = max(lower, min(maximum, value))
        color = '#f59e0b' if shown != value else '#38bdf8'
        parts.append(f'<circle cx="{x(chrom,pos):.1f}" cy="{245-155*(shown-lower)/(maximum-lower):.1f}" r="1.7" fill="{color}"/>')
    for (chrom, bucket), (value, count) in alleles.items():
        parts.append(f'<circle cx="{x(chrom,(bucket+.5)*1000000):.1f}" cy="{470-280*value/count:.1f}" r="1.5" fill="#c084fc"/>')
    text(30, 528, 'Coverage display: 1st–99th percentile range extended to zero; amber points are clipped outliers. Native values are unchanged.')
    text(30, 549, 'Chromosome widths: observed coverage-bin extents. Display aggregation does not change native calls.')
    parts.append('</g></svg>')
    Path(output).write_text(''.join(parts))


def arm_plot(arms, output):
    import re
    ordered = sorted(((k, v) for k, v in arms.items() if isinstance(v, (int, float))),
                     key=lambda r: (chromosome_key(re.sub('[pq]$', '', r[0])), r[0][-1]))
    maximum = max(4, max((v for _, v in ordered), default=4))
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 360"><rect width="1200" height="360" fill="#111827"/><g fill="#dbeafe" font-family="sans-serif" font-size="13"><text x="30" y="28">NASVAR native chromosome-arm copy estimates</text>']
    y = lambda v: 285 - 230 * v / maximum
    parts.append(f'<path d="M60 {y(2):.1f}H1170" stroke="#f59e0b" stroke-dasharray="6 4"/><text x="60" y="345">Amber line: two copies. Native estimates are retained without reinterpretation.</text>')
    for index, (arm, copies) in enumerate(ordered):
        x = 75 + index * 1090 / max(1, len(ordered))
        width = 850 / max(1, len(ordered))
        parts.append(f'<rect x="{x:.1f}" y="{y(copies):.1f}" width="{width:.1f}" height="{max(0,285-y(copies)):.1f}" fill="#38bdf8"/><text x="{x:.1f}" y="305" transform="rotate(45 {x:.1f} 305)">{escape(arm)}</text>')
    for value in range(int(maximum)+1):
        parts.append(f'<text x="30" y="{y(value):.1f}">{value}</text>')
    parts.append('</g></svg>')
    Path(output).write_text(''.join(parts))
