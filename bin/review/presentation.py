"""Presentation-only summaries; never alter analytical calls or selection."""
import csv
import json
from collections import defaultdict
from html import escape
from pathlib import Path
from adapters import classifiers, chromosome_key


def enrich(data, root, results_root=None):
    from variant_presentation import enrich_variants
    enrich_variants(data)
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
            summaries.append(dict(metric='Mean target depth (' + {'primary': 'primary alignments', 'mapq60': 'high-confidence alignments', 'nofilter': 'all alignments'}[key] + ')', value=round(sum(r[key] * (r['end'] - r['start']) for r in rows) / total, 2)))
    for r in data.get('coverage', []):
        if r.get('gene') and isinstance(r.get('primary'), (int, float)):
            length = r['end'] - r['start']
            genes[r['gene']][0] += length * r['primary']
            genes[r['gene']][1] += length
    data['coverage_genes'] = sorted([dict(gene=g, depth=round(v[0] / v[1], 2)) for g, v in genes.items() if v[1] > 0], key=lambda r: r['depth'])[:20]
    qc = data.get('nasvar', {}).get('qc', {})
    for key in ('reads_aligned', 'reads_on_target', 'target_coverage', 'median_read_quality_enriched'):
        if isinstance(qc.get(key), (int, float)):
            summaries.append(dict(metric={'reads_aligned':'Aligned reads', 'reads_on_target':'On-target reads', 'target_coverage':'Mean target depth (NASVAR)', 'median_read_quality_enriched':'Median read quality (NASVAR)'}[key], value=qc[key]))
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
                            summaries.append(dict(metric=scope.capitalize() + ' reads: ' + {'num_seqs':'count', 'sum_len':'total bases', 'N50':'N50 length (bases)', 'AvgQual':'mean Phred quality', 'Q20(%)':'bases at Q20 (%)', 'Q30(%)':'bases at Q30 (%)'}[key], value=row[key]))
    data['qc_summary'] = summaries
    data['qc_plots'] = qc_plots(data, seqkit)
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
        baf = genome_plot(coverage, maf, plot)
        data['figures'] = [f for f in data.get('figures', []) if f.get('name') not in ('NASVAR genome coverage and folded allele fraction', 'GC-corrected coverage', 'B-allele frequency (BAF)')]
        data['figures'].append(dict(category='copy_number', name='GC-corrected coverage', path=str(plot.resolve()), reformatted=True, presentation_order=1))
        data['figures'].append(dict(category='copy_number', name='B-allele frequency (BAF)', path=str(baf.resolve()), reformatted=True, presentation_order=2))


def genome_plot(coverage, maf, output):
    """Separate raw GC-corrected counts and alternate allele fractions."""
    import math
    cov, alleles, lengths = [], defaultdict(int), {}
    with open(coverage) as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            chrom = 'chr' + row['chromosome'].removeprefix('chr')
            if chromosome_key(chrom)[0] > 24:
                continue
            start, end = int(row['bin_start']), int(row['bin_end'])
            value = float(row.get('n_reads_gc_adj') or row['n_reads'])
            cov.append((chrom, (start + end) / 2, value))
            lengths[chrom] = max(lengths.get(chrom, 0), end)
    with open(maf) as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            chrom = 'chr' + row['chrom'].removeprefix('chr')
            ref, alt = int(row['ref']), int(row['alt'])
            if chrom not in lengths or ref + alt < 20:
                continue
            fraction = alt / (ref + alt)
            alleles[(chrom, int(row['pos']) // 1000000, min(99, int(fraction * 100)))] += 1
    chromosomes = sorted(lengths, key=chromosome_key)
    offsets = {}; total = 0
    for chrom in chromosomes:
        offsets[chrom] = total; total += lengths[chrom]
    if not total:
        raise ValueError('No canonical NASVAR coverage bins')
    x = lambda chrom, pos: 90 + 1050 * (offsets[chrom] + pos) / total
    def canvas(title, unit, ticks, maximum):
        parts = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 420"><rect width="1200" height="420" fill="#111827"/><g fill="#dbeafe" font-family="sans-serif" font-size="14">']
        parts.append(f'<text x="30" y="27">{escape(title)}</text><text x="30" y="52">{escape(unit)}</text>')
        for tick in ticks:
            yp = 330 - 250 * tick / maximum
            parts.append(f'<path d="M85 {yp:.1f}H1145" stroke="#334155"/><text x="80" y="{yp+5:.1f}" text-anchor="end">{tick:,.2f}</text>')
        for chrom in chromosomes:
            xp = x(chrom, 0)
            parts.append(f'<path d="M{xp:.1f} 80V330" stroke="#334155"/><text x="{x(chrom,lengths[chrom]/2)-5:.1f}" y="355">{escape(chrom.removeprefix("chr"))}</text>')
        return parts
    values = sorted(v for _, _, v in cov)
    upper = max(1, values[int((len(values)-1)*.99)])
    magnitude = 10 ** math.floor(math.log10(upper))
    upper = math.ceil(upper / magnitude) * magnitude
    parts = canvas('GC-corrected coverage', 'Read counts per bin — not log ratio or copy number', [upper*i/4 for i in range(5)], upper)
    for chrom, pos, value in cov:
        shown = max(0, min(upper, value))
        color = '#f59e0b' if shown != value else '#38bdf8'
        parts.append(f'<circle cx="{x(chrom,pos):.1f}" cy="{330-250*shown/upper:.1f}" r="1.7" fill="{color}"/>')
    parts.append('<text x="30" y="390">Amber points exceed the displayed coverage range. Original counts are unchanged.</text></g></svg>')
    Path(output).write_text(''.join(parts))
    parts = canvas('B-allele frequency (BAF)', 'Alternate allele reads / total reads; depth ≥20', [0, .25, .5, .75, 1], 1)
    parts.append('<path d="M90 205H1140" stroke="#f59e0b" stroke-dasharray="5 4"/>')
    peak = max(alleles.values(), default=1)
    for (chrom, bucket, fraction_bin), count in alleles.items():
        opacity = .55 + .45 * math.log1p(count) / math.log1p(peak)
        width = max(.7, 1050 * 1000000 / total)
        parts.append(f'<rect x="{x(chrom,bucket*1000000):.1f}" y="{330-250*(fraction_bin+.5)/100:.1f}" width="{width:.2f}" height="2.5" fill="#e9d5ff" opacity="{opacity:.3f}"/>')
    if not alleles:
        parts.append('<text x="300" y="210">No sites meet the depth requirement.</text>')
    parts.append('<text x="30" y="390">Brighter bands contain more sites. Display bins: 1 Mb × 0.01 allele fraction.</text></g></svg>')
    baf = Path(output).with_name('nasvar.baf.svg')
    baf.write_text(''.join(parts))
    return baf


def arm_plot(arms, output):
    import re
    ordered = sorted(((k, v) for k, v in arms.items() if isinstance(v, (int, float))),
                     key=lambda r: (chromosome_key(re.sub('[pq]$', '', r[0])), r[0][-1]))
    maximum = max(4, max((v for _, v in ordered), default=4))
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 360"><rect width="1200" height="360" fill="#111827"/><g fill="#dbeafe" font-family="sans-serif" font-size="13"><text x="30" y="28">Estimated copies by chromosome arm</text>']
    y = lambda v: 285 - 230 * v / maximum
    parts.append(f'<path d="M60 {y(2):.1f}H1170" stroke="#f59e0b" stroke-dasharray="6 4"/><text x="60" y="345">Amber line: two copies. Estimates are reported by NASVAR.</text>')
    for index, (arm, copies) in enumerate(ordered):
        x = 75 + index * 1090 / max(1, len(ordered))
        width = 850 / max(1, len(ordered))
        parts.append(f'<rect x="{x:.1f}" y="{y(copies):.1f}" width="{width:.1f}" height="{max(0,285-y(copies)):.1f}" fill="#38bdf8"/><text x="{x:.1f}" y="305" transform="rotate(45 {x:.1f} 305)">{escape(arm)}</text>')
    for value in range(int(maximum)+1):
        parts.append(f'<text x="30" y="{y(value):.1f}">{value}</text>')
    parts.append('</g></svg>')
    Path(output).write_text(''.join(parts))


def qc_plots(data, seqkit):
    plots = []
    genes = data.get('coverage_genes', [])
    if genes:
        plots.append(dict(title='Lowest-coverage genes', unit='Mean depth', bars=[dict(label=r['gene'], value=r['depth']) for r in genes]))
    depths = sorted(r['primary'] for r in data.get('coverage', []) if isinstance(r.get('primary'), (int, float)) and r['primary'] >= 0)
    if depths:
        limit = max(1, depths[int((len(depths)-1)*.99)])
        counts = [0] * 12
        for depth in depths:
            counts[min(11, int(depth / limit * 12))] += 1
        plots.append(dict(title='Target coverage distribution', unit='Number of intervals', bars=[dict(label=(f'≥{11*limit/12:.0f}×' if i == 11 else f'{i*limit/12:.0f}–{(i+1)*limit/12:.0f}×'), value=count) for i, count in enumerate(counts)]))
    for key, title, unit in [('num_seqs', 'Read counts', 'Reads'), ('N50', 'Read length N50', 'Bases'), ('AvgQual', 'Mean read quality', 'Phred score')]:
        bars = []
        for scope, path in seqkit:
            if path.exists():
                with path.open() as stream:
                    for index, row in enumerate(csv.DictReader(stream, delimiter='\t')):
                        try:
                            value = float(row[key])
                        except (KeyError, ValueError):
                            continue
                        bars.append(dict(label=scope.capitalize() + (f' {index+1}' if index else ''), value=value))
        if bars:
            plots.append(dict(title=title, unit=unit, bars=bars))
    return plots
