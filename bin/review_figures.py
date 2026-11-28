#!/usr/bin/env python3
"""Generate readable, event-specific fusion panels using the pinned Figeno tool."""
import argparse
import csv
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def build(input, bam, output, jobs=4):
    data = json.loads(Path(input).read_text())
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    import pysam
    with pysam.AlignmentFile(bam,'rb') as alignments:
        lengths=dict(zip(alignments.references,alignments.lengths))
    selected = [e for e in data['structural_events'] if e['selected']]

    def draw(event):
        row = event['evidence'][0]
        endpoints = [(row['chrom'], row['pos']), (row['partner']['chrom'], row['partner']['pos'])]
        if any(c is None or not isinstance(p, int) for c, p in endpoints):
            raise ValueError('Selected event has unresolved plot endpoints')
        prefix = out / event['id']
        sv = prefix.with_suffix('.tsv')
        with sv.open('w') as stream:
            writer = csv.writer(stream, delimiter='\t')
            writer.writerow(['chr1', 'pos1', 'chr2', 'pos2', 'strand1', 'strand2', 'color', 'svtype'])
            writer.writerow([endpoints[0][0], endpoints[0][1], endpoints[1][0], endpoints[1][1],
                             '', '', '#086f83', row.get('evidence', {}).get('info', {}).get('SVTYPE', 'BND')])
        common = dict(margin_above=2, bounding_box=True, fontscale=1.3, label_rotate=False)
        config = dict(
            general=dict(layout='horizontal', reference='hg38'),
            output=dict(file=str(prefix.with_suffix('.png')), dpi=200, width=280),
            regions=[dict(chr=c, start=max(0, p - 5000), end=min(lengths[c],p + 5000), color=color)
                     for (c, p), color in zip(endpoints, ['#086f83', '#b16832'])],
            highlights=[],
            tracks=[
                dict(common, type='sv', height=20, label='Junction', file=str(sv)),
                dict(common, type='coverage', height=25, label='Coverage', file=str(bam),
                     n_bins=500, scale='auto', color='#748a96'),
                dict(common, type='alignments', height=45, label='Split reads', file=str(bam),
                     only_show_splitreads=True, link_splitreads=True, group_by='none',
                     color_by='none', read_color='#709eb2', splitread_color='#086f83',
                     link_color='#955715', min_splitreads_breakpoints=1),
                dict(common, type='genes', height=20, label='Genes', style='default',
                     collapsed=True, genes='auto', show_gene_names=True, exon_color='#086f83'),
                dict(common, type='chr_axis', height=12, label='Genomic position', unit='kb'),
            ],
        )
        path = prefix.with_suffix('.json')
        path.write_text(json.dumps(config))
        subprocess.run(['figeno', 'make', str(path)], check=True)
        pairs = '; '.join('::'.join(p) for p in event['gene_pairs'])
        return dict(category='fusions', path=prefix.name + '.png', generated=True,
                    event_id=event['id'], name=pairs + ' | ' + ', '.join(event['callers']) + ' | ' + event['id'])

    with ThreadPoolExecutor(max_workers=max(1, min(4, jobs))) as pool:
        figures = list(pool.map(draw, selected))
    (out / 'figures.json').write_text(json.dumps(figures, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('input', 'bam', 'output'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--jobs', type=int, default=4)
    build(**vars(parser.parse_args()))
