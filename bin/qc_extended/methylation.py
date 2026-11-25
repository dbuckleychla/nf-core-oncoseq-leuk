"""CpG dyad summaries; valid depth counted once across modification rows."""
from collections import Counter
import csv
import gzip
import itertools
from pathlib import Path

import pysam

from .contracts import checksum, identifier, state
from .regions import Regions, indexes, primary_contig, read_bed, text_open

SCOPES = ('genome_primary', 'other_contigs', 'targets', 'enrichment', 'off_enrichment')


def panels(manifest, fasta, lengths, assembly):
    result = []
    if not manifest:
        return result
    with open(manifest) as source, pysam.FastaFile(fasta) as ref:
        reader = csv.DictReader(source, delimiter='\t')
        if reader.fieldnames != ['panel_id', 'assembly', 'model_name', 'model_version', 'bed']:
            raise ValueError('Invalid classifier CpG manifest columns')
        seen = set()
        for row in reader:
            key = identifier(row['panel_id'])
            if key in seen or row['assembly'] != assembly:
                raise ValueError('Duplicate panel or incompatible assembly')
            seen.add(key)
            path = Path(row['bed'])
            path = path if path.is_absolute() else Path(manifest).parent / path
            coords = read_bed(path, lengths, features=True)
            features = {(c, a) for c, rows in coords.items() for a, _b, _key, _label in rows}
            if any(ref.fetch(c, a, a + 2).upper() != 'CG' for c, a in features):
                raise ValueError('Panel feature is not a reference CpG C')
            result.append(dict(panel_id=key, assembly=assembly, model_name=row['model_name'],
                               model_version=row['model_version'], sha256=checksum(path), features=features,
                               input_rows=sum(len(r) for r in coords.values())))
    return result


def collect(bed, reference, cpgs, enrichment, targets, supported, feature_panels=()):
    es, ts = indexes(read_bed(enrichment, reference['lengths'])), indexes(read_bed(targets, reference['lengths']))
    empty = Regions([])

    def scopes(chrom, pos):
        if not primary_contig(chrom):
            return ['other_contigs']
        result = ['genome_primary', 'enrichment' if es.get(chrom, empty).contains(pos) else 'off_enrichment']
        if ts.get(chrom, empty).contains(pos):
            result.append('targets')
        return result

    denominators = Counter()
    with gzip.open(cpgs, 'rt') as source:
        for line in source:
            chrom, pos = line.split()
            for scope in scopes(chrom, int(pos)):
                denominators[scope] += 1
    metrics = {s: {m: dict(sites=0, valid=0, modified=0, beta_sum=0., depth=Counter(), beta=Counter(), beta10=Counter(),
                          filter_fail=0, no_call=0, deletion=0, differing_base=0)
                   for m in ['combined', '5mC', '5hmC']} for s in SCOPES}
    feature_depths = {p['panel_id']: Counter() for p in feature_panels}
    coordinate_order = {c: i for i, c in enumerate(reference['lengths'])}
    previous = None
    with text_open(bed) as source:
        lines = (line.split() for line in source if line.strip() and not line.startswith('#'))
        for (chrom, pos), group in itertools.groupby(lines, lambda f: (f[0], int(f[1]))):
            order = (coordinate_order.get(chrom, -1), pos)
            if order[0] < 0 or not 0 <= pos < reference['lengths'][chrom] or previous is not None and order <= previous:
                raise ValueError('bedMethyl is unsorted or incompatible with reference')
            previous, records = order, {}
            for row in group:
                if len(row) != 18 or row[5] != '.' or row[3] not in ('m', 'h') or row[3] in records:
                    raise ValueError('Expected unique strand-combined 18-column m/h records')
                counts = [int(row[i]) for i in [9, 11, 12, 13, 14, 15, 16, 17]]
                if min(counts) < 0 or counts[0] != sum(counts[1:4]):
                    raise ValueError('Inconsistent bedMethyl counts')
                records[row[3]] = counts
            first = next(iter(records.values()))
            valid, modified, canonical, other, deletion, fail, different, nocall = first
            if any(r[0] != valid or r[2] != canonical or r[4:] != first[4:] for r in records.values()):
                raise ValueError('Inconsistent modification denominators')
            if len(records) == 2 and records['m'][1] + records['h'][1] != valid - canonical:
                raise ValueError('Inconsistent combined modification counts')
            values = {'combined': valid - canonical}
            values.update({'5mC' if code == 'm' else '5hmC': row[1] for code, row in records.items() if code in supported})
            for scope in scopes(chrom, pos):
                for mod, count in values.items():
                    d = metrics[scope][mod]
                    d['filter_fail'] += fail
                    d['no_call'] += nocall
                    d['deletion'] += deletion
                    d['differing_base'] += different
                    if valid:
                        beta = count / valid
                        d['sites'] += 1
                        d['valid'] += valid
                        d['modified'] += count
                        d['beta_sum'] += beta
                        d['depth'][valid] += 1
                        d['beta'][min(99, int(beta * 100))] += 1
                        if valid >= 10:
                            d['beta10'][min(99, int(beta * 100))] += 1
            for panel in feature_panels:
                if (chrom, pos) in panel['features']:
                    feature_depths[panel['panel_id']][valid] += 1
    for scope, modifications in metrics.items():
        for mod, data in modifications.items():
            data['reference_cpgs'] = denominators[scope]
            data['available'] = mod == 'combined' or ('m' if mod == '5mC' else 'h') in supported
    features = []
    for panel in feature_panels:
        hist = feature_depths[panel['panel_id']]
        size = len(panel['features'])
        hist[0] += size - sum(hist.values())
        features.append(dict(**{k: v for k, v in panel.items() if k != 'features'}, size=size, depth=dict(hist)))
    return dict(status=state(), scopes=metrics, features=features,
                filtering=dict(percentile=.1, seed=42, sampled_reads=10042, max_depth=8000,
                               policy='Pinned modkit primary-alignment policy; not MAPQ20 coverage',
                               thresholds_status='threshold_values_unavailable'))
