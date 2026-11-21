"""Native classification, coverage and NASVAR report adapters."""
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def classifiers(data):
    result = []
    order = ('marlin', 'alma', 'lamprey', 'crossnn', 'mpact', 'sturgeon', 'tucan', 'nanomix')
    for name in order:
        root = data.get(name)
        if not isinstance(root, dict):
            continue
        models = root.get('models', [root])
        for model in models:
            inference = model.get('inference', {})
            values = next((v for v in (model.get('probabilities'),
                                       inference.get('annotated_probabilities'),
                                       model.get('predictions'), model.get('votes'))
                           if isinstance(v, list) and v and isinstance(v[0], dict)), [])
            scores = []
            for row in values:
                score = row.get('probability', row.get('score'))
                label = row.get('class', row.get('class_name', row.get('label')))
                if isinstance(score, (int, float)) and label:
                    scores.append(dict(label=label, score=score))
            scores.sort(key=lambda r: (-r['score'], r['label']))
            top = model.get('top_class') or (scores[0]['label'] if scores else None)
            score = model.get('diagnostic_confidence', scores[0]['score'] if scores else None)
            thresholds = {k:model[k] for k in ('diagnostic_threshold','threshold','cutoff') if isinstance(model.get(k),(int,float))}
            cutoff = thresholds.get('diagnostic_threshold',next(iter(thresholds.values())) if len(thresholds)==1 else None)
            result.append(dict(name=name + (':' + model['training_set'] if model.get('training_set') else ''),
                               status='completed', label=top, score=score, cutoff=cutoff, thresholds=thresholds,
                               scores=scores, target_coverage=model.get('target_coverage'),
                               model=model.get('model_label', model.get('model_path', model.get('model_root')))))
    result.append(dict(name='optic', status='disabled', reason='Model assets unavailable in pinned image'))
    return result


def coverage(inputs):
    regions = {}
    for kind, path in inputs:
        label = kind.split(':', 1)[1]
        with open_text(path) as stream:
            for line in stream:
                if line.startswith('#') or not line.strip():
                    continue
                f = line.rstrip().split('\t')
                if len(f) < 4:
                    raise ValueError('Malformed coverage interval')
                start, end = int(f[1]), int(f[2])
                gene = f[3] if len(f) >= 5 else ''
                key = (f[0], start, end, gene)
                row = regions.setdefault(key, dict(chrom=f[0], start=start, end=end, gene=gene))
                if label in row:
                    raise ValueError('Duplicate coverage interval: ' + str(key))
                row[label] = float(f[-1])
    return sorted(regions.values(), key=lambda r: (r.get('primary', float('inf')), r['gene'], r['chrom'], r['start']))


def nasvar_events(data, sample, gene_database=None):
    section = data.get('fusions', {})
    if not isinstance(section, dict) or not isinstance(section.get('fusions'), list):
        return []
    result = []
    for ordinal, event in enumerate(section['fusions']):
        a, b = event.get('gene1', {}), event.get('gene2', {})
        clusters = event.get('breakpoints') or [None]
        for cluster_index, cluster in enumerate(clusters):
            p = cluster or {}
            side_a, side_b = p.get('gene0_dir'), p.get('gene1_dir')
            # Native left traversals end at te; right traversals start at ts.
            native_id = f'NASVAR:{ordinal}:{cluster_index}'
            result.append(dict(
                id=native_id, caller='NASVAR', sample=sample, native_id=native_id,
                chrom=p.get('gene0_chr', a.get('chr')), pos=p.get('gene0_pos', a.get('pos')),
                partner=dict(chrom=p.get('gene1_chr', b.get('chr')),
                             pos=p.get('gene1_pos', b.get('pos')), local_side={'left':'right','right':'left'}.get(side_a), remote_side={'left':'right','right':'left'}.get(side_b)),
                coordinate_convention='Native NASVAR; left boundary end-exclusive, right boundary inclusive',
                native_directions=[side_a, side_b],
                gene_pair=sorted([a['name'], b['name']]) if a.get('name') and b.get('name') else None,
                native_gene_pair=[a.get('name'), b.get('name')],
                genes1=[dict(name=a['name'])] if a.get('name') else [],
                genes2=[dict(name=b['name'])] if b.get('name') else [],
                filter=None, mate_status='native breakpoint cluster',
                evidence=dict(event, supporting_reads=p.get('n_reads', event.get('supporting_reads'))),
                source_records=[event],
            ))
    if gene_database:
        import sqlite3
        from breakpoint_annotation import overlap
        db = sqlite3.connect(gene_database[0])
        aliases = gene_database[1]['contig_mapping']
        for row in result:
            for index, (chrom, pos) in enumerate([(row['chrom'], row['pos']),
                                                (row['partner']['chrom'], row['partner']['pos'])]):
                direction = row['native_directions'][index]
                if chrom in aliases and isinstance(pos, int) and direction in ('left', 'right'):
                    hits = overlap(db, aliases[chrom], pos - (direction == 'left'))
                    row['genes' + str(index + 1)] = hits
            row['annotation_status'] = 'shared GFF breakpoint overlap'
        db.close()
    return result


def chromosome_key(value):
    name = value.removeprefix('chr')
    if name.isdigit():
        return int(name), ''
    return {'X': 23, 'Y': 24, 'M': 25, 'MT': 25}.get(name, 26), name


def interval_table(path):
    result = []
    with open_text(path) as stream:
        for line in stream:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            try:
                start, end = int(fields[1]), int(fields[2])
            except (IndexError, ValueError):
                continue
            result.append(dict(chrom=fields[0], start=start, end=end, native_values=fields[3:]))
    return sorted(result, key=lambda r: (chromosome_key(r['chrom']), r['start'], r['end']))
