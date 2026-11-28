"""Conservative complete-link consensus of oriented native junctions."""
import hashlib
import json
import math

from report_prioritization import fusion_candidate


def finite_count(value):
    return value if isinstance(value, (int, float)) and math.isfinite(value) and value >= 0 else None


def support(row):
    evidence = row.get('evidence', {})
    if row['caller'] == 'NASVAR':
        return finite_count(evidence.get('supporting_reads'))
    fmt = evidence.get('format', {}).get(row.get('vcf_sample',row['sample']), {})
    info = evidence.get('info', {})
    value = next((finite_count(fmt.get(k)) for k in ('DV', 'SR', 'AD')
                  if finite_count(fmt.get(k)) is not None), None)
    if value is None:
        value = next((finite_count(info.get(k)) for k in ('SUPPORT', 'SR')
                      if finite_count(info.get(k)) is not None), None)
    if row.get('mate_evidence'):
        other = support(row['mate_evidence'])
        value = min(value, other) if value is not None and other is not None else None
    return value


def junction(row):
    partner = row.get('partner', {})
    a, b = partner.get('local_side'), partner.get('remote_side')
    if not row.get('chrom') or not partner.get('chrom') or not isinstance(partner.get('pos'), int):
        return None
    kind = row.get('evidence', {}).get('info', {}).get('SVTYPE', 'BND')
    kind = 'BND' if kind in ('TRA', 'BND', 'native fusion') else kind
    if kind == 'DUP:TANDEM':
        kind = 'DUP'
    if a is None or b is None:
        return None
    left = (row['chrom'], row['pos'] - (1 if row['caller']=='NASVAR' else a == 'left'), a)
    right = (partner['chrom'], partner['pos'] - (1 if row['caller']=='NASVAR' else b == 'left'), b)
    return (kind, *sorted((left, right)))


def compatible(a, b, distance):
    ja, jb = junction(a), junction(b)
    return bool(ja and jb and a['caller'] != b['caller'] and ja[0] == jb[0]
                and all(x[0] == y[0] and x[2] == y[2] and abs(x[1] - y[1]) <= distance
                        for x, y in zip(ja[1:], jb[1:])))


def compatibility_edges(rows, distance):
    """Build the same graph using type/orientation partitions and a position window."""
    from collections import defaultdict
    groups = defaultdict(list)
    edges = {i: set() for i in range(len(rows))}
    for index, row in enumerate(rows):
        value = junction(row)
        if value:
            kind, left, right = value
            groups[(kind, left[0], left[2], right[0], right[2])].append((left[1], right[1], index))
    for values in groups.values():
        values.sort()
        first = 0
        for current, (left, right, index) in enumerate(values):
            while first < current and left - values[first][0] > distance:
                first += 1
            for previous in range(first, current):
                _, other_right, other = values[previous]
                if rows[index]['caller'] != rows[other]['caller'] and abs(right - other_right) <= distance:
                    edges[index].add(other)
                    edges[other].add(index)
    return edges


def consolidate(rows, known_pairs, distance=100):
    if distance < 0:
        raise ValueError('Fusion distance must be nonnegative')
    rows = sorted(rows, key=lambda r: (str(junction(r)), r['caller'], r['id']))
    edges = compatibility_edges(rows, distance)
    seen, events = set(), []
    for i in range(len(rows)):
        if i in seen:
            continue
        component, pending = set(), [i]
        while pending:
            n = pending.pop()
            if n in component:
                continue
            component.add(n)
            pending.extend(edges[n] - component)
        seen.update(component)
        unambiguous = (len({rows[n]['caller'] for n in component}) == len(component)
                       and all(component - {n} <= edges[n] for n in component))
        groups = [sorted(component)] if unambiguous else [[n] for n in sorted(component)]
        for group in groups:
            evidence = [rows[n] for n in group]
            selected = []
            for row in evidence:
                candidate = fusion_candidate(row, known_pairs)
                if candidate:
                    selected.append(candidate)
            identity = [(r['caller'], r['id'], junction(r)) for r in evidence]
            identifier = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:24]
            events.append(dict(
                id=identifier,
                representative=junction(evidence[0]),
                gene_pairs=sorted({tuple(c['pair']) for c in selected}),
                selected=bool(selected),
                selection=selected,
                consensus=len(evidence) > 1,
                ambiguous=not unambiguous,
                possible_matches=[rows[n]['id'] for n in sorted(component) if n not in group],
                callers=sorted({r['caller'] for r in evidence}),
                evidence=[dict(r, native_support=support(r)) for r in evidence],
            ))
    return sorted(events, key=lambda e: (not e['selected'], not e['consensus'], str(e['gene_pairs']), e['id']))
