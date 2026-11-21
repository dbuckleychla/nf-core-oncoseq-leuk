"""Exact COSMIC allele annotation with transcript-aware recurrence provenance."""
from collections import defaultdict

from annotate_small_variants import exact_records


def annotate(handles, record):
    sources = []
    for name, handle in sorted(handles.items()):
        if not name.startswith('cosmic'):
            continue
        observations = {}
        counts = defaultdict(set)
        for hit in exact_records(handle, record):
            ids = tuple(sorted((hit.id or '').split(';')))
            info = dict(hit.info)
            for identifier in ids:
                for field in ('TARGETED_SCREEN_SAMPLE_COUNT', 'GENOME_SCREEN_SAMPLE_COUNT'):
                    if field in info:
                        counts[(identifier, field)].add(info[field])
            key = (ids, tuple(sorted((k, str(v)) for k, v in info.items())))
            observations[key] = dict(ids=list(ids), info=info)
        sources.append(dict(
            source=name,
            observations=list(observations.values()),
            identifiers=sorted({i for row in observations.values() for i in row['ids']}),
            recurrence=[dict(identifier=i, field=f, values=sorted(v), inconsistent=len(v) > 1)
                        for (i, f), v in sorted(counts.items())],
        ))
    return dict(status='not_configured' if not sources else
                'matched' if any(s['observations'] for s in sources) else 'no_match',
                sources=sources)
