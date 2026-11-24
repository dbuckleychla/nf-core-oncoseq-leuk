"""Exact COSMIC allele annotation with targeted-screen and phenotype evidence."""
from collections import defaultdict
import re

from annotate_small_variants import exact_records


def annotate(handles, record, evidence):
    observations = {}
    counts = defaultdict(set)
    for hit in exact_records(handles['cosmic_targeted'], record):
        ids = tuple(sorted(set((hit.id or '').split(';'))))
        if any(not re.fullmatch(r'COSV\d+', i) for i in ids):
            raise ValueError('COSMIC allele is missing a valid COSV identifier')
        info = dict(hit.info)
        for identifier in ids:
            value = info.get('TARGETED_SCREEN_SAMPLE_COUNT')
            if value is not None:
                counts[identifier].add(value)
        key = (ids, tuple(sorted((k, str(v)) for k, v in info.items())))
        observations[key] = dict(ids=list(ids), info=info)
    identifiers = sorted({i for row in observations.values() for i in row['ids']})
    source = dict(source='cosmic_targeted', observations=list(observations.values()),
                  identifiers=identifiers,
                  recurrence=[dict(identifier=i, field='TARGETED_SCREEN_SAMPLE_COUNT',
                                   values=sorted(v), inconsistent=len(v) > 1)
                              for i, v in sorted(counts.items())])
    targeted = [evidence.lookup(identifier) for identifier in identifiers]
    return dict(status='matched' if targeted else 'no_match', sources=[source],
                targeted_evidence=targeted,
                blood_match=any(t['blood_match'] for t in targeted),
                warnings=sorted({t['status'] for t in targeted if t['status'] != 'matched'} |
                                {p['status'] for t in targeted for p in t['phenotypes'] if p['status'] != 'matched'}))
