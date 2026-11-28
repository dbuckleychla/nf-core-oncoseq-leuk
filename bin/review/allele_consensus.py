"""Derive exact PASS agreement with caller-native depth and AF provenance."""
import json
import sqlite3
from pathlib import Path

import pysam

from report_prioritization import sample_evidence


def add(catalogue):
    root = Path(catalogue)
    sources = json.loads((root / 'sources.json').read_text())
    index = {s['caller']: i for i, s in enumerate(sources)}
    if not {'clairsto', 'deepsomatic'} <= index.keys():
        return
    db = sqlite3.connect(root / 'catalog.sqlite')
    header = pysam.VariantFile(root / 'catalog.vcf').header.copy()
    for caller in ('CS', 'DS'):
        header.info.add(caller + '_DP', 1, 'Integer', 'Copied native caller depth')
        header.info.add(caller + '_AF', 1, 'Float', 'Copied native caller allele fraction')
    header.add_meta('source', value='Exact normalised ClairS-TO and DeepSomatic PASS intersection')
    ordinal = len(sources)
    (root / f'{ordinal}.header.vcf').write_text(str(header))
    for chrom in header.contigs:
        for key, pos, ref, alt in db.execute('SELECT key,pos,ref,alt FROM alleles WHERE chrom=? ORDER BY pos,ref,alt', (chrom,)):
            evidence = {}
            for caller in ('clairsto', 'deepsomatic'):
                records = [text for (text,) in db.execute('SELECT record FROM records WHERE source=? AND key=?', (index[caller], key))
                           if text.split('\t')[6] == 'PASS']
                if records:
                    # Multiple distinct source records remain in ordinary exports.
                    # Derived metrics are copied only when their values agree.
                    values = [sample_evidence(r, caller)[caller] for r in records]
                    evidence[caller] = {field: values[0].get(field) if all(v.get(field) == values[0].get(field) for v in values) else None
                                        for field in ('depth', 'af')}
            if len(evidence) != 2:
                continue
            record = header.new_record(contig=chrom, start=pos - 1, alleles=(ref, alt))
            record.filter.add('PASS')
            for caller, tag in [('clairsto', 'CS'), ('deepsomatic', 'DS')]:
                for field, suffix in [('depth', '_DP'), ('af', '_AF')]:
                    value = evidence[caller][field]
                    if value is not None:
                        record.info[tag + suffix] = value
            db.execute('INSERT INTO records VALUES (?,?,?)', (ordinal, key, str(record)))
    sources.append(dict(caller='consensus', kind='derived_exact_pass', sample=sources[0]['sample'],
                        target_filter_counts=dict(pass_records=db.execute('SELECT COUNT(*) FROM records WHERE source=?', (ordinal,)).fetchone()[0])))
    (root / 'sources.json').write_text(json.dumps(sources))
    db.commit()
    db.close()
