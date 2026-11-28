"""Read stage counts without changing preprocessing or exporting read records."""
import csv
import json

import pysam


def collect(descriptor):
    result = []
    for entry in json.loads(open(descriptor).read()):
        reads = bases = 0
        seen = set()
        for path in entry['files']:
            if path in seen:
                raise ValueError('Duplicate preprocessing file')
            seen.add(path)
            if entry['format'] == 'seqkit':
                with open(path) as source:
                    for row in csv.DictReader(source, delimiter='\t'):
                        reads += int(row['num_seqs'])
                        bases += int(row['sum_len'])
            else:
                with pysam.AlignmentFile(path, 'rb', check_sq=False) as source:
                    for read in source.fetch(until_eof=True):
                        if not read.is_secondary and not read.is_supplementary:
                            reads += 1
                            if read.query_length is None:
                                raise ValueError('Unknown length in preprocessing input')
                            bases += read.query_length
        result.append(dict(scope_type=entry['scope_type'], scope_id=entry['scope_id'], stage=entry['stage'],
                           source_kind=entry['format'], reads=reads, bases=bases, exhaustive=True))
    return result
