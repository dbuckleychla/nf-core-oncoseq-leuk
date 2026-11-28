"""Half-open interval and reference helpers.

Interval indexing adapted from FrankONTstein 4655311; MIT, see
scripts/LICENSE.FrankONTstein. No source report/header export is reused.
"""
from bisect import bisect_left
from collections import defaultdict
import gzip
import re


def text_open(path):
    with open(path, 'rb') as source:
        compressed = source.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if compressed else open(path)


def primary_contig(name):
    if re.fullmatch(r'(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)', name):
        return True
    match = re.fullmatch(r'NC_0*(\d+)\.\d+', name)
    return bool(match and (1 <= int(match[1]) <= 24 or 60925 <= int(match[1]) <= 60948))


def merged(rows):
    result = []
    for start, end, *_ in sorted(rows):
        if result and start <= result[-1][1]:
            result[-1][1] = max(result[-1][1], end)
        else:
            result.append([start, end])
    return result


class Regions:
    def __init__(self, rows):
        self.rows = sorted(rows)
        self.starts = [r[0] for r in self.rows]
        self.maxends = []
        maximum = 0
        for row in self.rows:
            maximum = max(maximum, row[1])
            self.maxends.append(maximum)

    def hits(self, start, end):
        i = bisect_left(self.starts, end) - 1
        while i >= 0 and self.maxends[i] > start:
            if self.rows[i][1] > start:
                yield self.rows[i]
            i -= 1

    def contains(self, pos):
        return next(self.hits(pos, pos + 1), None) is not None


def read_bed(path, lengths, features=False):
    result = defaultdict(list)
    with text_open(path) as source:
        for number, line in enumerate(source, 1):
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 3:
                raise ValueError('BED requires at least three columns')
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if chrom not in lengths or not 0 <= start < end <= lengths[chrom]:
                raise ValueError('BED coordinates do not match reference')
            if features and end != start + 1:
                raise ValueError('CpG features require one-base C coordinates')
            label = fields[3] if len(fields) > 3 else f'{chrom}:{start}-{end}'
            if re.search(r'[\x00-\x1f]', label):
                raise ValueError('Control character in BED label')
            result[chrom].append((start, end, f'{number}:{label}', label))
    if not result:
        raise ValueError('Empty BED union')
    return dict(result)


def indexes(rows):
    return {c: Regions(merged(r)) for c, r in rows.items()}


def reference(fasta, output, cpgs):
    import pysam
    from .contracts import checksum, write_json
    with pysam.FastaFile(fasta) as ref, gzip.open(cpgs, 'wt') as out:
        lengths = dict(zip(ref.references, ref.lengths))
        counts = {}
        for chrom, size in lengths.items():
            count = 0
            for start in range(0, size, 1_000_000):
                seq = ref.fetch(chrom, start, min(size, start + 1_000_001)).upper()
                for match in re.finditer('CG', seq):
                    pos = start + match.start()
                    out.write(f'{chrom}\t{pos}\n')
                    count += 1
            counts[chrom] = count
    write_json(output, dict(lengths=lengths, cpg_counts=counts, reference_sha256=checksum(fasta)))
