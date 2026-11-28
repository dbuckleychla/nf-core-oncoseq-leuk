#!/usr/bin/env python3
"""Generate explicitly synthetic QC material; never copy clinical output."""
import argparse
from array import array
import gzip
import json
from pathlib import Path
import pysam


def build(root):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    (root / 'reference.fa').write_text('>chr1\n' + 'CG' * 50 + '\n')
    pysam.faidx(str(root / 'reference.fa'))
    (root / 'enrichment.bed').write_text('chr1\t0\t10\tone\nchr1\t5\t15\ttwo\n')
    (root / 'targets.bed').write_text('chr1\t5\t10\tSYNTHETIC_TARGET\n')
    header = dict(HD=dict(SO='coordinate'), SQ=[dict(SN='chr1', LN=100)])
    for variant in ('plain', 'tagged', 'invalid'):
        with pysam.AlignmentFile(str(root / (variant + '.bam')), 'wb', header=header) as out:
            for name, start, length, mapq, flag, cigar in [('A', 0, 10, 60, 0, '10M'), ('B', 10, 10, 0, 0, '10M'),
                                                        ('SUPP', 10, 10, 60, 2048, '10M'), ('C', 20, 20, 60, 0, '10M5D10M'),
                                                        ('UNMAPPED', -1, 5, 0, 4, None)]:
                r = pysam.AlignedSegment()
                r.query_name = 'SYNTHETIC_' + name
                r.query_sequence = 'C' * length
                r.query_qualities = array('B', [30] * length)
                r.flag = flag
                if start >= 0:
                    r.reference_id, r.reference_start, r.mapping_quality, r.cigarstring = 0, start, mapq, cigar
                if variant != 'plain':
                    r.set_tag('MM', 'C+m,0;')
                    r.set_tag('ML', array('B', [240, 240] if variant == 'invalid' else [240]))
                    r.set_tag('MN', length)
                out.write(r)
        pysam.index(str(root / (variant + '.bam')))
    with gzip.open(root / 'mods.bed.gz', 'wt') as out:
        for code, modified, other in [('m', 3, 2), ('h', 2, 3)]:
            out.write(f'chr1\t0\t1\t{code}\t10\t.\t0\t1\t0\t10\t{modified*10}\t{modified}\t5\t{other}\t0\t0\t0\t0\n')
    inputs = dict(qc_fixture_dir=str(root.resolve()), outdir=str((root / 'results').resolve()))
    (root / 'inputs.json').write_text(json.dumps(inputs) + '\n')
    return root


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', required=True)
    build(parser.parse_args().outdir)
