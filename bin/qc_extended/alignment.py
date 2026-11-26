"""Streamed alignment/coverage metrics. FrankONTstein interval/depth ideas,
revision 4655311, MIT (scripts/LICENSE.FrankONTstein); primary reads scanned once.
"""
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import math
import re
import subprocess
import sys

import pysam

from .regions import Regions, indexes, merged, primary_contig, read_bed

THRESHOLDS = (1, 5, 10, 20, 30, 50, 100)
GROUPS = ('all', 'mapped', 'unmapped', 'on_target', 'off_target', 'on_enrichment', 'off_enrichment')


def length_summary(counts, missing=0):
    n = sum(counts.values())
    bases = sum(k * v for k, v in counts.items())
    result = dict(records=n + missing, known_length_records=n, missing_length_records=missing,
                  bases=bases, mean=None, median=None, n50=None, minimum=None, maximum=None,
                  histogram={}, histogram_bases={})
    if not n:
        return result
    ranks, middle, acc = ((n - 1) // 2, n // 2), [], 0
    for length, count in sorted(counts.items()):
        middle.extend(length for rank in ranks if acc <= rank < acc + count)
        acc += count
    acc, n50 = 0, None
    if bases:
        for length, count in sorted(counts.items(), reverse=True):
            acc += length * count
            if acc * 2 >= bases:
                n50 = length
                break
    hist, histbases = Counter(), Counter()
    for length, count in counts.items():
        hist[length // 100 * 100] += count
        histbases[length // 100 * 100] += length * count
    result.update(mean=bases / n, median=sum(middle) / 2, n50=n50,
                  minimum=min(counts), maximum=max(counts), histogram=dict(hist), histogram_bases=dict(histbases))
    result['fractions'] = {}
    for label, predicate in [('le_500', lambda x: x <= 500), ('le_1000', lambda x: x <= 1000), ('ge_10000', lambda x: x >= 10000)]:
        result['fractions'][label] = dict(reads=sum(v for k, v in counts.items() if predicate(k)),
                                         bases=sum(k * v for k, v in counts.items() if predicate(k)))
    return result


def depth_summary(hist, size):
    hist = Counter(hist)
    hist[0] += size - sum(hist.values())
    if hist[0] < 0:
        raise ValueError('Coverage denominator smaller than observed positions')
    aligned = sum(k * v for k, v in hist.items())
    mean = aligned / size if size else None
    median, acc = None, 0
    if size:
        ranks, values = ((size - 1) // 2, size // 2), []
        for d, count in sorted(hist.items()):
            values.extend(d for rank in ranks if acc <= rank < acc + count)
            acc += count
        median = sum(values) / 2
    return dict(bases=size, aligned_bases=aligned, mean=mean, median=median,
                cv=math.sqrt(max(0, sum(k * k * v for k, v in hist.items()) / size - mean * mean)) / mean if mean else None,
                breadth={str(t): sum(v for k, v in hist.items() if k >= t) / size if size else None for t in THRESHOLDS},
                fraction_ge_0_2_mean=sum(v for k, v in hist.items() if k >= .2 * mean) / size if mean else None,
                histogram=dict(sorted(hist.items())))


def read_groups(read, enrichment, targets):
    if read.is_secondary or read.is_supplementary:
        return []
    if read.is_unmapped:
        return ['all', 'unmapped']
    blocks = read.get_blocks()
    target = any(next(targets.hits(a, b), None) is not None for a, b in blocks)
    enriched = any(next(enrichment.hits(a, b), None) is not None for a, b in blocks)
    return ['all', 'mapped', 'on_target' if target else 'off_target', 'on_enrichment' if enriched else 'off_enrichment']


def modifications(read):
    """Check legal encoding separately from usable mapped modification evidence."""
    tags = dict(read.get_tags())
    names = [x for x in ('MM', 'Mm') if x in tags]
    probabilities = [x for x in ('ML', 'Ml') if x in tags]
    result = dict(mm=int(bool(names)), ml=int(bool(probabilities)), mn=int('MN' in tags),
                  usable=0, invalid=0, m=0, h=0, unsupported=0, missing_probabilities=0,
                  implicit=0, explicit=0, empty=0)
    if not names:
        result['invalid'] = int(bool(probabilities))
        return result
    if len(names) > 1 or len(probabilities) > 1:
        result['invalid'] = 1
        return result
    value = tags[names[0]]
    if not isinstance(value, str):
        result['invalid'] = 1
        return result
    if not value:
        result['empty'] = 1
    try:
        total = 0
        for group in value.split(';'):
            if not group:
                continue
            match = re.fullmatch(r'([ACGTUN])([+-])([A-Za-z]+|\d+)([.?]?)(?:,([0-9]+(?:,[0-9]+)*))?', group)
            if not match:
                raise ValueError('Invalid MM syntax')
            base, strand, codes, skip, deltas = match.groups()
            result['explicit' if skip == '?' else 'implicit'] += 1
            if base == 'C':
                result['m'] |= int('m' in codes)
                result['h'] |= int('h' in codes)
                result['unsupported'] |= int(any(c not in 'mh' for c in codes))
            count = len(deltas.split(',')) if deltas else 0
            total += count * (1 if codes.isdigit() else len(codes))
        if 'MN' in tags and read.query_length is not None and tags['MN'] != read.query_length:
            raise ValueError('MN mismatch')
        if probabilities:
            array = tags[probabilities[0]]
            if len(array) != total or any(not 0 <= x <= 255 for x in array):
                raise ValueError('ML mismatch')
        else:
            result['missing_probabilities'] = 1
        # HTSlib validates canonical-base positions, orientation, and query length.
        if total:
            decoded = read.modified_bases
            if decoded is None:
                raise ValueError('Modification positions invalid')
        result['usable'] = int(bool(probabilities) and bool(result['m'] or result['h']) and not read.is_unmapped)
    except (ValueError, TypeError, IndexError, OverflowError):
        result['invalid'] = 1
        result['usable'] = 0
    return result


def scan_reads(bam, erows, trows):
    es, ts = indexes(erows), indexes(trows)
    empty = Regions([])
    counts, mods, mapq = Counter(), Counter(), Counter()
    lengths, missing, sources = defaultdict(Counter), Counter(), Counter()
    soft = 0
    quality_sum, quality_count = 0., 0
    error_sum, error_bases = 0., 0
    with pysam.AlignmentFile(bam, 'rb') as source:
        models = set()
        for rg in source.header.to_dict().get('RG', []):
            for model in re.findall(r'(?:basecall_model|modbase_models)=([^\s]+)', rg.get('DS', '')):
                if re.fullmatch(r'[a-zA-Z0-9_@.,+-]+', model):
                    models.add(model)
        for read in source.fetch(until_eof=True):
            counts['alignment_records'] += 1
            if read.is_secondary or read.is_supplementary:
                counts['secondary_records'] += int(read.is_secondary)
                counts['supplementary_records'] += int(read.is_supplementary)
                counts['both_nonprimary_flags'] += int(read.is_secondary and read.is_supplementary)
                continue
            counts['primary_records'] += 1
            counts['unmapped_primary' if read.is_unmapped else 'mapped_primary'] += 1
            counts['duplicate_primary'] += int(read.is_duplicate)
            counts['qcfail_primary'] += int(read.is_qcfail)
            groups = read_groups(read, es.get(read.reference_name, empty), ts.get(read.reference_name, empty))
            length = read.query_length
            origin = 'sequence'
            if length is None:
                length = read.infer_query_length()
                origin = 'cigar_inferred' if length is not None else 'unavailable'
            sources[origin] += 1
            for group in groups:
                if length is None:
                    missing[group] += 1
                else:
                    lengths[group][length] += 1
            counts['hard_clipped_primary'] += int(any(op == 5 for op, _ in read.cigartuples or []))
            if not read.is_unmapped:
                mapq[read.mapping_quality] += 1
                soft += sum(n for op, n in read.cigartuples or [] if op == 4)
            if read.has_tag('qs'):
                q = read.get_tag('qs')
                if isinstance(q, (float, int)) and math.isfinite(q) and q >= 0:
                    quality_sum += q
                    quality_count += 1
            if read.query_qualities is not None and len(read.query_qualities):
                if all(q != 255 for q in read.query_qualities):
                    error_sum += sum(10 ** (-q / 10) for q in read.query_qualities)
                    error_bases += len(read.query_qualities)
            mod = modifications(read)
            # Mapped evidence drives pileup; all-primary presence is still reported.
            for key, value in mod.items():
                mods[key] += value
                if not read.is_unmapped:
                    mods['mapped_' + key] += value
    for key in ('alignment_records', 'primary_records', 'mapped_primary', 'unmapped_primary', 'secondary_records',
                'supplementary_records', 'both_nonprimary_flags', 'duplicate_primary', 'qcfail_primary', 'hard_clipped_primary'):
        counts[key] += 0
    for key in ('mm', 'ml', 'mn', 'usable', 'invalid', 'm', 'h', 'unsupported', 'missing_probabilities', 'implicit', 'explicit', 'empty'):
        mods[key] += 0
        mods['mapped_' + key] += 0
    return dict(counts=dict(counts), modifications=dict(mods), models=sorted(models),
                read_lengths={g: length_summary(lengths[g], missing[g]) for g in GROUPS},
                length_sources=dict(sources), mapq=dict(mapq), soft_clipped_bases=soft,
                quality=dict(qs_records=quality_count, qs_mean=quality_sum / quality_count if quality_count else None,
                             error_probability=error_sum / error_bases if error_bases else None,
                             error_phred=-10 * math.log10(error_sum / error_bases) if error_sum and error_bases else None,
                             error_bases=error_bases))


def depth_contig(job):
    bam, chrom, size, erows, trows = job
    e, t = Regions(merged(erows)), Regions(merged(trows))
    sizes = dict(genome_primary=size, enrichment=sum(b - a for a, b in e.rows), targets=sum(b - a for a, b in t.rows))
    sizes['off_enrichment'] = size - sizes['enrichment']
    sizes['padding'] = sum(b - a - sum(min(b, y) - max(a, x) for x, y in t.hits(a, b)) for a, b in e.rows)
    results = {}
    for threshold in (0, 20):
        hist = {k: Counter() for k in sizes}
        thist = {r[2]: Counter() for r in trows}
        events = defaultdict(list)
        for key, rows in [('enrichment', e.rows), ('targets', t.rows)]:
            for a, b in rows:
                events[a].append((key, 1))
                events[b].append((key, -1))
        for a, b, key, _label in trows:
            events[a].append((key, 1))
            events[b].append((key, -1))
        boundaries, cursor, active = sorted(events), 0, set()
        # Stream the Samtools implementation shipped with the active pysam.
        # save_stdout avoids materializing a chromosome's depth as a Python string.
        depth_driver = "import pysam,sys; pysam.samtools.depth(*sys.argv[1:], catch_stdout=False, save_stdout='/dev/stdout')"
        cmd = [sys.executable, '-c', depth_driver, '-r', chrom, '-q', '0', '-Q', str(threshold), '-g', '0x600', '-G', '0x900', str(bam)]
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True) as proc:
            for line in proc.stdout:
                _chrom, pos, depth = line.split()
                pos, depth = int(pos) - 1, int(depth)
                while cursor < len(boundaries) and boundaries[cursor] <= pos:
                    for key, delta in events[boundaries[cursor]]:
                        if delta == 1:
                            active.add(key)
                        else:
                            active.discard(key)
                    cursor += 1
                hist['genome_primary'][depth] += 1
                hist['enrichment' if 'enrichment' in active else 'off_enrichment'][depth] += 1
                if 'targets' in active:
                    hist['targets'][depth] += 1
                elif 'enrichment' in active:
                    hist['padding'][depth] += 1
                for key in active - {'enrichment', 'targets'}:
                    thist[key][depth] += 1
            if proc.wait():
                raise RuntimeError('Samtools depth failed')
        results[str(threshold)] = dict(scopes={k: depth_summary(hist[k], size) for k, size in sizes.items()},
                                      targets={r[2]: depth_summary(thist[r[2]], r[1] - r[0]) for r in trows})
    return dict(contig=chrom, coverage=results)


def collect(bam, reference, enrichment, targets, threads=8):
    lengths = reference['lengths']
    erows, trows = read_bed(enrichment, lengths), read_bed(targets, lengths)
    with pysam.AlignmentFile(bam, 'rb') as source:
        if dict(zip(source.references, source.lengths)) != lengths:
            raise ValueError('BAM/FASTA sequence dictionary mismatch')
        if not source.has_index() or source.header.to_dict().get('HD', {}).get('SO') != 'coordinate':
            raise ValueError('QC requires coordinate-sorted indexed BAM')
    reads = scan_reads(bam, erows, trows)
    jobs = [(bam, c, size, erows.get(c, []), trows.get(c, [])) for c, size in lengths.items()]
    with ProcessPoolExecutor(max_workers=min(4, max(1, threads // 2))) as pool:
        contigs = list(pool.map(depth_contig, jobs))
    aggregate = {}
    for threshold in ('0', '20'):
        aggregate[threshold] = {}
        for scope in ('genome_primary', 'enrichment', 'targets', 'off_enrichment', 'padding'):
            hist, size = Counter(), 0
            for contig in contigs:
                if primary_contig(contig['contig']):
                    data = contig['coverage'][threshold]['scopes'][scope]
                    hist.update(data['histogram'])
                    size += data['bases']
            aggregate[threshold][scope] = depth_summary(hist, size)
    outside = sum(b - a - sum(min(b, y) - max(a, x) for x, y in indexes(erows).get(c, Regions([])).hits(a, b))
                  for c, rows in trows.items() for a, b in merged(rows))
    return dict(**reads, coverage=aggregate, contigs=contigs,
                target_intervals=[dict(target_id=key, label=label, contig=c, start=a, end=b)
                                  for c, rows in trows.items() for a, b, key, label in rows],
                targets_outside_enrichment_bases=outside,
                target_union_bases=sum(b - a for rows in trows.values() for a, b in merged(rows)))
