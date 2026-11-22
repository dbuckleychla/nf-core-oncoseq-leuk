#!/usr/bin/env python3
"""Local synthetic benchmarks of unchanged review code versus contained prototypes.

No clinical inputs or AWS calls. Candidate implementations live only here;
production review code and resource requests are not modified.
"""
import argparse
from collections import OrderedDict, defaultdict
import hashlib
import inspect
import json
import platform
from pathlib import Path
import random
import resource
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'bin/review'))
import pysam
import annotate_small_variants as small
import structural_consensus as structural


class WindowLookup:
    """Bounded, exact-key cache; preserve native record order and duplicates."""
    def __init__(self, handle, window, max_windows=4):
        self.handle = handle
        self.window = window
        self.max_windows = max_windows
        self.cache = OrderedDict()
        self.fetches = 0

    def exact(self, query):
        if query.contig not in self.handle.header.contigs:
            return []
        start = query.start // self.window * self.window
        location = (query.contig, start)
        if location not in self.cache:
            grouped = defaultdict(list)
            self.fetches += 1
            for record in self.handle.fetch(query.contig, start, start + self.window):
                grouped[(record.pos, record.ref, record.alts)].append(record)
            self.cache[location] = grouped
            if len(self.cache) > self.max_windows:
                self.cache.popitem(last=False)
        self.cache.move_to_end(location)
        return self.cache[location].get((query.pos, query.ref, query.alts), [])


def indexed_edges(rows, distance):
    """Same compatibility graph, precomputed junctions and bounded candidates."""
    groups = defaultdict(list)
    edges = {i: set() for i in range(len(rows))}
    for index, row in enumerate(rows):
        junction = structural.junction(row)
        if junction:
            kind, left, right = junction
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


def indexed_consolidate():
    # Reuse the exact existing event/ambiguity/selection code. Replace only the
    # all-pairs edge-building block in the trusted repository function source.
    source = inspect.getsource(structural.consolidate)
    start = source.index('    edges = ')
    end = source.index('    seen, events = ', start)
    source = source[:start] + '    edges = indexed_edges(rows, distance)\n' + source[end:]
    namespace = dict(vars(structural), indexed_edges=indexed_edges)
    exec(compile(source, '<synthetic-indexed-consensus-prototype>', 'exec'), namespace)
    return namespace['consolidate']


def make_vcf(root, records):
    header = pysam.VariantHeader()
    header.contigs.add('chr1', length=records * 10 + 2000)
    header.contigs.add('chr2', length=10000)
    header.info.add('AC', 'A', 'Integer', 'Synthetic allele count')
    header.info.add('AN', 1, 'Integer', 'Synthetic allele denominator')
    header.info.add('FILLER', 1, 'String', 'Synthetic annotation payload')
    path = root / 'synthetic-reference.vcf.gz'
    rng = random.Random(719)
    with pysam.VariantFile(str(path), 'wz', header=header) as writer:
        for i in range(records):
            # Rich rows approximate expensive BGZF parsing; data are invented.
            payload = ''.join(rng.choices('ACGT', k=2048))
            record = writer.new_record(contig='chr1', start=10 * i + 1,
                                       alleles=('A', 'C'), id='SYNTHETIC_' + str(i))
            record.info['AC'] = (1,)
            record.info['AN'] = 2000
            record.info['FILLER'] = payload
            writer.write(record)
            if i % 997 == 0:
                writer.write(record)  # Preserve duplicate catalogue observations.
                record.alts = ('C', 'G')
                record.info['AC'] = (1, 2)
                writer.write(record)  # Exact single-alt queries must not match.
    pysam.tabix_index(str(path), preset='vcf', force=True)
    return path


def lookup_worker(path, pattern, count, window):
    with pysam.VariantFile(str(path)) as handle:
        cache = WindowLookup(handle, window) if window else None
        digest = hashlib.sha256()
        found = 0
        began = time.perf_counter()
        for i in range(count):
            index = i if pattern == 'dense' else i * 100
            query = handle.new_record(contig='chr1', start=10 * index + 1, alleles=('A', 'C'))
            records = cache.exact(query) if cache else small.exact_records(handle, query)
            digest.update(str(query.pos).encode() + b'\0')
            for record in records:
                digest.update(str(record).encode())
                found += 1
        elapsed = time.perf_counter() - began
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return dict(seconds=elapsed, queries=count, matched_records=found,
                    indexed_fetches=cache.fetches if cache else count,
                    logical_digest=digest.hexdigest(),
                    peak_rss_mib=peak / (1024 * 1024 if sys.platform == 'darwin' else 1024))


def make_junctions(count, dense=False):
    rows = []
    for i in range(count):
        group, member = divmod(i, 3)
        caller = ('sniffles', 'severus', 'stellerator')[member]
        offset = (0, 80, 160)[member] if group % 11 == 0 else (0, 4, 8)[member]
        position = 1000 + (0 if dense else group * 1000) + offset
        rows.append(dict(id='SYNTHETIC_JUNCTION_' + str(i), sample='SYNTHETIC', caller=caller,
                         chrom='chr1' if dense else 'chr' + str(1 + group % 5), pos=position,
                         partner=dict(chrom='chr9', pos=position + 10000,
                                      local_side='right', remote_side='left' if i % 31 else None),
                         evidence=dict(info=dict(SVTYPE='BND', SUPPORT=10), format={}),
                         filter=['PASS'], genes1=[dict(name='GENE_A')], genes2=[dict(name='GENE_B')],
                         gene_pair=['GENE_A', 'GENE_B'], source_records=[],
                         coordinate_convention='Synthetic VCF 1-based POS'))
    return rows


def structural_worker(count, method, dense):
    rows = make_junctions(count, dense)
    function = structural.consolidate if method == 'baseline' else indexed_consolidate()
    began = time.perf_counter()
    events = function(rows, {('GENE_A', 'GENE_B')}, 100)
    elapsed = time.perf_counter() - began
    encoded = json.dumps(events, sort_keys=True, separators=(',', ':')).encode()
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return dict(seconds=elapsed, input_records=count, events=len(events),
                logical_digest=hashlib.sha256(encoded).hexdigest(),
                peak_rss_mib=peak / (1024 * 1024 if sys.platform == 'darwin' else 1024))


def invoke(*args):
    result = subprocess.run([sys.executable, __file__, '--worker', *map(str, args)],
                            check=True, text=True, capture_output=True)
    return json.loads(result.stdout)


def benchmark(output):
    report = dict(environment=dict(python=platform.python_version(), machine=platform.machine(),
                                   system=platform.system(), pysam=pysam.__version__),
                  synthetic=True, production_changed=False, lookup=[], junctions=[])
    with tempfile.TemporaryDirectory(prefix='oncoseq-synthetic-benchmark-') as tmp:
        path = make_vcf(Path(tmp), 30000)
        report['reference_compressed_mib'] = path.stat().st_size / 2**20
        for pattern, count in [('dense', 20000), ('sparse', 300)]:
            expected = None
            for window in [0, 128, 512, 1024]:
                # Each case gets its own process for comparable peak RSS.
                result = invoke('lookup', path, pattern, count, window)
                result.update(pattern=pattern, window_bp=window)
                if expected is None:
                    expected = result['logical_digest']
                    baseline = result['seconds']
                assert result['logical_digest'] == expected, 'Exact lookup records changed'
                result['speedup'] = baseline / result['seconds']
                report['lookup'].append(result)
                print('Lookup', pattern, window, round(result['seconds'], 3),
                      'seconds;', round(result['speedup'], 2), 'x', flush=True)
        for count, dense in [(600, False), (1800, False), (3600, False), (180, True)]:
            baseline = invoke('junction', count, 'baseline', int(dense))
            candidate = invoke('junction', count, 'indexed', int(dense))
            assert baseline['logical_digest'] == candidate['logical_digest'], 'Consensus decisions changed'
            row = dict(count=count, dense_collisions=dense, baseline=baseline, indexed=candidate,
                       speedup=baseline['seconds'] / candidate['seconds'])
            report['junctions'].append(row)
            print('Junctions', count, 'dense', dense, round(row['speedup'], 2), 'x', flush=True)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + '\n')
    return report


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--worker':
        kind, *args = sys.argv[2:]
        result = (lookup_worker(Path(args[0]), args[1], int(args[2]), int(args[3]))
                  if kind == 'lookup' else structural_worker(int(args[0]), args[1], bool(int(args[2]))))
        print(json.dumps(result))
    else:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('--output', required=True)
        benchmark(parser.parse_args().output)
