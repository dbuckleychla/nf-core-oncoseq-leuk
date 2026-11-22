"""Bounded adaptive prefetch for exact-allele annotation, without filtering rows."""
from collections import OrderedDict, defaultdict

import pysam


class IndexedReference:
    def __init__(self, path, window=512, budget_bytes=16 * 1024 * 1024, max_windows=4):
        if window <= 0 or budget_bytes <= 0 or max_windows <= 0:
            raise ValueError('Reference cache limits must be positive')
        self.reader = pysam.VariantFile(str(path))
        self.window = window
        self.budget_bytes = budget_bytes
        self.max_windows = max_windows
        self.cache = OrderedDict()
        self.seen = OrderedDict()
        self.blocked = OrderedDict()
        self.cached_bytes = 0
        self.peak_cached_bytes = 0
        self.direct_fetches = 0
        self.window_fetches = 0
        self.cache_hits = 0
        self.oversized_windows = 0

    @property
    def header(self):
        return self.reader.header

    def _direct(self, query):
        self.direct_fetches += 1
        return [r for r in self.reader.fetch(query.contig, query.start, query.start + 1)
                if r.pos == query.pos and r.ref == query.ref and r.alts == query.alts]

    def _evict(self):
        _, (_, charge) = self.cache.popitem(last=False)
        self.cached_bytes -= charge

    def exact(self, query):
        if query.contig not in self.header.contigs:
            return []
        start = query.start // self.window * self.window
        location = (query.contig, start)
        allele = (query.pos, query.ref, query.alts)
        if location in self.cache:
            self.cache_hits += 1
            self.cache.move_to_end(location)
            return self.cache[location][0].get(allele, [])
        if location in self.blocked:
            return self._direct(query)
        if location not in self.seen:
            self.seen[location] = True
            if len(self.seen) > 32:
                self.seen.popitem(last=False)
            return self._direct(query)
        self.seen.move_to_end(location)
        grouped, charge = defaultdict(list), 0
        self.window_fetches += 1
        for record in self.reader.fetch(query.contig, start, start + self.window):
            # Conservatively account for native storage, keys, lists and objects.
            charge += 4 * len(str(record).encode()) + 512
            while self.cache and self.cached_bytes + charge > self.budget_bytes:
                self._evict()
            if charge > self.budget_bytes:
                self.oversized_windows += 1
                self.blocked[location] = True
                if len(self.blocked) > 32:
                    self.blocked.popitem(last=False)
                # Never use an incomplete window: exact query remains authoritative.
                return self._direct(query)
            grouped[(record.pos, record.ref, record.alts)].append(record)
        while len(self.cache) >= self.max_windows:
            self._evict()
        self.cache[location] = (grouped, charge)
        self.cached_bytes += charge
        self.peak_cached_bytes = max(self.peak_cached_bytes, self.cached_bytes)
        return grouped.get(allele, [])

    def statistics(self):
        return dict(direct_fetches=self.direct_fetches, window_fetches=self.window_fetches,
                    cache_hits=self.cache_hits, oversized_windows=self.oversized_windows,
                    peak_cache_accounted_bytes=self.peak_cached_bytes,
                    cache_budget_bytes=self.budget_bytes)

    def close(self):
        self.cache.clear()
        self.seen.clear()
        self.blocked.clear()
        self.cached_bytes = 0
        self.reader.close()
