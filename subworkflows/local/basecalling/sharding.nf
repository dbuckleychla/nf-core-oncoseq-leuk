def planPod5Shards(rows, requestedShards) {
    if (!(requestedShards?.toString() ==~ /[1-9][0-9]*/)) {
        throw new IllegalArgumentException('basecall_shards must be a positive integer')
    }
    if (!rows) {
        throw new IllegalArgumentException('No POD5 inputs supplied')
    }

    def count = requestedShards.toString().toInteger()
    def seen = [] as Set
    def groups = rows.groupBy { row -> row[0].id }
    def result = []

    groups.each { sampleId, sampleRows ->
        def meta = sampleRows[0][0]
        def assets = [2, 3, 4].collect { column ->
            def values = sampleRows.collect { row -> row[column] }.unique()
            if (values.size() != 1) {
                throw new IllegalArgumentException("Inconsistent resume/model inputs for '${sampleId}'")
            }
            values[0]
        }
        if (sampleRows.any { row -> row[0] != meta }) {
            throw new IllegalArgumentException("Inconsistent metadata for '${sampleId}'")
        }

        def resolved = sampleRows
            .collectMany { row ->
                def inputs = row[1] instanceof Collection ? row[1] : [row[1]]
                inputs.collectMany { value ->
                    def input = file(value, checkIfExists: true)
                    def paths = input.name.endsWith('.pod5') && !input.isDirectory()
                        ? [input]
                        : files("${input.toUriString().replaceAll('/+$', '')}/*.pod5", checkIfExists: true)
                    paths.collect { path -> path.toUri().scheme == 'file' ? path.toRealPath() : path.normalize() }
                }
            }
            .sort { left, right -> left.toUriString() <=> right.toUriString() }

        if (!resolved) {
            throw new IllegalArgumentException("No POD5 files for '${sampleId}'")
        }
        resolved.each { path ->
            if (!seen.add(path.toUriString())) {
                throw new IllegalArgumentException("Duplicate POD5 input: ${path}")
            }
        }
        def effective = Math.min(count, resolved.size())
        if (effective > 1 && assets[0].name != 'NOFILE') {
            throw new IllegalArgumentException('A shared resume uBAM cannot be used across shards; use basecall_shards=1 or Nextflow -resume')
        }
        def batches = (0..<effective).collect { [] }
        resolved.eachWithIndex { path, index -> batches[index % effective] << path }
        batches.eachWithIndex { batch, index ->
            def shardMeta = meta + [id: "${sampleId}__bc${String.format('%04d', index + 1)}", basecall_parent: new LinkedHashMap(meta), basecall_shard: index + 1, basecall_count: effective]
            result << tuple(shardMeta, batch, assets[0], assets[1], assets[2])
        }
    }
    return result
}

def orderedShardBams(key, indexes, bams) {
    def expected = key.getGroupSize()
    if (bams.size() != expected || indexes.sort(false) != (1..expected).toList()) {
        throw new IllegalArgumentException('Incomplete or duplicate basecalling shards')
    }
    def order = (0..<indexes.size()).toList().sort { left, right -> indexes[left] <=> indexes[right] }
    tuple(key.getGroupTarget(), order.collect { index -> bams[index] })
}
