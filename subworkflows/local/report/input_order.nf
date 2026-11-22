def orderedReviewArtifacts(kinds, paths) {
    if (kinds.size() != paths.size()) {
        throw new IllegalArgumentException('Review artifact kinds and paths must have matching lengths')
    }
    def indexes = (0..<kinds.size()).findAll { index -> kinds[index] != 'sample' }
    indexes.sort { left, right ->
        (kinds[left] <=> kinds[right]) ?: (paths[left].toUriString() <=> paths[right].toUriString())
    }
    return [indexes.collect { index -> kinds[index] }, indexes.collect { index -> paths[index] }]
}
