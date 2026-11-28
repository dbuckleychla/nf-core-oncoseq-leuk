include { DORADO_BASECALL } from '../../../modules/local/dorado_sharded/main.nf'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_BASECALL } from '../../../modules/nf-core/samtools/merge/main'
include { planPod5Shards ; orderedShardBams } from './sharding'

workflow BASECALL_SHARDED {
    take:
    ch_samplesheet

    main:
    ch_shards = ch_samplesheet
        .collect(flat: false)
        .flatMap { rows -> planPod5Shards(rows, params.basecall_shards) }

    DORADO_BASECALL(ch_shards)

    ch_grouped = DORADO_BASECALL.out.ubam
        .map { meta, bam -> tuple(groupKey(meta.basecall_parent, meta.basecall_count), meta.basecall_shard, bam) }
        .groupTuple(remainder: true)
        .map { key, indexes, bams -> orderedShardBams(key, indexes, bams) }
        .branch { _meta, bams ->
            single: bams.size() == 1
            multiple: true
        }

    SAMTOOLS_MERGE_BASECALL(
        ch_grouped.multiple.map { meta, bams -> tuple(meta, bams, []) },
        channel.value([[id: 'unaligned'], [], [], []]),
        '',
    )

    ch_ubam = ch_grouped.single
        .map { meta, bams -> tuple(meta, bams[0]) }
        .mix(SAMTOOLS_MERGE_BASECALL.out.bam)

    emit:
    ubam = ch_ubam
    versions = DORADO_BASECALL.out.versions
}
