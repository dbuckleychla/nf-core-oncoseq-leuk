include { MINIMAP2_ALIGN } from '../../../../../modules/local/clinical_threads/minimap2_align/main'

process CHECK_UBAM_TAGS {
    cpus 1
    memory 1.GB
    container 'community.wave.seqera.io/library/oncoseq_review_1_0_0@sha256:b143e350e1dae33ad26758fc6f9ccf2392cb8b49881de73960f42525602d5eb5'
    input:
    tuple val(meta), path(bam), path(index)
    output:
    path 'tags.checked.txt'
    script:
    """
    python3 - <<'PY' > tags.checked.txt
    import pysam
    with pysam.AlignmentFile('${bam}') as source:
        reads=list(source)
    assert len(reads)==1004, 'SUP quality filtering changed'
    assert all(not read.is_unmapped for read in reads), 'Alignment failed'
    assert all(read.has_tag('MM') and read.has_tag('ML') and read.has_tag('MN') and read.has_tag('qs') and read.has_tag('st') for read in reads), 'Tags were lost'
    assert all(read.get_tag('MM')=='C+m,0;' and list(read.get_tag('ML'))==[220] and read.get_tag('MN')==1000 for read in reads), 'Modification evidence changed'
    assert all(read.query_length==1000 for read in reads), 'Read lengths changed'
    print('Aligned 1004 synthetic reads with methylation, timing and quality tags preserved')
    PY
    """
}

workflow UBAM_SMOKE {
    take:
    reads
    reference
    main:
    MINIMAP2_ALIGN(reads, reference, true, 'bai', false, false)
    CHECK_UBAM_TAGS(MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.index))
    emit:
    validation = MINIMAP2_ALIGN.out.ubam_validation
    checked = CHECK_UBAM_TAGS.out
}
