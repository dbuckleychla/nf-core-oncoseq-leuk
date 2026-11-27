include { CLASSY_BLOOD_WORKFLOW } from '../../../../subworkflows/local/methylation_analysis/classy_blood'

workflow CHECK_BLOOD {
    main:
    def root = "${projectDir}/test_data/review"
    CLASSY_BLOOD_WORKFLOW(
        channel.of(tuple([id: 'SYNTHETIC'], file("${root}/sample.bam"), file("${root}/sample.bam.bai"))),
        channel.value(tuple([id: 'hg38'], file("${root}/genome.fa"))),
    )

    emit:
    json = CLASSY_BLOOD_WORKFLOW.out.json
    section = CLASSY_BLOOD_WORKFLOW.out.sections
    versions = CLASSY_BLOOD_WORKFLOW.out.versions
}

workflow {
    CHECK_BLOOD()
}
