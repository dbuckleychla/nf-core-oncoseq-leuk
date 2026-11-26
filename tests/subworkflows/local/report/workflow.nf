include { REVIEW_REPORT } from '../../../../subworkflows/local/report/review'

workflow CHECK_REVIEW {

    main:
    def root = file(params.review_fixture_root)
    def inputs = new groovy.json.JsonSlurperClassic().parseText(file("${root}/inputs.json").text)
    ch_inputs = channel.of(inputs)
        .flatMap { entries -> entries.collect { entry -> tuple('SYNTHETIC', entry.kind, file("${root}/${params.review_empty && entry.kind.startsWith('small:') ? entry.kind.split(':')[1] + '.empty.vcf.gz' : entry.path}")) } }
        .mix(channel.of(tuple('SYNTHETIC', 'sample', [])))
    // Version topics close after rendering; the attachment must not emit one.
    ch_versions = channel.topic('versions')
        .map { process, tool, version -> "${process}:\n  ${tool}: ${version}\n" }
        .toList().map { it + ['SYNTHETIC_TOOL:\n  synthetic: 1.20\n'] }
    REVIEW_REPORT(
        ch_inputs,
        channel.of(tuple([id: 'SYNTHETIC'], file("${root}/sample.bam"), file("${root}/sample.bam.bai"))),
        channel.value(tuple([id: 'reference'], file("${root}/genome.fa"))),
        channel.value(file("${root}/cache")),
        channel.value([workflow: [pipeline: 'SYNTHETIC_WORKFLOW'], parameters: [basecall_model: 'dna_synthetic@v1'],
            software_version_texts: []]),
        ch_versions,
    )

    emit:
    html = REVIEW_REPORT.out.html
    data = REVIEW_REPORT.out.data
}

workflow {
    CHECK_REVIEW()
}
