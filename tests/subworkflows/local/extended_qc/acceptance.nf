include { EXTENDED_QC } from '../../../../subworkflows/local/extended_qc/main'

workflow {
    def root = file(params.qc_fixture_dir, checkIfExists: true)
    def reference = channel.value([[id: 'reference'], file("${root}/reference.fa")])
    def fai = channel.value([[id: 'reference'], file("${root}/reference.fa.fai")])
    def bams = channel.of([[id: 'SYNTHETIC_A'], file("${root}/plain.bam"), file("${root}/plain.bam.bai")],
                          [[id: 'SYNTHETIC_B'], file("${root}/plain.bam"), file("${root}/plain.bam.bai")])
    def regions = channel.of([[id: 'SYNTHETIC_A'], file("${root}/enrichment.bed"), 0, []],
                             [[id: 'SYNTHETIC_B'], file("${root}/enrichment.bed"), 0, []])
    EXTENDED_QC(bams, reference, fai, regions, channel.empty(), 'hg38')
}
