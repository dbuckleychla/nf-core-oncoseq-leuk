def alignmentQCShellQuote(value) {
    "'" + value.toString().replace("'", "'\\''") + "'"
}

process SAMTOOLS_ALIGNMENT_QC {
    tag "${meta.id}"
    label 'process_single_cpu'
    label 'process_low_memory'
    label 'process_low_time'
    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/htslib_samtools@sha256:a55ddea590e567a91df592300a960aa534cfc1bd16e7623e3938ec21f4f3df15'
    input:
    tuple val(meta), path(bam)
    output:
    tuple val(meta), path('alignment_qc.json'), emit: qc
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | head -1 | cut -d' ' -f2"), topic: versions, emit: versions_samtools
    script:
    """
    samtools flagstat -@ ${Math.max(0, task.cpus - 1)} -O json ${alignmentQCShellQuote(bam)} > alignment_qc.json
    """
}
