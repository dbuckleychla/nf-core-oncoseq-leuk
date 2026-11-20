// Stage each deterministic file group as one Dorado input directory.
process DORADO_BASECALL {
    container 'ghcr.io/chusj-pigu/dorado:851b0a37ecbff6b3b952c811b31dad48ca6bf995'

    // nf-core resource label
    label "process_high"

    // MPGI DRAC resource labels
    label "process_medium_low_cpu"
    label params.wgs ? "process_highest_memory" : "process_higher_memory"
    label params.wgs ? "process_high_time" : "process_medium_high_time"
    label "process_gpu"

    tag "${meta.id}"

    input:
    tuple val(meta), path(pod5, stageAs: 'pod5_inputs/input????????.pod5'), path(ubam), path(model), path(model_mh)

    output:
    tuple val(meta), path("*.bam"), emit: ubam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def device = params.device != null ? "-x ${params.device}" : ""
    def mod = model_mh.name != "NOMOD"
        ? "--modified-bases-models ${model_mh}"
        : ""
    def multi = params.demux ? "--no-trim" : ""
    def resume = ubam.name != 'NOFILE'
        ? "--resume-from ${ubam} > ${prefix}_unaligned_final.bam"
        : "> ${prefix}_unaligned.bam"
    """
    dorado basecaller \\
        ${args} \\
        ${device} \\
        ${model} \\
        ${mod} \\
        pod5_inputs \\
        ${multi} \\
        ${resume}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | \
            sed 's/^.*dorado //; s/Using.*\$//')
    END_VERSIONS
    """
}
