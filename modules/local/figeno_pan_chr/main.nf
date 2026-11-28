// Adapted from the release-pinned figeno module; bounded chromosome concurrency.
process FIGENO_PAN_CHR {
    // TODO : SET FIXED VERSION WHEN PIPELINE IS STABLE
    container 'ghcr.io/chusj-pigu/figeno:e42833528621000fed543c5a6d999262ab751bcf'

    tag "$meta.id"
    label 'local'

    input:
    tuple val(meta),
        path(delly_cov),
        path(delly_bed)

    output:
    tuple val(meta),
        path("*.png"),
        emit: figure
    path "versions.yml",
        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    //def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = Math.min(task.cpus as int, (task.ext.threads ?: 4) as int)
    """
    pan_chr.sh \\
        ${prefix} \\
        ${delly_cov} \\
        ${delly_bed} \\
        --jobs ${threads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        figeno: \$(awk '{print \$NF}' /opt/version.txt)
    END_VERSIONS
    """
}
