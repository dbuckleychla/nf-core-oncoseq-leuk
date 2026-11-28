def reviewShellQuote(value) {
    "'" + value.toString().replace("'", "'\\''") + "'"
}

process REVIEW_PREPARE {
    tag "${meta.id}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/oncoseq_review_1_0_0@sha256:b143e350e1dae33ad26758fc6f9ccf2392cb8b49881de73960f42525602d5eb5'

    input:
    tuple val(meta), val(kinds), path(artifacts, stageAs: 'inputs/*')
    path reference_files, stageAs: 'references/*'
    val reference_manifest
    path panel
    val run_provenance
    val _code_digest

    output:
    tuple val(meta), path('catalogue'), path('payload'), emit: prepared
    tuple val("${task.process}"), val('oncoseq-review'), eval('review_report.py --version'), topic: versions, emit: versions

    script:
    def files = artifacts instanceof Collection ? artifacts : [artifacts]
    def inputs = groovy.json.JsonOutput.toJson(
        (0..<files.size()).collect { i ->
            [kind: kinds[i], path: files[i].toString(), name: files[i].name]
        }
    )
    def refs = groovy.json.JsonOutput.toJson(reference_manifest)
    def sample = reviewShellQuote(meta.id)
    def provenance = groovy.json.JsonOutput.toJson(run_provenance)
    """
    cat > inputs.json <<'REVIEW_INPUTS'
    ${inputs}
    REVIEW_INPUTS
    cat > references.json <<'REVIEW_REFERENCES'
    ${refs}
    REVIEW_REFERENCES
    cat > run.provenance.json <<'REVIEW_PROVENANCE'
    ${provenance}
    REVIEW_PROVENANCE
    review_report.py prepare --sample ${sample} --inputs inputs.json --references references.json \\
        --panel ${reviewShellQuote(panel)} --provenance run.provenance.json --output catalogue --payload payload
    """
}

process REVIEW_VEP {
    tag "${meta.id}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/ensembl-vep_perl-math-cdf@sha256:442c7140afaaa36f17671eb42656d6412b5196bbe742751b88d7eae5c7cc4d4b'

    input:
    tuple val(meta), path(catalogue), path(payload)
    tuple val(ref_meta), path(fasta)
    path cache
    val _code_digest
    val fasta_sha256

    output:
    tuple val(meta), path(catalogue), path(payload), path('union.vep.vcf.gz'), emit: annotated
    tuple val("${task.process}"), val('ensemblvep'), eval("${reviewShellQuote(task.ext.executable ?: 'vep')} --help | sed -n '/ensembl-vep/s/.*: //p'"), topic: versions, emit: versions

    script:
    def threads = task.ext.threads ?: task.cpus
    def executable = reviewShellQuote(task.ext.executable ?: 'vep')
    def command = meta.review_candidates > 0
        ? """
        ${executable} --input_file ${reviewShellQuote("${catalogue}/catalog.vcf")} --output_file union.vep.vcf.gz \\
            --vcf --compress_output bgzip --everything --offline --cache --cache_version 115 \\
            --assembly GRCh38 --species homo_sapiens --dir_cache ${reviewShellQuote(cache)} \\
            --fasta ${reviewShellQuote(fasta)} --fork ${threads} --force_overwrite
    """
        : "gzip -c ${reviewShellQuote("${catalogue}/empty.vep.vcf")} > union.vep.vcf.gz"
    """
    if [ -n ${reviewShellQuote(fasta_sha256)} ]; then
        review_fasta_sha=\$(sha256sum < ${reviewShellQuote(fasta)})
        test "\${review_fasta_sha:0:64}" = ${reviewShellQuote(fasta_sha256)} || { echo 'Analytical FASTA checksum mismatch' >&2; exit 1; }
    fi
    ${command}
    """
}

process REVIEW_COMPLETE {
    tag "${meta.id}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/oncoseq_review_1_0_0@sha256:b143e350e1dae33ad26758fc6f9ccf2392cb8b49881de73960f42525602d5eb5'

    input:
    tuple val(meta), path(catalogue), path(payload), path(vep)
    path reference_files, stageAs: 'references/*'
    val reference_manifest
    path panel
    path known_fusions
    path image_lock
    val code_digest

    output:
    tuple val(meta), path('review_data'), emit: data
    tuple val("${task.process}"), val('oncoseq-review'), eval('review_report.py --version'), topic: versions, emit: versions

    script:
    def refs = groovy.json.JsonOutput.toJson(reference_manifest)
    def sample = reviewShellQuote(meta.id)
    """
    cat > references.json <<'REVIEW_REFERENCES'
    ${refs}
    REVIEW_REFERENCES
    review_report.py complete --sample ${sample} --catalogue ${reviewShellQuote(catalogue)} --payload ${reviewShellQuote(payload)} \\
        --vep ${reviewShellQuote(vep)} --references references.json --panel ${reviewShellQuote(panel)} \\
        --known-fusions ${reviewShellQuote(known_fusions)} --image-lock ${reviewShellQuote(image_lock)} --code-digest ${code_digest} --distance ${params.review_fusion_distance} --output review_data
    """
}

process REVIEW_RENDER {
    tag "${meta.id}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/oncoseq_review_1_0_0@sha256:b143e350e1dae33ad26758fc6f9ccf2392cb8b49881de73960f42525602d5eb5'

    input:
    tuple val(meta), path(data), path(figures)
    path template
    val _code_digest

    output:
    tuple val(meta), path('*.review.html'), emit: html
    tuple val("${task.process}"), val('oncoseq-review'), eval('review_report.py --version'), topic: versions, emit: versions

    script:
    def output = reviewShellQuote("${meta.id}.review.html")
    // Anchor generated figures to the task directory, not the review JSON parent.
    def additional = figures ? "--figures \"\$(pwd)/\"${reviewShellQuote(figures)}" : ''
    """
    review_report.py render --input ${reviewShellQuote("${data}/review.json")} --output ${output} --template ${reviewShellQuote(template)} ${additional}
    """
}

process REVIEW_PROVENANCE_ATTACH {
    tag "${meta.id}"
    label 'process_low'
    container 'community.wave.seqera.io/library/oncoseq_review_1_0_0@sha256:b143e350e1dae33ad26758fc6f9ccf2392cb8b49881de73960f42525602d5eb5'

    input:
    tuple val(meta), path(html, stageAs: 'native_input/*')
    val software_versions
    val _code_digest

    output:
    tuple val(meta), path('*.review.html'), emit: html

    // No versions topic: this step runs after version producers finish.
    script:
    def versions = groovy.json.JsonOutput.toJson(software_versions)
    """
    cat > software_versions.json <<'REVIEW_VERSIONS'
    ${versions}
    REVIEW_VERSIONS
    review_report.py attach-provenance --input ${reviewShellQuote(html)} \\
        --versions software_versions.json --output ${reviewShellQuote(meta.id + '.review.html')}
    """
}

process REVIEW_FIGURES {
    tag "${meta.id}"
    label 'process_medium'
    container 'ghcr.io/chusj-pigu/figeno@sha256:e9a90169d979043902f3584990cc05d989e229b60fb90a019e7c65724503e98c'

    input:
    tuple val(meta), path(data), path(bam), path(bai)
    val _code_digest

    output:
    tuple val(meta), path(data), path('figures'), emit: figures
    tuple val("${task.process}"), val('figeno'), eval("awk '{print \$NF}' /opt/version.txt"), topic: versions, emit: versions

    script:
    def threads = task.ext.threads ?: 4
    """
    review_figures.py --input ${reviewShellQuote("${data}/review.json")} --bam ${reviewShellQuote(bam)} --output figures --jobs ${threads}
    """
}

process REVIEW_VALIDATE {
    tag 'Review references'
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/oncoseq_review_1_0_0@sha256:b143e350e1dae33ad26758fc6f9ccf2392cb8b49881de73960f42525602d5eb5'

    input:
    path reference_files, stageAs: 'references/*'
    val reference_manifest
    val _code_digest
    tuple val(_reference_meta), path(analytical_fasta)

    output:
    path 'validated.json', emit: validated
    tuple val("${task.process}"), val('oncoseq-review'), eval('review_report.py --version'), topic: versions, emit: versions

    script:
    def refs = groovy.json.JsonOutput.toJson(reference_manifest)
    """
    cat > references.json <<'REVIEW_REFERENCES'
    ${refs}
    REVIEW_REFERENCES
    review_report.py validate --references references.json --analytical-fasta ${reviewShellQuote(analytical_fasta)} --output validated.json
    """
}
