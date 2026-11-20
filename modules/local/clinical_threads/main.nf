// Preserve upstream process interfaces while separating tool threads from worker reservations.

process CLAIR3_CALL {

    container "ghcr.io/chusj-pigu/clair3:5c1b224addf317c2cfe5ab5e03f6cdb23da13303"

    label "process_gpu_clair3"
    label 'process_high'
    // nf-core labels
    label "process_high_cpu"
    // Label for mpgi drac cpu alloc
    label "process_higher_memory"
    // Label for mpgi drac memory alloc
    label "process_medium_low_time"
    label "singleton"

    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai), path(bed)
    tuple val(meta2), path(ref)
    tuple val(meta3), path(fai)
    val model

    output:
    tuple val(meta), path("${meta.id}/merge_output.vcf.gz"), emit: vcf
    tuple val(meta), path("${meta.id}/merge_output.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp -L ${ref} ref_local.fa
    cp -L ${fai} ref_local.fa.fai

    run_clair3.sh \\
        ${args} \\
        --threads=${task.ext.threads ?: task.cpus} \\
        --sample_name=${prefix} \\
        --bed_fn=${bed} \\
        --platform='ont' \\
        --model_path="/opt/models/${model}" \\
        --bam_fn=\${PWD}/${bam} \\
        --ref_fn=\${PWD}/ref_local.fa \\
        -o \${PWD}/${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Clair3:  \$(echo \$(run_clair3.sh --version 2>&1) | awk '{print \$NF}' )
    END_VERSIONS
    """
}

process CLAIRS_TO_CALL {

    container "ghcr.io/chusj-pigu/clairsto:32908b7a702353b94cb595160a88ee05620015eb"

    label 'process_high'
    // nf-core labels
    label "process_medium_high_cpu"
    // Label for mpgi drac cpu alloc
    label "process_highest_memory"
    // Label for mpgi drac memory alloc
    label "process_medium_low_time"
    label "singleton"

    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai), path(bed)
    tuple val(meta2), path(ref)
    tuple val(meta3), path(fai)
    val model

    output:
    tuple val(meta), path("${meta.id}/snv*.vcf.gz"), emit: snv
    tuple val(meta), path("${meta.id}/indel*.vcf.gz"), emit: indel
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    run_clairs_to \\
        ${args} \\
        --threads ${task.ext.threads ?: task.cpus} \\
        -s ${prefix} \\
        --tumor_bam_fn ${bam} \\
        --ref_fn ${ref} \\
        --platform ${model} \\
        --bed_fn ${bed} \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ClairS-TO: \$(echo \$(run_clairs_to -v 2>&1) | awk '{ print \$2}')
    END_VERSIONS
    """
}

process DEEPSOMATIC_TUMORONLY {

    container "docker.io/google/deepsomatic:1.10.0"

    label 'process_high'
    // nf-core labels
    label "process_high_cpu"
    // Label for mpgi drac cpu alloc
    label "process_higher_memory"
    // Label for mpgi drac memory alloc
    label "process_medium_low_time"
    label "singleton"

    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai), path(bed)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: vcf_tbi
    tuple val(meta), path("${prefix}.g.vcf.gz"), emit: gvcf
    tuple val(meta), path("${prefix}.g.vcf.gz.tbi"), emit: gvcf_tbi
    tuple val("${task.process}"), val("deepsomatic"), val("1.10.0"), emit: versions_deepsomatic, topic: versions

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("DEEPSOMATIC module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: '--model_type=ONT_TUMOR_ONLY'
    prefix = task.ext.prefix ?: "${meta.id}"
    def regions = bed ? "--regions=${bed}" : ""
    """
    run_deepsomatic \\
        --ref=${fasta} \\
        --reads_tumor=${bam} \\
        --output_vcf=${prefix}.vcf.gz \\
        --output_gvcf=${prefix}.g.vcf.gz \\
        --sample_name_tumor="${prefix}" \\
        ${args} \\
        ${regions} \\
        --intermediate_results_dir=tmp \\
        --num_shards=${task.ext.threads ?: task.cpus}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("DEEPSOMATIC module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    echo "" | gzip > ${prefix}.g.vcf.gz
    touch ${prefix}.g.vcf.gz.tbi
    """
}

process CLASSY_COMBINED {
    container "${params.classy_container ?: 'ghcr.io/chusj-pigu/classy-slim:sha-71c699b'}"
    label "classy"
    label "process_low"
    label "process_medium_cpu"
    label "process_medium_mid_memory"
    label "process_low_time"

    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai), val(genome), path(ref)

    output:
    tuple val(meta), path("*combined_blood.svg"), emit: svg_blood
    tuple val(meta), path("*combined_brain.svg"), emit: svg_brain
    tuple val(meta), path("*combined_solid.svg"), emit: svg_solid
    tuple val(meta), path("*combined_classification.json"), emit: json_combined
    tuple val(meta), path("*combined_task_pies.svg"), emit: svg_combined
    tuple val(meta), path("*nanomix_pie.svg"), emit: svg_nanomix
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def threads = task.ext.threads ?: task.cpus
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ref_id = genome.id == "hs1" ? 't2t' : genome.id
    """
    classy combined \\
        -i ${bam} \\
        -o ${prefix}_combined_classification.json \\
        --sample ${prefix} \\
        --reference ${ref} \\
        --use-pileup \\
        --motif "CpG:CG" \\
        --genome ${ref_id} \\
        --pileup-threads ${threads} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        classy: Pre-release
        modkit: \$( modkit --version  | awk '{print \$2}')
        Combined: bundled
    END_VERSIONS
    """

    stub:
    """
    test -d /opt/classy/models/alma
    test -d /opt/classy/models/MPACT
    test -f /opt/classy/models/MPACT/runtime/model.json
    echo '{}' > ${meta.id}_combined_classification.json
    echo '{}' > ${meta.id}_combined_classification_lamprey.json
    echo '{}' > ${meta.id}_combined_classification_mpact.json
    touch ${meta.id}_combined_classification_mpact_mpact_pie.svg
    touch ${meta.id}_combined_classification_mpact_mpact_pie.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        classy: stub
    END_VERSIONS
    """
}

process SEVERUS_TUMOR_SINGLE {
    container "ghcr.io/chusj-pigu/severus:a294582458464a429891218e2be1c8c0db07b199"
    // TO DO: SET CONTAINER TO FIXED VERSION

    tag "${meta.id}"
    label 'process_medium'
    // nf-core labels
    label "process_medium_low_cpu"
    // Label for mpgi drac cpu alloc
    label "process_medium_memory"
    // Label for mpgi drac memory alloc
    label "process_low_time"

    input:
    tuple val(meta), path(bam), path(bai), path(snp_vcf)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(vntr)
    tuple val(meta4), path(pon)

    output:
    tuple val(meta), path("somatic_SVs/severus_somatic.vcf"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def threads = task.ext.threads ?: task.cpus
    def phased_vcf = snp_vcf ? "--phasing-vcf ${snp_vcf}" : ""
    """
    severus \\
        ${args} \\
        --target-bam ${bam} \\
        --out-dir '.' \\
        -t ${threads} \\
        ${phased_vcf} \\
        --vntr-bed ${vntr} \\
        --PON ${pon}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Severus: \$(severus -v)
    END_VERSIONS
    """
}

process STELLERATOR {

    container "ghcr.io/chusj-pigu/stellerator:5723fc7c4383bc34f07a9e7c159f9205ce1f8e08"

    label 'process_low'
    // nf-core labels
    label "process_low_cpu"
    // Label for mpgi drac cpu alloc
    label "process_medium_low_memory"
    // Label for mpgi drac memory alloc
    label "process_medium_low_time"

    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai)
    val genome
    path fusion_list

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val(meta), path("*.fasta.gz"), emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = task.ext.threads ?: task.cpus
    def gtf = "${genome}.ncbiRefSeq.gtf"
    """
    stellerator \\
        --bam ${bam} \\
        --annotation /opt/data/${gtf} \\
        --loci ${fusion_list} \\
        --output-vcf ${prefix}.vcf \\
        --threads ${threads} --verbose

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Stellerator : \$(stellerator --version | awk '{print \$2}')
    END_VERSIONS
    """
}
