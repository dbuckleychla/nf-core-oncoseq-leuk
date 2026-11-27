def bloodQuote(value) {
    "'" + value.toString().replace("'", "'\\''") + "'"
}

process CLASSY_BLOOD {
    tag "${meta.id}"
    label 'classy_blood'
    label 'process_medium_cpu'
    label 'process_medium_mid_memory'
    label 'process_low_time'
    container '288525090463.dkr.ecr.us-west-2.amazonaws.com/oncoseq-leuk/classy-blood-nextflow@sha256:77aed6f284d94eda23f27a75080919a7b6888accb75cbbe73d7fd08c72121561'

    input:
    tuple val(meta), path(bam), path(bai), val(genome), path(ref)

    output:
    tuple val(meta), path('*_combined_classification.json'), val(task.container), emit: json
    tuple val(meta), path('*_combined_blood.svg'), emit: plot
    tuple val(meta), path('*_combined_task_pies.svg'), emit: summary
    path 'versions.yml', emit: versions

    script:
    def ref_id = genome.id == 'hs1' ? 't2t' : genome.id
    def threads = task.ext.threads ?: Math.min(8, task.cpus as int)
    def short_reads = params.cfdna ? '--min-read-length 20' : ''
    """
    classy combined --input ${bloodQuote(bam)} --output ${bloodQuote(meta.id + "_combined_classification.json")} \\
        --sample ${bloodQuote(meta.id)} --reference ${bloodQuote(ref)} --genome ${bloodQuote(ref_id)} \\
        --use-pileup --motif CpG:CG --pileup-threads ${threads} \\
        --exclude crossnn,tucan,sturgeon,alma,mpact,nanomix \\
        --enable-alma3 ${params.classy_blood_enable_alma3} \\
        --alma3-modification-mode ${params.classy_blood_alma3_modification_mode} \\
        --threshold-figure-colors ${short_reads}
    cat > versions.yml <<END_VERSIONS
    "${task.process}":
        classy: \$(classy --version | awk '{print \$2}')
        modkit: \$(modkit --version | awk '{print \$2}')
    END_VERSIONS
    """
}

process CLASSY_BLOOD_TABLES {
    tag "${meta.id}"
    label 'process_low'
    container 'community.wave.seqera.io/library/oncoseq_review_1_0_0@sha256:b143e350e1dae33ad26758fc6f9ccf2392cb8b49881de73960f42525602d5eb5'

    input:
    tuple val(meta), path(json), path(plots, stageAs: 'plots/*'), val(classifier_image)

    output:
    tuple val(meta), path('blood_metrics/review_classification.json'), emit: review_json
    tuple val(meta), path('blood_metrics/tables/*.tsv'), emit: tables
    tuple val(meta), path('blood_metrics/classification.metadata.json'), emit: metadata
    tuple val(meta), val('Methylation'), path('blood_metrics/native_report_inputs'), val("${meta.id}-Methylation.qmd"), emit: section

    script:
    """
    classy_blood_tables.py --input ${bloodQuote(json)} --sample ${bloodQuote(meta.id)} --plots plots --output blood_metrics \\
        --enable-alma3 ${params.classy_blood_enable_alma3} \\
        --modification-mode ${params.classy_blood_alma3_modification_mode} \\
        --image ${bloodQuote(classifier_image)}
    """
}
