def qcQuote(value) {
    "'" + value.toString().replace("'", "'\\''") + "'"
}

def qcPy(value) {
    groovy.json.JsonOutput.toJson(value.toString())
}

process QC_EXT_REFERENCE {
    tag 'reference'
    label 'qc_extended_small'
    container 'ghcr.io/dbuckleychla/frankontstein@sha256:15c315c121d37b5c115d3d965836e60d56756498ca2deebe066a373d41e96785'
    input:
    tuple val(meta), path(fasta), path(fai)
    val code_digest
    output:
    tuple path('reference.json'), path('reference.cpgs.gz'), emit: summary
    script:
    """
    extended_qc.py reference --fasta ${qcQuote(fasta)} --cpgs reference.cpgs.gz --output reference.json
    """
}

process QC_EXT_ALIGNMENT {
    tag "${meta.id}"
    label 'qc_extended_bam'
    container 'ghcr.io/dbuckleychla/frankontstein@sha256:15c315c121d37b5c115d3d965836e60d56756498ca2deebe066a373d41e96785'
    input:
    tuple val(meta), path(bam), path(bai), path(enrichment, stageAs: 'enrichment/*'), path(targets, stageAs: 'targets/*')
    tuple path(reference), path(cpgs)
    val code_digest
    output:
    tuple val(meta), path(bam), path(bai), path(enrichment), path(targets), path('alignment.json'), emit: alignment
    script:
    def disabled = params.qc_methylation ? '' : '--disable-methylation'
    """
    extended_qc.py alignment --bam ${qcQuote(bam)} --reference ${qcQuote(reference)} \\
        --enrichment ${qcQuote(enrichment)} --targets ${qcQuote(targets)} --threads ${Math.min(8, task.cpus as int)} \\
        ${disabled} --output alignment.json
    """
}

process QC_EXT_MODKIT_PILEUP {
    tag "${meta.id}"
    label 'qc_extended_bam'
    container 'ghcr.io/chusj-pigu/classy-slim@sha256:7ec5fc1e4607be13d30dc65a1b07f2844f8289edbec0ca67a70e9b88ef8beb4c'
    input:
    tuple val(meta), path(bam), path(bai), path(enrichment, stageAs: 'enrichment/*'), path(targets, stageAs: 'targets/*'), path(alignment)
    tuple val(ref_meta), path(fasta), path(fai)
    val code_digest
    output:
    tuple val(meta), path('pileup.bed.gz'), path('modkit.version.txt'), emit: pileup
    script:
    """
    python3 - <<'PY'
    import gzip, json, subprocess
    from pathlib import Path
    data = json.loads(Path(${qcPy(alignment)}).read_text())
    if data['methylation_status']['status'] in ('disabled', 'unavailable'):
        with gzip.open('pileup.bed.gz', 'wt'): pass
        Path('modkit.version.txt').write_text('not_run\\n')
    else:
        subprocess.run(['modkit', 'pileup', ${qcPy(bam)}, 'pileup.bed.gz', '--ref', ${qcPy(fasta)},
                        '--cpg', '--modified-bases', '5mC', '5hmC', '--combine-strands', '--bgzf',
                        '--threads', '${Math.min(8, task.cpus as int)}', '--num-reads', '10042', '--filter-percentile', '0.1',
                        '--seed', '42', '--max-depth', '8000', '--log-filepath', 'modkit.log'], check=True)
        Path('modkit.version.txt').write_text(subprocess.check_output(['modkit', '--version'], text=True))
    PY
    """
}

process QC_EXT_METHYLATION {
    tag "${meta.id}"
    label 'qc_extended_bam'
    container 'ghcr.io/dbuckleychla/frankontstein@sha256:15c315c121d37b5c115d3d965836e60d56756498ca2deebe066a373d41e96785'
    input:
    tuple val(meta), path(enrichment, stageAs: 'enrichment/*'), path(targets, stageAs: 'targets/*'), path(alignment), path(pileup), path(modkit_version)
    tuple path(reference), path(cpgs)
    tuple val(ref_meta), path(fasta), path(fai)
    val panel_rows
    path panel_files, stageAs: 'panels/*'
    val assembly
    val code_digest
    output:
    tuple val(meta), path('methylation.json'), emit: metrics
    script:
    def manifest = (['panel_id\tassembly\tmodel_name\tmodel_version\tbed'] + panel_rows.collect { row ->
        [row.panel_id, row.assembly, row.model_name, row.model_version, "panels/${row.bed_name}"].join('\t')
    }).join('\n')
    def versionPath = modkit_version ? qcPy(modkit_version) : 'None'
    """
    cat > panels.tsv <<'QC_PANELS'
    ${manifest}
    QC_PANELS
    extended_qc.py methylation --alignment ${qcQuote(alignment)} --reference ${qcQuote(reference)} \\
        --cpgs ${qcQuote(cpgs)} --fasta ${qcQuote(fasta)} --bedmethyl ${qcQuote(pileup)} \\
        --enrichment ${qcQuote(enrichment)} --targets ${qcQuote(targets)} --assembly ${qcQuote(assembly)} \\
        --panel-manifest panels.tsv --output methylation.json
    python3 - <<'PY'
    import json
    from pathlib import Path
    data = json.loads(Path('methylation.json').read_text())
    version = ${versionPath}
    data['modkit_version'] = Path(version).read_text().strip() if version else 'not_run'
    Path('methylation.json').write_text(json.dumps(data, allow_nan=False) + '\\n')
    PY
    """
}

process QC_EXT_PREPROCESSING {
    tag "${scope.id}:${scope.stage}"
    label 'qc_extended_bam'
    container 'ghcr.io/dbuckleychla/frankontstein@sha256:15c315c121d37b5c115d3d965836e60d56756498ca2deebe066a373d41e96785'
    input:
    tuple val(scope), path(bams, stageAs: 'stage_inputs/*')
    val code_digest
    output:
    tuple val(scope), path('stage.json'), emit: metrics
    script:
    def inputs = bams instanceof Collection ? bams : [bams]
    def descriptor = groovy.json.JsonOutput.toJson([[scope_type: scope.scope_type, scope_id: scope.id,
        stage: scope.stage, format: scope.format, files: inputs.collect { it.toString() }]])
    """
    cat > stage_inputs.json <<'QC_INPUTS'
    ${descriptor}
    QC_INPUTS
    extended_qc.py preprocessing --input stage_inputs.json --output stage.json
    """
}

process QC_EXT_ASSEMBLE {
    tag "${meta.id}"
    label 'qc_extended_small'
    container 'ghcr.io/dbuckleychla/frankontstein@sha256:15c315c121d37b5c115d3d965836e60d56756498ca2deebe066a373d41e96785'
    input:
    tuple val(meta), path(enrichment, stageAs: 'enrichment/*'), path(targets, stageAs: 'targets/*'), path(alignment), path(methylation)
    tuple path(reference), path(cpgs)
    val stages
    val assembly
    val code_digest
    output:
    tuple val(meta), path('metrics/metrics.json'), emit: metrics
    tuple val(meta), path('metrics/*.tsv'), emit: tables
    script:
    def preprocessing = groovy.json.JsonOutput.toJson(stages.findAll { it.scope_type == 'sample' && it.scope_id == meta.id })
    """
    cat > preprocessing.json <<'QC_STAGES'
    ${preprocessing}
    QC_STAGES
    extended_qc.py assemble --sample ${qcQuote(meta.id)} --alignment ${qcQuote(alignment)} \\
        --methylation ${qcQuote(methylation)} --reference ${qcQuote(reference)} --assembly ${qcQuote(assembly)} \\
        --enrichment ${qcQuote(enrichment)} --targets ${qcQuote(targets)} --preprocessing preprocessing.json \\
        --code-digest ${qcQuote(code_digest)} --output metrics
    """
}

process QC_EXT_REPORT {
    tag "${meta.id}"
    label 'qc_extended_small'
    container 'ghcr.io/dbuckleychla/frankontstein@sha256:15c315c121d37b5c115d3d965836e60d56756498ca2deebe066a373d41e96785'
    input:
    tuple val(meta), path(metrics)
    val code_digest
    output:
    tuple val(meta), val('Extended_QC'), path('qc_report_inputs'), val('Extended_QC.qmd'), emit: section
    tuple val(meta), path('content.html'), emit: content
    script:
    """
    mkdir qc_report_inputs
    extended_qc.py report --input ${qcQuote(metrics)} --output content.html
    python3 - <<'PY'
    from pathlib import Path
    Path('qc_report_inputs/Extended_QC.qmd').write_text('## Extended sequencing QC {#extended-qc}\\n\\n```{=html}\\n' + Path('content.html').read_text() + '\\n```\\n')
    PY
    """
}

process QC_EXT_ATTACH_REPORT {
    tag "${meta.id}"
    label 'qc_extended_small'
    container 'ghcr.io/dbuckleychla/frankontstein@sha256:15c315c121d37b5c115d3d965836e60d56756498ca2deebe066a373d41e96785'
    input:
    tuple val(meta), path(native_html, stageAs: 'native_input/*'), path(fragment, stageAs: 'qc_input/*')
    output:
    tuple val(meta), path("${meta.id}.html"), emit: html
    script:
    """
    extended_qc.py attach-report --input ${qcQuote(native_html)} --fragment ${qcQuote(fragment)} --output ${qcQuote(meta.id + '.html')}
    """
}

process QC_EXT_RUN_SUMMARY {
    tag 'run'
    label 'qc_extended_small'
    container 'ghcr.io/dbuckleychla/frankontstein@sha256:15c315c121d37b5c115d3d965836e60d56756498ca2deebe066a373d41e96785'
    input:
    path tables, stageAs: 'samples/??/*'
    val stages
    val code_digest
    output:
    path('run_summary.tsv'), emit: summary
    path('run_preprocessing.tsv'), emit: preprocessing
    script:
    def rows = groovy.json.JsonOutput.toJson(stages.findAll { it.scope_type == 'run' })
    """
    extended_qc.py run-summary --input samples --output run_summary.tsv
    python3 - <<'PY'
    import csv, json
    rows = json.loads(${qcPy(rows)})
    with open('run_preprocessing.tsv', 'w') as out:
        writer = csv.DictWriter(out, fieldnames=['scope_type','scope_id','stage','source_kind','reads','bases','exhaustive'], delimiter='\\t')
        writer.writeheader()
        writer.writerows(rows)
    PY
    """
}
