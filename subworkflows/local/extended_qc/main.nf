include { QC_EXT_REFERENCE; QC_EXT_ALIGNMENT; QC_EXT_MODKIT_PILEUP; QC_EXT_METHYLATION; QC_EXT_PREPROCESSING; QC_EXT_ASSEMBLE; QC_EXT_REPORT; QC_EXT_RUN_SUMMARY } from '../../../modules/local/extended_qc/main'

def extendedQcDigest() {
    def root = file("${moduleDir}/../../..").normalize()
    def digest = java.security.MessageDigest.getInstance('SHA-256')
    (files("${root}/bin/qc_extended/*.py") + [file("${root}/bin/extended_qc.py")]).sort { it.name }.each { digest.update(it.bytes) }
    digest.digest().encodeHex().toString()
}

def qcManifest(path, columns) {
    if (!path) return []
    def lines = file(path, checkIfExists: true).readLines().findAll { it.trim() }
    if (!lines || lines[0].split('\t', -1).toList() != columns) {
        throw new IllegalArgumentException('Invalid QC manifest columns: ' + columns.join(','))
    }
    lines.drop(1).collect { line ->
        def values = line.split('\t', -1).toList()
        if (values.size() != columns.size() || values.any { !it }) throw new IllegalArgumentException('Incomplete QC manifest row')
        [columns, values].transpose().collectEntries { key, value -> [key, value] }
    }
}

def qcAsset(value, manifest) {
    // Typed remote parameters already carry their filesystem provider.
    if (value instanceof java.nio.file.Path) return file(value, checkIfExists: true)
    def uri = value.toString()
    file(uri.contains('://') || uri.startsWith('/') || !manifest ? uri : file(manifest).parent.resolve(uri), checkIfExists: true)
}

workflow EXTENDED_QC {
    take:
    bams
    reference
    fai
    regions
    preprocessing
    assembly

    main:
    def overrides = qcManifest(params.qc_sample_manifest, ['sample','targets_bed'])
    if (overrides.collect { it.sample }.unique(false).size() != overrides.size()) throw new IllegalArgumentException('Duplicate QC target override')
    def panels = qcManifest(params.qc_cpg_panels, ['panel_id','assembly','model_name','model_version','bed'])
    def panelPaths = panels.collect { qcAsset(it.bed, params.qc_cpg_panels) }
    if (panelPaths.collect { it.name }.unique(false).size() != panelPaths.size()) throw new IllegalArgumentException('QC panel filenames must be unique')
    panels = panels.withIndex().collect { row, i -> row + [bed_name: panelPaths[i].name] }
    def code = channel.value(extendedQcDigest())
    def ref = reference.join(fai).first()
    QC_EXT_REFERENCE(ref, code)
    def selected = bams.collect(flat: false).flatMap { rows ->
        def ids = rows.collect { it[0].id }
        if (ids.unique(false).size() != ids.size() || overrides.any { !ids.contains(it.sample) }) throw new IllegalArgumentException('Duplicate sample or unknown QC override')
        rows.collect { row ->
            def selected = overrides.find { it.sample == row[0].id }
            def target = selected ? qcAsset(selected.targets_bed, params.qc_sample_manifest) : params.qc_targets_bed ? qcAsset(params.qc_targets_bed, null) : null
            if (!target) throw new IllegalArgumentException('Extended QC requires an explicit exact-target BED: --qc_targets_bed or --qc_sample_manifest')
            if (row[0].id ==~ /.*[\/\\\x00-\x1f].*/ || row[0].id in ['.','..']) throw new IllegalArgumentException('Unsafe QC sample identifier')
            tuple(row[0].id, row[0], row[1], row[2], target)
        }
    }
    def inputs = selected.join(regions.map { meta, bed, _padding, _lf -> tuple(meta.id, bed) })
        .map { _id, meta, bam, index, targets, enrichment -> tuple(meta, bam, index, enrichment, targets) }
    QC_EXT_ALIGNMENT(inputs, QC_EXT_REFERENCE.out.summary, code)
    def eligible = QC_EXT_ALIGNMENT.out.alignment.branch { _meta, _bam, _index, _enrichment, _targets, alignment ->
        def status = new groovy.json.JsonSlurperClassic().parseText(alignment.text).methylation_status.status
        ready: status in ['complete', 'partial']
        absent: true
    }
    QC_EXT_MODKIT_PILEUP(eligible.ready, ref, code)
    def pileups = QC_EXT_MODKIT_PILEUP.out.pileup.mix(eligible.absent.map { meta, _bam, _index, _enrichment, _targets, _alignment -> tuple(meta, [], []) })
    def methylation = QC_EXT_ALIGNMENT.out.alignment.map { meta, _bam, _index, enrichment, targets, alignment -> tuple(meta.id, meta, enrichment, targets, alignment) }
        .join(pileups.map { meta, pileup, version -> tuple(meta.id, pileup, version) })
        .map { _id, meta, enrichment, targets, alignment, pileup, version -> tuple(meta, enrichment, targets, alignment, pileup, version) }
    QC_EXT_METHYLATION(methylation, QC_EXT_REFERENCE.out.summary, ref, channel.value(panels), channel.value(panelPaths), channel.value(assembly), code)
    QC_EXT_PREPROCESSING(preprocessing.groupTuple().map { scope, stageFiles -> tuple(scope, stageFiles.flatten()) }, code)
    def stages = QC_EXT_PREPROCESSING.out.metrics.map { _scope, json ->
        new groovy.json.JsonSlurperClassic().parseText(json.text)
    }.toList().map { it.flatten() }
    def assembled = QC_EXT_ALIGNMENT.out.alignment.map { meta, _bam, _index, enrichment, targets, alignment -> tuple(meta.id, meta, enrichment, targets, alignment) }
        .join(QC_EXT_METHYLATION.out.metrics.map { meta, json -> tuple(meta.id, json) })
        .map { _id, meta, enrichment, targets, alignment, methylationJson -> tuple(meta, enrichment, targets, alignment, methylationJson) }
    QC_EXT_ASSEMBLE(assembled, QC_EXT_REFERENCE.out.summary, stages, channel.value(assembly), code)
    QC_EXT_REPORT(QC_EXT_ASSEMBLE.out.metrics, code)
    QC_EXT_RUN_SUMMARY(QC_EXT_ASSEMBLE.out.tables.map { _meta, files -> files }.collect().map { it.flatten() }, stages, code)

    emit:
    metrics = QC_EXT_ASSEMBLE.out.metrics
    sections = QC_EXT_REPORT.out.section
    report_content = QC_EXT_REPORT.out.content
    tables = QC_EXT_ASSEMBLE.out.tables
}
