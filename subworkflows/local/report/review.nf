include { orderedReviewArtifacts } from './input_order'
include { REVIEW_PREPARE ; REVIEW_VEP ; REVIEW_COMPLETE ; REVIEW_FIGURES ; REVIEW_RENDER ; REVIEW_PROVENANCE_ATTACH } from '../../../modules/local/review/prepare/main'

def resolveReviewReferences(source) {
    source = source ?: "${projectDir}/assets/institutional_hg38.review.json"
    def manifestFile = file(source, checkIfExists: true)
    def manifest = new groovy.json.JsonSlurperClassic().parseText(manifestFile.text)
    if (manifest.schema_version != 1 || manifest.assembly != 'GRCh38') {
        throw new IllegalArgumentException('Review references require schema_version=1 and assembly=GRCh38')
    }
    def required = ['fasta', 'fai', 'targets', 'gff', 'gnomad_genome', 'gnomad_exome', 'clinvar', 'clinvar_common_pathogenic', 'cosmic_targeted', 'cosmic_evidence']
    if (!manifest.resources || required.any { name -> !manifest.resources.containsKey(name) }) {
        throw new IllegalArgumentException('Review references are missing required resources')
    }
    def paths = []
    def names = [] as Set
    manifest.resources.each { _name, entry ->
        def refUri = entry.uri.toString()
        def ref = file(refUri.contains('://') || refUri.startsWith('/') ? refUri : manifestFile.parent.resolve(refUri), checkIfExists: true)
        entry.uri = ref.toUriString()
        if (!names.add(ref.name)) {
            throw new IllegalArgumentException('Review reference filenames must be unique: ' + ref.name)
        }
        paths << ref
        entry.path = "references/${ref.name}"
        if (entry.index) {
            def indexUri = entry.index.uri.toString()
            def index = file(indexUri.contains('://') || indexUri.startsWith('/') ? indexUri : manifestFile.parent.resolve(indexUri), checkIfExists: true)
            entry.index.uri = index.toUriString()
            if (!names.add(index.name)) {
                throw new IllegalArgumentException('Duplicate review index filename: ' + index.name)
            }
            paths << index
            entry.index.path = "references/${index.name}"
        }
    }
    return [manifest, paths]
}

def reviewCodeDigest() {
    def paths = files("${projectDir}/bin/review/*.py") + [file("${projectDir}/bin/qc_extended/report.py"), file("${projectDir}/bin/review_report.py"), file("${projectDir}/bin/review_figures.py")]
    def digest = java.security.MessageDigest.getInstance('SHA-256')
    paths.sort { left, right -> left.name <=> right.name }.each { path -> digest.update(path.bytes) }
    digest.digest().encodeHex().toString()
}

workflow REVIEW_REPORT {
    take:
    artifacts
    bams
    ref
    vep_cache
    run_provenance
    software_versions

    main:
    def resolved = resolveReviewReferences(params.review_references)
    def refs = channel.value(resolved[0])
    def code = channel.value(reviewCodeDigest())
    def paths = channel.value(resolved[1])
    def panel = channel.value(file(params.review_panel, checkIfExists: true))
    def known = channel.value(file("${projectDir}/assets/fusion-list.txt", checkIfExists: true))
    def template = channel.value(file("${projectDir}/assets/review/report.html", checkIfExists: true))
    def reference = ref.first()

    ch_samples = artifacts
        .groupTuple()
        .map { sample, kinds, files ->
            def ordered = orderedReviewArtifacts(kinds, files)
            tuple([id: sample], ordered[0], ordered[1])
        }

    REVIEW_PREPARE(ch_samples, paths, refs, panel, run_provenance, code)
    ch_prepared = REVIEW_PREPARE.out.prepared.map { meta, catalogue, payload ->
        def counts = new groovy.json.JsonSlurperClassic().parseText(file("${catalogue.toUriString()}/summary.json").text)
        tuple(meta + [review_candidates: counts.candidate_alleles], catalogue, payload)
    }
    REVIEW_VEP(ch_prepared, reference, vep_cache, code, channel.value(resolved[0].resources.fasta.sha256 ?: ''))
    REVIEW_COMPLETE(REVIEW_VEP.out.annotated, paths, refs, panel, known, channel.value(file("${projectDir}/images.lock.json")), code)

    ch_data = REVIEW_COMPLETE.out.data.branch { _meta, data ->
        figures: new groovy.json.JsonSlurperClassic().parseText(file("${data.toUriString()}/review.json").text).structural_events.any { event -> event.selected }
        plain: true
    }
    ch_figures = ch_data.figures
        .map { meta, data -> tuple(meta.id, meta, data) }
        .join(bams.map { meta, bam, bai -> tuple(meta.id, bam, bai) })
        .map { _id, meta, data, bam, bai -> tuple(meta, data, bam, bai) }
    REVIEW_FIGURES(ch_figures, code)
    ch_render = ch_data.plain
        .map { meta, data -> tuple(meta, data, []) }
        .mix(REVIEW_FIGURES.out.figures)
    REVIEW_RENDER(ch_render, template, code)
    REVIEW_PROVENANCE_ATTACH(REVIEW_RENDER.out.html, software_versions, code)

    emit:
    html = REVIEW_PROVENANCE_ATTACH.out.html
    data = REVIEW_COMPLETE.out.data
}
