include { CLASSY_BLOOD; CLASSY_BLOOD_TABLES } from '../../../modules/local/classy_blood/main'

workflow CLASSY_BLOOD_WORKFLOW {
    take:
    bams
    reference

    main:
    CLASSY_BLOOD(bams.combine(reference))
    ch_tables = CLASSY_BLOOD.out.json.join(CLASSY_BLOOD.out.plot).join(CLASSY_BLOOD.out.summary)
        .map { meta, json, image, plot, summary -> tuple(meta, json, [plot, summary], image) }
    CLASSY_BLOOD_TABLES(ch_tables)

    emit:
    json = CLASSY_BLOOD_TABLES.out.review_json
    native_json = CLASSY_BLOOD.out.json
    plots = CLASSY_BLOOD.out.plot.mix(CLASSY_BLOOD.out.summary)
    sections = CLASSY_BLOOD_TABLES.out.section
    metadata = CLASSY_BLOOD_TABLES.out.metadata
    versions = CLASSY_BLOOD.out.versions
}
