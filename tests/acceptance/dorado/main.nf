include { DORADO_BASECALL } from '../../../modules/local/dorado_sharded/main'

process VALIDATE_DORADO {
    container 'community.wave.seqera.io/library/htslib_samtools@sha256:a55ddea590e567a91df592300a960aa534cfc1bd16e7623e3938ec21f4f3df15'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(bam)
    path versions

    output:
    path 'validated.tsv', emit: validation

    script:
    """
    grep -Eq 'dorado:[[:space:]]+2\\.2\\.0([+[:space:]]|\$)' ${versions}
    samtools quickcheck -u ${bam}
    reads=\$(samtools view -c ${bam})
    test "\$reads" -gt 0 && test "\$reads" -le 32
    samtools view -H ${bam} > header.sam
    grep -F 'dna_r10.4.1_e8.2_400bps_sup@v5.2.0' header.sam
    samtools view ${bam} | awk 'index(\$0, "\\tMM:Z:") && index(\$0, "\\tML:B:C,") { n++ } END { if (n == 0) exit 1 }'
    printf 'dorado_version\\treads\\tmodification_tags\\n2.2.0\\t%s\\tpresent\\n' "\$reads" > validated.tsv
    """
}

workflow DORADO_ACCEPTANCE {
    main:
    if (!(params.image ==~ /^ghcr\.io\/dbuckleychla\/dorado-2\.2\.0@sha256:[a-f0-9]{64}$/)) {
        error 'Supply the published Dorado 2.2.0 Linux amd64 digest with --image; tags are not accepted'
    }
    if (!params.pod5 || !params.model || !params.mod_model || !params.outdir) {
        error 'Supply --pod5, --model, --mod_model and --outdir'
    }
    if (workflow.profile.tokenize(',').contains('aws') && (!params.aws_basecall_queue || !params.aws_queue)) {
        error 'AWS acceptance requires --aws_basecall_queue and --aws_queue from the existing deployment'
    }
    model = file(params.model, checkIfExists: true)
    mod_model = file(params.mod_model, checkIfExists: true)
    if (model.name != 'dna_r10.4.1_e8.2_400bps_sup@v5.2.0' || mod_model.name != 'dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mCG_5hmCG@v2') {
        error 'Acceptance requires the existing SUP v5.2.0 and matching CpG v2 model directories'
    }
    ch_input = channel.of(tuple(
        [id: 'SYNTHETIC_DORADO_SMOKE'],
        file(params.pod5, checkIfExists: true),
        file("${moduleDir}/../../../assets/NOFILE", checkIfExists: true),
        model,
        mod_model,
    ))
    DORADO_BASECALL(ch_input)
    VALIDATE_DORADO(DORADO_BASECALL.out.ubam, DORADO_BASECALL.out.versions)

    emit:
    bam = DORADO_BASECALL.out.ubam
    validation = VALIDATE_DORADO.out.validation
}

workflow {
    DORADO_ACCEPTANCE()
}
