include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_SNP } from './processes'

workflow CHECK_ALIAS {
    main:
    BCFTOOLS_FILTER_SNP()
    emit:
    settings = BCFTOOLS_FILTER_SNP.out.settings
}
