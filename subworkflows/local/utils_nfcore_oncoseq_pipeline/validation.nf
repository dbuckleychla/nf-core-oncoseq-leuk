def validateExecutionOptions(options, profiles, workDirectory) {
    if (profiles.tokenize(',').contains('aws')) {
        ['aws_region', 'aws_queue', 'aws_job_role'].each { name ->
            if (!options[name]) {
                throw new IllegalArgumentException("AWS execution requires --${name} from the Terraform outputs")
            }
        }
        def outputDirectory = options.outdir instanceof Path ? options.outdir.toUriString() : options.outdir?.toString()
        if (!(workDirectory ==~ /s3:\/\/[^\/]+\/.+/) || !(outputDirectory ==~ /s3:\/\/[^\/]+\/.+/)) {
            throw new IllegalArgumentException('AWS execution requires S3 work and output prefixes. Load the Terraform environment and set --run_id, or supply -work-dir and --outdir')
        }
        if (!options.skip_basecalling && !options.skip_mapping && !options.aws_basecall_queue) {
            throw new IllegalArgumentException('Dorado basecalling requires --aws_basecall_queue')
        }
        if (options.aws_bam_low_queue && !options.aws_bam_high_queue) {
            throw new IllegalArgumentException('The BAM low queue requires --aws_bam_high_queue for larger tasks')
        }
    }

    if (options.clinical) {
        if (!options.ref_cache || !options.vep_cache) {
            throw new IllegalArgumentException('Clinical execution requires explicit --ref_cache and --vep_cache release assets')
        }
        if (!options.skip_basecalling && !options.skip_mapping) {
            if (!(options.basecall_model ==~ /dna_.+@v[0-9]+(\.[0-9]+)*/)) {
                throw new IllegalArgumentException('Clinical basecalling requires a complete versioned --basecall_model name')
            }
            if (options.m_bases && !(options.m_bases ==~ /.+@v[0-9]+(\.[0-9]+)*/)) {
                throw new IllegalArgumentException('Clinical basecalling requires a versioned --m_bases model')
            }
        }
    }
    return true
}
