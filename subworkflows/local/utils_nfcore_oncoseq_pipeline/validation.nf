def validateExecutionOptions(options, profiles, workDirectory) {
    if (options.input_ubam) {
        if (options.skip_mapping || options.demux) {
            throw new IllegalArgumentException('--input_ubam currently requires sample-separated inputs and cannot be combined with --skip_mapping or --demux')
        }
        if (options.ubam_min_qscore == null || options.ubam_min_qscore < 0) {
            throw new IllegalArgumentException('--ubam_min_qscore must be non-negative')
        }
    }
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
        if (!options.input_ubam && !options.skip_basecalling && !options.skip_mapping && !options.aws_basecall_queue) {
            throw new IllegalArgumentException('Dorado basecalling requires --aws_basecall_queue')
        }
        if (options.aws_bam_low_queue && !options.aws_bam_high_queue) {
            throw new IllegalArgumentException('The BAM low queue requires --aws_bam_high_queue for larger tasks')
        }
        if (options.clinical) {
            if (!options.aws_bam_high_queue) {
                throw new IllegalArgumentException('Clinical AWS execution requires --aws_bam_high_queue for isolated BAM staging')
            }
            if ((options.max_cpus ?: 0) < 32 || options.aws_bam_high_vcpus != 32 || (options.aws_bam_low_queue && options.aws_bam_low_vcpus != 8)) {
                throw new IllegalArgumentException('Clinical AWS execution requires max_cpus >= 32 and full BAM reservations of 32 high / 8 low vCPUs')
            }
            if (options.max_memory && options.max_memory > 224.GB) {
                throw new IllegalArgumentException('Clinical AWS execution requires max_memory <= 224 GB to leave host headroom on the BAM high workers')
            }
        }
    }

    if (options.clinical) {
        if (!options.ref_cache || !options.vep_cache) {
            throw new IllegalArgumentException('Clinical execution requires explicit --ref_cache and --vep_cache release assets')
        }
        if (!options.input_ubam && !options.skip_basecalling && !options.skip_mapping) {
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
