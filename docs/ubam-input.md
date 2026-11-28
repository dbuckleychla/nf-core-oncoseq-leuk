Use `--input_ubam` to start from sample-separated, already basecalled unaligned BAMs. Put the BAM in the samplesheet **input** column; the existing **ubam** column remains reserved for resuming Dorado basecalling.

```csv
sample,input,tumor_type
SAMPLE_A,/path/to/sampleA.unaligned.bam,blood
SAMPLE_B,s3://your-input-bucket/sampleB.unaligned.bam,blood
```

The standard samplesheet metadata is preserved. Use `tumor_type` (not `sample_type`) to specify `blood` or another supported tissue type for classification; if omitted, it defaults to `other`. Other standard fields, including per-sample BED settings, remain available.

Each sample must have exactly one `.bam` input. Merge multiple same-sample uBAMs first. Inputs should already be demultiplexed and trimmed; this mode cannot be combined with `--demux` or `--skip_mapping`. Do not add `--skip_basecalling`: uBAM mode bypasses Dorado and its model downloads itself.

Add the new flag to your existing configured invocation. For example, after loading the usual Terraform/Nextflow environment:

```bash
nextflow run . -profile aws,clinical \
  --input ubam-inputs.csv \
  --input_ubam \
  --basecall_model 'dna_r10.4.1_e8.2_400bps_sup@v5.2.0' \
  --ref_cache "$REF_CACHE" \
  --vep_cache "$VEP_CACHE" \
  --run_id ubam-validation-001 \
  -resume
```

Set `basecall_model` to the model that generated the BAM; downstream caller model selection still needs this metadata. The Dorado GPU queue is unused, but GPU variant callers still use the general GPU queue. Existing mapping resource requests and CPU/BAM queues are unchanged.

Before alignment starts, the same task runs `samtools quickcheck -u` and checks the first 1,000 BAM records, or every record if fewer are present. Sampled records must be unmapped primary reads with sequences and base qualities. Every sampled record must contain:

- `MM` and `ML` for methylation analyses;
- `qs` for the default SUP read-quality filter;
- `st` when `--realtime` or `--time_series` is requested.

Empty BAMs and missing required tags fail before alignment. This checks tag presence, not tag calibration or the full contents of the BAM. If all sampled records pass, the rest of the BAM is assumed to have the same required tags, as requested. No read names or sequences are logged by the sampler.

The successful check writes `alignments/<sample>.ubam.validation.json`, recording the sample size, required tags and sampling assumption. AWS stages the input once for this mapping task; validation does not create another large BAM copy.

Alignment streams `samtools view → samtools reset → samtools fastq -T '*' → minimap2 -y → samtools sort`. No intermediate FASTQ is written for alignment. The existing optional realtime time-range reporting path still uses its legacy time-only FASTQ conversion; ordinary and time-series alignment remain streamed.

The default `qs ≥10` filter matches the clinical SUP path. Set `--ubam_min_qscore 9` for another appropriate cutoff, or `--ubam_min_qscore 0` to disable read-quality filtering and the `qs` requirement. `MM`/`ML` requirements remain. cfDNA length filtering reuses the existing aligned-BAM filtering path after mapping.

A small synthetic check of streaming alignment and tag retention can be run with:

```bash
nf-test test tests/subworkflows/local/mapping/ubam.nf.test --profile test,docker
```
