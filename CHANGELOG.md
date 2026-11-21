# nf-core/oncoseq: Changelog

## 1.0.2-leuk.1-dev

- Add an optional offline molecular review report with exact allele annotation, COSMIC evidence, conservative breakpoint consensus and readable event plots.
- Add immutable reference preflight and complete native caller evidence exports.
- Retain both aligned BAMs and consolidated all-read uBAMs while suppressing shard, FASTQ and diagnostic BAM publication by default.

- Add deterministic POD5 file-group sharding, default 16 in clinical mode, and per-sample/project unaligned BAM consolidation before existing downstream processing.
- Add an independent institutional AWS Batch Terraform stack and launch profile.
- Align example capacity with the institutional stack and make GHCR authentication optional for public images.
- Extend exclusive BAM routing, add clinical step budgets and bounded resource retries, and set ClairS-TO/DeepSomatic indel minimum AF to 0.05.
- Reuse institutional FASTA, panel and Dorado model paths through separate cache parameters.
- Pin clinical candidate container digests and runtime versions, require versioned cached models, and correct S3 path handling.
- Document FrankONTstein settings and reporting differences and Dorado SUP v5.2.0 compatibility.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of nf-core/oncoseq, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
