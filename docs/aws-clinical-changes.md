# AWS Batch adaptation and clinical validation changes

This fork adapts Oncoseq release `v1.0.2` for AWS Batch execution and reproducible clinical validation. The changes address worker disk contention, software version drift, inconsistent failure handling and reference staging. The current implementation is a validation candidate; production acceptance still requires representative data and result concordance.

## Deployment and worker isolation

A standalone Terraform stack provisions CPU, BAM high/low, basecalling GPU and general GPU queues. Independent compute environments and capacity limits allow each workload to scale according to its requirements. Existing S3 storage holds work and results under a dedicated prefix. Deployment uses explicit CPU/NVIDIA AMI IDs, IAM roles and CloudWatch diagnostics.

Large BAM inputs and intermediates can exhaust a worker's disk when several jobs share it. The BAM queues therefore use fixed worker sizes and reserve every CPU for each assigned task:

POD5 basecalling can also stage large input collections. The clinical profile partitions existing POD5 files into up to 16 deterministic shards per sample/project, with independent Dorado tasks and subsequent unaligned BAM consolidation. The original downstream sample processing is retained. See [POD5 sharding](pod5-sharding.md) for partitioning and resume behaviour.

| Queue    | Worker                | Task reservation | Purpose                                                                                   |
| -------- | --------------------- | ---------------- | ----------------------------------------------------------------------------------------- |
| BAM high | 32 vCPUs, 256 GiB RAM | 32 vCPUs         | Alignment, sorting/merging, phasing, classification, large callers and VEP cache staging. |
| BAM low  | 8 vCPUs, 64 GiB RAM   | 8 vCPUs          | Lighter BAM/read processing, QC, CNA and plotting tasks.                                  |

Low-tier tasks move to the high queue when their actual memory request exceeds 48 GB, including on retries. Clinical configuration rejects missing high queues, partial CPU reservations and resource caps that defeat placement. Multi-shard Dorado jobs request 8 CPUs, 30 GB and one GPU each, allowing concurrent shards on separate GPUs in the same worker. Single-shard runs retain their larger reservation.

These reservations control placement. Local tool thread counts are configured independently so a full-worker reservation does not unnecessarily increase parsing or inference parallelism. Six process interfaces retain the pinned upstream inputs, outputs and image selections while supporting this thread control. Their source revisions are recorded in `modules/local/clinical_threads/upstream.json`; the original submodule pins remain intact.

High and low BAM disks default to 2000 GiB. Exclusive placement prevents concurrent staging, but each task's inputs, temporary files, outputs and Docker layers must still fit its disk. Fleet capacity increases the number of workers rather than the disk available to a single task.

## Resource budgets and failure handling

Clinical overrides increase memory for sorting, merging, phasing, Classy, Severus and VEP. Larger existing caller budgets are retained. Expensive steps receive the full configured time limit on every attempt, default 24 hours. Memory increases with retry attempts, honours the `huge` multiplier and remains within worker limits. Sorting memory is explicitly bounded per thread.

Failures remain visible. Exit 137 permits a memory retry only while the allocation can increase; exit 143 permits a bounded interruption retry. Other exits finish the workflow, including typical disk-full errors returning exit 1. Most steps allow one retry, and expensive steps allow two. Retries do not increase disk capacity. The trace records resource requests, queue, exit status, attempts and I/O counters to support diagnosis and tuning.

## Software and reference reproducibility

The `aws,clinical` profiles pin Nextflow 26.04.6, nf-schema 2.7.2 and nf-amazon 3.9.2. The release lock selects 36 Linux AMD64 container digests across 87 process names. Terraform 1.16.4 and AWS provider 6.13.0 are also pinned. Registry authentication is optional because the current pinned GHCR images permit anonymous pulls.

The `clinical` profile selects exact container builds through immutable SHA-256 image references. `images.lock.json` records the original release/build tags and resolved digests. Sources using moving tags such as `latest` are also frozen at their recorded digests, so later tag changes cannot alter a run. Changing an image requires a deliberate lock/config update and release validation. Runs must select the `clinical` profile to apply these pins.

Clinical runs require explicit, pre-staged reference and model assets. Separate `ref_cache`, `dorado_model_cache`, `reference_aux_cache` and `vep_cache` parameters accommodate existing storage layouts. Missing clinical assets cause an error rather than a download fallback. Existing bundled classifier, NASVAR and CNA resources are retained where applicable.

The hg38 preset uses an institutional mirror of the public VEP 115 cache. Set `vep_cache` to the parent of `homo_sapiens`, preserving `homo_sapiens/115_GRCh38/info.txt` below it. The public source metadata was verified anonymously; the mirror is read using the Batch job role. The auxiliary reference manifest records fixed source URLs for PharmCAT positions/index, tandem repeats and the structural-variant panel of normals. Operators must record reference/model checksums and immutable object versions; runtime checks currently establish presence and version selection rather than complete content identity.

## Explicit indel sensitivity

ClairS-TO receives `--indel_min_af 0.05`, and DeepSomatic receives `--make_examples_extra_args=vsc_min_fraction_indels=0.05`. These explicit candidate-generation thresholds avoid relying on container defaults. They can be adjusted through `clairsto_indel_min_af` and `deepsomatic_indel_min_af`. Their effect on sensitivity requires validation independently of downstream annotation and display filters.

## Validation status

Focused checks cover worker routing, input/alias coverage, retry promotion and ceilings, cache lookup, resource scaling, image locks and wrapper argument generation. The recorded suite comprises 20 mocked Terraform cases, 15 Python tests and 35 nf-test cases. Formatting and Nextflow syntax checks pass. Mock callers establish command contracts rather than inference correctness.

Remaining acceptance work includes deployment and concurrent-job placement checks, largest-sample disk/RAM measurements, GPU model execution, reference content verification and biological concordance. Dorado's pinned release catalogue supports SUP v5.2.0; GPU execution, including Blackwell compatibility, remains unverified. The current Clair3 selection still uses v5.0 weights and omits move emission, so complete SUP v5.2.0 model pairing remains unresolved. Missing integration fixtures and broader pipeline lint failures also remain open release gates.
