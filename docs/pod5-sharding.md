# POD5 basecalling shards

`basecall_shards` divides each sample's existing POD5 files into deterministic file groups. The `clinical` profile defaults to 16 groups; other profiles default to 1. The command-line value overrides either default:

```bash
nextflow run . -profile aws,clinical \
  -params-file assets/institutional_hg38.params.json \
  --input samplesheet.csv \
  --reference_aux_cache s3://YOUR_BUCKET/release/aux/ \
  --basecall_shards 16 \
  --run_id validation-001
```

The planner sorts resolved file URIs and distributes files round-robin across up to N nonempty groups. A directory input selects its direct-child `*.pod5` files. A single POD5 file remains one shard. Nested directories can be supplied as additional rows with the same sample identity. Rows for one sample must use consistent metadata, model and resume inputs. Duplicate resolved paths are rejected.

Each group becomes an independent Dorado task with a stable shard ID. POD5 files are staged under unique names in a shard-local input directory, supporting identical source basenames from different directories. This divides file staging and inference across jobs without rewriting POD5 data. It does not split reads inside a file or guarantee equal byte/read counts between groups; a large individual file must still fit its worker.

Shard tasks use the existing basecalling model, modification, device and trimming options. The AWS profile assigns each task to the basecalling queue with the existing CPU/GPU reservations. Shard count controls partitioning; fleet capacity controls simultaneous execution.

After every shard for an identity succeeds, unaligned BAMs are ordered by shard index and merged with SAMtools. Read-group/program headers are combined and per-read tags are retained. One-shard inputs pass through without a merge. The original metadata is restored before downstream processing.

For simplex inputs, consolidation is per sample before quality filtering and FASTQ conversion. For multiplex inputs, consolidation is per project before the existing demultiplexing, trimming, quality filtering and sample assignment. Those downstream operations retain their existing logic.

The merge uses the exclusive BAM high worker so all shard BAMs and the consolidated output are staged without concurrent BAM jobs on that instance. Its disk must fit that complete merge workload. Sharding reduces per-GPU input staging but does not reduce the sample's total BAM volume.

## Failure and resume

Incomplete or duplicate shard results cannot produce a valid consolidated BAM. Normal clinical failure/retry policies apply to individual shard tasks. `-resume` reuses completed shard and merge tasks when their inputs and commands match. Changing the partition changes the task inputs and cache identities.

A single partial resume uBAM cannot be supplied to several shards because it could duplicate existing reads in the merged output. Multi-shard plans reject that input. Use `basecall_shards = 1` for manually supplied Dorado resume BAMs, or use Nextflow task resume for sharded runs.

Focused tests check deterministic planning, duplicate input/resume rejection, sample separation, fewer-files-than-shards handling, incomplete-shard failure, and real SAMtools merging of valid fixture BAMs with read-group and modification tags. Dorado is mocked in these tests; representative GPU basecalling and result concordance remain required.
