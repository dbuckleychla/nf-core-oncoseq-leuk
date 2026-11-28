# Validation-run resource review

Evaluated on 21 November 2026. No task requests, tool thread settings, queue
allocations or instance lists were changed. Performance changes require user
approval.

## Evidence and limits

The evidence is the four execution reports and associated traces under
`/Volumes/DNB_External_1/temp_SSD_data/validation-shards-001/pipeline_info/`,
from 20 November at 13:42, 19:14, 20:25 and 21:20. The last report contains
202 tasks, including cached tasks from earlier sessions. There are 213 unique
work hashes across the reports. Cached measurements describe the original
execution; they were not counted as additional executions.

This is one specimen with interrupted/resumed runs. Do not interpret the time
between sessions as uninterrupted pipeline runtime. Table times are measured
tool runtime, excluding queueing and most wrapper/staging overhead. Average
CPU usage is expressed in cores (100% = one core), not percentage of the full
allocation. RSS excludes GPU memory and is not a reliable measurement of all
filesystem cache or transient host/container memory pressure. Physical disk
bytes do not include every S3 transfer or cached read. There is no GPU
utilisation, disk-latency, scratch-occupancy or detailed stage timing in these
reports. The new optional review modules are absent from this run.

## Main findings

| Step                     | Tool runtime | Mean CPU cores | Peak RSS GiB | Request           | Interpretation                                                                      |
| ------------------------ | -----------: | -------------: | -----------: | ----------------- | ----------------------------------------------------------------------------------- |
| DeepSomatic              |     47.8 min |           10.3 |         12.1 | 32 CPUs / 100 GiB | Runs on CPU BAM queue, no GPU requested; tool has only 8 example-generation shards. |
| WhatsHap haplotag        |     46.0 min |            2.1 |          1.3 | 32 CPUs / 36 GiB  | Already has 32 output-compression threads; main scan/assignment dominates.          |
| Minimap2 + samtools sort |     34.4 min |           22.4 |         59.3 | 32 CPUs / 72 GiB  | Best CPU/memory benchmark candidate; tool threads are 32 + 31, sharing one worker.  |
| ClairS-TO                |     33.3 min |           11.5 |         16.7 | 32 CPUs / 180 GiB | Uses 16 tool threads; more memory is not justified by this specimen.                |
| ONTIME range filter      |     35.2 min |            1.0 |          1.8 | 8 CPUs / 10 GiB   | Serial timestamp filtering reads the large aligned BAM.                             |
| NASVAR                   |     28.6 min |            1.0 |          3.5 | 8 CPUs / 20 GiB   | Current pipeline command exposes no thread setting.                                 |
| WhatsHap phase           |     28.2 min |            1.0 |          1.4 | 32 CPUs / 32 GiB  | Core phasing is serial; 32-thread bgzip at the end is not core parallelism.         |
| Delly CNV                |     27.6 min |            0.8 |          1.9 | 8 CPUs / 20 GiB   | More reserved CPUs alone unlikely to help.                                          |
| Severus                  |     22.4 min |            0.6 |         34.1 | 32 CPUs / 64 GiB  | Uses 4 tool threads, but low average utilisation suggests a serial/blocked phase.   |
| QDNAseq                  |     17.5 min |            1.0 |          4.0 | 8 CPUs / 20 GiB   | Mostly serial; more RAM/CPUs is not an evidenced remedy.                            |
| Clair3                   |      9.6 min |           11.3 |         25.7 | 32 CPUs / 100 GiB | Relatively short; lower priority than DeepSomatic and mapping.                      |

The phase → haplotag → Severus dependency contains about 97 minutes of tool
runtime before staging and indexing. Accelerating just fast plot/report steps
will not remove this path.

### Basecalling: measure the replacement before resizing

This recorded run used **32** shards, rather than the current clinical default
of 16. All 32 completed in the initial session and were subsequently cached.
Tool runtimes range from 134.3 to 142.7 minutes (median 135.5); the span from
first submission to last completion is about 149 minutes. They requested
8 CPUs, 30 GiB and one GPU each, averaging about 1.5 CPU cores; maximum shard
RSS is 9.2 GiB. Shard balance looks good.

These reports use the old Dorado image. The replacement preprocessing digest
is not represented. Its OCI metadata declares the same Dorado 1.3.0 release;
any improvement from the replacement needs measurement on the actual GPU.
The user expects the replacement to partially address basecalling performance.
Keep it out of immediate resizing recommendations: compare reads/second,
GPU duty cycle/VRAM, CPU stalls and stage time on the replacement first. Do not
increase N or CPU/RAM based on these historical averages alone.

## Recommended experiments, in order

1. **DeepSomatic example generation: 8 → 16 shards.** Keep the existing
   32-CPU/100-GiB reservation and immutable CPU image for the first benchmark.
   This changes tool parallelism, not the Batch request. Compare make-examples,
   call-variants and postprocessing timings; aggregate task utilisation alone
   does not prove which phase dominates. Test 32 shards only if 16 helps.
   Keep the model, target BED and indel VAF settings fixed and compare outputs.

2. **Mapping on the current worker.** Benchmark 72 → 96 GiB memory and a
   balanced mapper/sort thread split before buying larger workers. The present
   peak is 59.3 GiB, with substantial disk traffic (~299 GiB read, 155 GiB
   written). Memory may reduce sort spill; it has not been shown to be the
   bottleneck. CPU averages below 32 do not prove CPU starvation, and asking
   both minimap2 and sort for almost all 32 threads can create contention.
   Change one factor at a time.

3. **Chromosome-parallel phasing/haplotagging.** Additional generic threads
   will not solve the approximately one-core phasing workload. A workflow
   extension could phase independent chromosomes concurrently and haplotag
   coordinate partitions, then merge and index once. Keep all contigs,
   off-target/unphased reads, original read groups and tags; check phase-set
   stability and BAM read counts. This requires implementation and regression
   validation, not just larger requests. One whole-sample BAM staged per
   chromosome would multiply storage/network work; prepare compact partitions
   once and benchmark a small bounded fan-out (e.g. four workers).

4. **DeepSomatic GPU inference.** A separate benchmark could use a verified,
   digest-pinned GPU build and one L4 GPU. The recorded task is CPU-only; the
   presence of a general GPU queue does not accelerate it. Example generation
   remains CPU work. This requires an image/tool-compatibility check, explicit
   accelerator request, routing, memory budgets and task metrics, followed by
   variant-output comparison. Avoid Blackwell for this caller until validated.

5. **Alignment partitioning if the single-worker experiment is insufficient.**
   Split reads into a small number of independent alignment batches, retain
   the current alignment preset and tags, then coordinate-merge/index one
   final BAM. This is separate from the POD5 shard feature: mapping currently
   consumes one consolidated FASTQ. Four workers could improve wall time, but
   add splitting, reference staging and merge cost. Validate deterministic read
   retention, supplementary alignments and caller outputs.

6. **Serial and I/O-heavy work.** ONTIME, NASVAR, Delly CNV, QDNAseq and Severus
   need stage profiling or tool-level changes before larger allocations.
   ONTIME's 0–1-hour filter is used by the existing high-coverage classifier
   path; it is not simply a time-series feature that can be disabled without
   changing classifier inputs. Preserve that sampling policy. Reading/filtering
   timestamps earlier in the read-processing path could avoid scanning the
   full aligned BAM, but requires careful sample and sort-equivalence checks.

## Queue, storage and instance options

Current Terraform uses exclusive `r8i.8xlarge` (32 CPUs / 256 GiB) high BAM
workers and `r8i.2xlarge` (8 CPUs / 64 GiB) low BAM workers, each with 2-TB gp3
scratch. Reservations intentionally occupy a whole worker. Low CPU utilisation
is **not** a reason to pack extra large-BAM tasks onto the same scratch disk.

Minimap2 spends about 17.1 minutes outside its tool runtime: approximately
3.3 minutes from submission to start and 13.8 minutes in the remaining wrapper
interval. Haplotag has about 7.4 minutes of nonqueue wrapper time; FASTQ
conversion and SeqKit have roughly 11 minutes. These intervals include staging,
startup, transfers, task monitoring and finalisation; they are not direct
measurements of S3 copy duration. Queue waits for the listed long tasks are
mostly below six minutes. Raising global fleet max_vCPUs is not the first
remedy for this single-specimen run.

The configured volume throughput is 500 MiB/s and 6000 IOPS. Increasing capacity
in GiB does not itself increase gp3 throughput. The current `r8i.8xlarge`
instance has EBS bandwidth headroom, but these reports do not establish volume
saturation. Measure host volume latency/queue depth, IOPS, throughput, CPU
I/O-wait, scratch free space, and transfer timings before purchasing more I/O.

AWS `DescribeInstanceTypes` was queried read-only in `us-west-2` on
21 November 2026. These are valid instance specifications; the query does not
guarantee capacity, quota or placement availability.

| Candidate                     | CPUs / RAM   | Storage/GPU                              | When useful                                                                                                               |
| ----------------------------- | ------------ | ---------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| Current `r8i.8xlarge`         | 32 / 256 GiB | EBS                                      | First mapping/DeepSomatic thread benchmarks fit on this worker.                                                           |
| `r8i.16xlarge`                | 64 / 512 GiB | EBS; 30-Gbps network                     | Single-task mapping benchmark if current-worker tuning helps but still leaves a compute bottleneck.                       |
| `g6.8xlarge`                  | 32 / 128 GiB | One L4, ~22.4-GiB VRAM; 900-GB local SSD | Matches the current 32-CPU caller reservation for a GPU DeepSomatic experiment.                                           |
| Existing general `g6.4xlarge` | 16 / 64 GiB  | One L4; 600-GB local SSD                 | Cannot satisfy the existing 32-CPU/100-GiB DeepSomatic request. Requires explicitly approved smaller requests.            |
| `i4i.8xlarge`                 | 32 / 256 GiB | 7.5-TB local NVMe; 18.75-Gbps network    | Storage experiment only if measurements show random I/O or scratch as the bottleneck; older CPU may be worse for compute. |

A 64-CPU mapping worker should use a dedicated queue with a full-worker
reservation, not silently raise the shared high-fleet reservation for all tools.
For an initial comparison, 64 CPUs / 96 GiB and at least the existing 2-TB
scratch would isolate CPU scaling. This proposal needs approval. Total sort
memory must remain bounded as thread counts change.

GPU workers need the same scratch-isolation discipline as CPU BAM workers.
The current general fleet lacks `g6.8xlarge`; adding it alone is insufficient
without explicit GPU routing. Its built-in local disk is smaller than the
existing 2-TB BAM scratch. Preserve a suitably sized task volume; mount/setup
local NVMe deliberately if using it. Merely selecting an instance with NVMe
does not move Nextflow or Docker scratch onto those disks.

## Approval boundary

The preferred first experiment is DeepSomatic with 16 tool shards on its
unchanged 32-CPU/100-GiB worker. The next is a mapping memory/thread benchmark
on the current instance. Ask before editing either configuration. Larger
mapping instances, GPU inference, parallel phasing and alignment partitioning
are separate decisions. No production infrastructure deployment or benchmark
analysis was launched during this evaluation.
