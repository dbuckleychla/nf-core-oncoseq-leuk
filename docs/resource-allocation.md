# Clinical resources and failure handling

This audit compares Oncoseq release `v1.0.2` with the current FrankONTstein resource configuration, including its local QDNAseq memory override. The settings are production candidates. Representative Batch runs must establish their memory, time and disk requirements.

## Worker placement

The dedicated BAM high fleet contains only `r8i.8xlarge` workers, with 32 vCPUs and 256 GiB RAM. Every routed task reserves all 32 vCPUs. The low fleet contains only `r8i.2xlarge` workers, with eight vCPUs and 64 GiB RAM. Every routed task reserves all eight. Low-tier tasks request the high queue when their actual memory exceeds 48 GB, including on retries. A high-only deployment preserves isolation when the low queue is absent.

Clinical AWS execution requires a high queue, `max_cpus >= 32`, exact high/low reservations of 32/8, and `max_memory <= 224 GB`. CPU reservations control Batch placement. An aggregate fleet capacity such as 256 or 512 vCPUs allows more workers; it does not change per-worker exclusivity. These rules assume the queues are those created by this Terraform stack, with no shared general CPU compute environment.

High routing covers alignment, merge, sort, phasing/haplotagging, BAM-based Classy modules, the three small-variant callers, Severus and VEP. VEP also stages a large reference-cache directory, so it uses the exclusive high worker. Low routing covers demux/trim, BAM index/view/filter/count/FASTQ conversion, Sniffles, Stellerator, NASVAR, CNA tools, mosdepth/cramino QC, Ontime, read preprocessing/statistics and BAM-based Figeno plots. An automated input/alias audit prevents BAM and read modules from silently bypassing these routes.

Dorado uses the separate basecalling queue with one GPU per job. When `basecall_shards > 1`, AWS requests 8 CPUs and clinical execution requests a fixed 30 GB per shard, allowing two jobs to share a two-GPU worker. Single-shard runs retain 32 CPUs and 100 GB, or 180 GB in WGS mode. Sharded memory remains 30 GB on retries; out-of-memory failures stop rather than retry at the same limit. General GPU routing remains separate. Blackwell inference compatibility remains unverified.

Dorado requests one GPU by default, while the example g7e/g6 workers contain two/four GPUs. After SUP model acceptance, benchmark a two-GPU process override with `device = cuda:all` to assess throughput and memory. Both instance types can satisfy that request; the worker CPU reservation still prevents concurrent staging. A capacity fallback to g6 does not automatically move a failed inference task away from g7e.

## Step budgets

The table gives initial clinical memory. Memory grows linearly with attempt, retains the threefold `huge` multiplier, and is capped by `max_memory`; Dorado also has a 180-GB ceiling for g6 host headroom. Other steps retain Oncoseq's existing labels and resource formulas.

| Step               | Initial memory                              | Tool parallelism                               | Reason                                                                                    |
| ------------------ | ------------------------------------------- | ---------------------------------------------- | ----------------------------------------------------------------------------------------- |
| Minimap2 alignment | 72 GB                                       | 32 reserved CPUs on Batch                      | Retain Oncoseq's larger budget; cap sort memory per thread and leave half for alignment.  |
| VEP                | 128 GB                                      | 32 forks from the upstream wrapper             | Isolate reference-cache staging and increase memory for the full-worker fork count.       |
| SAMtools sort      | 64 GB                                       | Reserved task CPUs                             | Increase the former 36-GB budget for sorting and cap `-m` per thread.                     |
| SAMtools merge     | 32 GB                                       | Reserved task CPUs                             | Increase the former 12-GB budget; staging and output copies remain the main disk concern. |
| WhatsHap phase     | 32 GB                                       | Existing wrapper behaviour                     | Increase the former 12-GB budget; phase and haplotag use exclusive high workers.          |
| Classy combined    | 64 GB                                       | 8 threads                                      | Match FrankONTstein's larger budget instead of the former 20 GB.                          |
| Severus            | 64 GB                                       | 4 threads                                      | Match the legacy increase after a large ONT input exceeded its parsing budget.            |
| Clair3             | 100 GB                                      | 16 threads                                     | Retain Oncoseq's larger memory allocation while separating threads from placement.        |
| ClairS-TO          | 180 GB                                      | 16 threads                                     | Retain Oncoseq's larger allocation; a memory retry reaches at most 224 GB.                |
| DeepSomatic        | 100 GB                                      | 8 shards                                       | Retain Oncoseq's larger memory budget and FrankONTstein's bounded shard count.            |
| Stellerator        | 12 GB                                       | 4 threads                                      | Raise the former 10 GB to the legacy budget.                                              |
| Dorado             | 30 GB sharded; 100 GB unsharded, 180 GB WGS | 8 CPUs sharded, 32 unsharded on Batch; one GPU | Permit concurrent shards while retaining larger single-shard budgets.                     |

The expensive alignment, sorting/merging, phasing, classification, caller and QC steps receive the full `max_time` on every attempt, default 24 hours. These are allocation choices rather than measured requirements. Tool thread caps use `task.ext.threads` in six repository-owned process interfaces under `modules/local/clinical_threads`. They derive from pinned Oncoseq submodules, retain their inputs/outputs and image pins, and record source revisions in `upstream.json`. The original gitlinks remain intact. Imported nf-core modules retain their upstream interfaces.

## Failure policy

Clinical jobs never ignore failures. Exit 137 can retry with increased memory; all steps stop immediately when their memory ceiling is already reached. Exit 143 permits a bounded retry for interruption. Other exits, including tool/argument errors, segmentation faults and typical disk-full failures returning exit 1, finish the workflow without repeating the same command. Most steps allow one retry; the expensive group allows two. The pipeline records exit status, attempts, memory, queue and I/O counters in the trace.

Retries do not enlarge EBS disks. CPU isolation prevents simultaneous BAM tasks from sharing a worker, but one task's inputs, temporary files, outputs and Docker layers must still fit its disk. High and low fleets default to 2000 GiB each; the GPU/general fleets use `scratch_gb`. Diagnose a single-task disk failure from CloudWatch disk/inode metrics and task logs, then increase the appropriate disk or reduce input/chunk size. Submission, image-pull and S3 staging failures may occur before tool execution and still require executor/host diagnostics.

Acceptance requires concurrent samples with verified exclusive Batch placement, the largest expected merge/sort input, measured disk/RAM peaks, and a successful resumed run. The automated routing and retry tests do not establish clinical performance or result equivalence.

`FIGENO_PAN_CHR` runs at most four chromosome plotting jobs, capped by task CPUs. Clinical execution requests 24 GB initially and 48 GB on its single retry. AWS reserves an isolated low worker (eight CPUs), or a high worker when no low queue exists. Exit 11 is eligible for this bounded retry because AWS reported the observed plotting out-of-memory kill with that status; exit 137 is also eligible below the memory ceiling. Other plot failures remain fatal. This replaces the script's default 20 concurrent jobs under a 6-GB task budget.
