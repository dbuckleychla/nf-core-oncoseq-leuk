# Local review-code benchmark

Measured on 21 November 2026 using synthetic data only. No clinical files or
AWS jobs were used. Production lookup/consensus implementations and resource
requests were unchanged; optimized implementations are benchmark prototypes.

## Method

The benchmark opens an indexed BGZF VCF containing 30,000 invented records,
with 2-KiB annotation payloads and deliberate duplicate/multiallelic observations
(18 MiB compressed). Each measurement runs in a separate process. It compares
the current per-allele fetch with a bounded four-window exact-key cache at
128, 512 and 1024 bases. It also compares the existing all-pairs junction graph
with a partitioned, sorted-position graph builder, reusing the same event,
selection and ambiguity code.

Environment: macOS arm64, Python 3.14.6, pysam 0.23.3. These are single-run
measurements of already local files, not cold EBS or S3 staging benchmarks.
Both implementations use the same environment. Native Linux/amd64 production
performance and full-review wall time will differ.

## Indexed reference lookup

| Query pattern            | Implementation  | Seconds | Indexed fetches | Peak process RSS MiB | Relative speed |
| ------------------------ | --------------- | ------: | --------------: | -------------------: | -------------: |
| 20,000 clustered queries | Current         |  44.283 |          20,000 |                 31.1 |          1.00× |
| 20,000 clustered queries | 128-base cache  |   3.615 |           1,563 |                 32.2 |         12.25× |
| 20,000 clustered queries | 512-base cache  |   0.995 |             391 |                 32.8 |         44.52× |
| 20,000 clustered queries | 1024-base cache |   0.543 |             196 |                 33.6 |         81.53× |
| 300 sparse queries       | Current         |   0.667 |             300 |                 31.1 |          1.00× |
| 300 sparse queries       | 128-base cache  |   0.674 |             300 |                 31.8 |          0.99× |
| 300 sparse queries       | 512-base cache  |   0.716 |             300 |                 32.7 |          0.93× |
| 300 sparse queries       | 1024-base cache |   0.746 |             292 |                 33.5 |          0.89× |

All alternatives returned byte-identical logical record streams, including
20,021 matching records for the clustered pattern and 301 for the sparse
pattern. Duplicate observations remain duplicated; single-alt queries do not
match multiallelic records. The cache does not change exact-allele matching.

Caching has substantial potential where candidate positions cluster. It has
little benefit and up to about 12% elapsed-time overhead in this sparse case.
The tradeoff argues for a small/adaptive cache rather than blindly prefetching
large intervals. The measured RSS applies to one synthetic reference handle;
real annotation-row size and multiple catalogues require separate checks.

## Junction consensus

| Input junctions | Distribution               | Current seconds | Indexed seconds | Relative speed |
| --------------: | -------------------------- | --------------: | --------------: | -------------: |
|             600 | Separate clusters          |           0.201 |           0.005 |         38.63× |
|           1,800 | Separate clusters          |           1.781 |           0.015 |        117.76× |
|           3,600 | Separate clusters          |           7.262 |           0.032 |        223.79× |
|             180 | Dense colliding candidates |           0.025 |           0.004 |          5.54× |

Complete event JSON digests match in every case, including stable identifiers,
selection, caller evidence and ambiguity decisions. The data include missing
orientation and transitive-chain cases. Indexing removes irrelevant comparisons;
truly dense compatible regions can still require quadratic work and large
ambiguity output. It does not relax the compatibility rules.

## Implications

These results support contained changes inside the review code, without file
splitting, separate parallel tasks or larger resource requests. They do not
predict an 80× improvement for the whole REVIEW_COMPLETE process: catalogue
staging, VEP-result processing, SQLite writes, export compression and evidence
serialization remain. Progress logging is needed to measure those phases.

The reported missing STRANDS-header crash is a separate correctness defect;
its regression fix does not incorporate these optimization prototypes.

## Reproduce

Use Python with pysam 0.23.3 installed, or the pinned review runtime:

```bash
docker run --rm --platform linux/amd64 \
    -v "$PWD:/workspace" -w /workspace \
    community.wave.seqera.io/library/oncoseq_review_1_0_0@sha256:b143e350e1dae33ad26758fc6f9ccf2392cb8b49881de73960f42525602d5eb5 \
    python3 tests/review/benchmark_review.py --output /workspace/review-benchmark.json
```

On Apple Silicon, that Docker command uses amd64 emulation; the measurements
above used native arm64 Python. Synthetic VCFs are created in temporary storage
and removed afterward. The output contains counts, timings and logical digests,
not patient records. Raw measurements are in [review-local.json](review-local.json).
