# Optional review report and output retention

## Purpose and agreed behaviour

Add an optional, self-contained HTML report for leukemia review. The existing
report remains available. The new report uses native analytical outputs and
explicit selection rules. It retains useful fusion, chromosome, coverage and
classification figures in a readable layout.

Retain the original aligned BAM and haplotagged BAM with their indexes. Retain
one consolidated per-sample uBAM containing quality-passing and quality-failing
reads by default. Stop publishing shard BAMs. Make FASTQ and diagnostic
intermediate BAM publication explicit options. Publication changes do not
remove work files or existing results.

The new report does not change caller thresholds, basecalling, multiplexing,
trimming, mapping, phasing or classification. VEP stays at release 115. The
existing 0.05 indel candidate thresholds remain in effect.

## Evaluation

The examined execution contains about 231.5 GiB of uncompressed FASTQ, 89 GiB
of shard BAMs and two aligned BAMs of about 78.2 GiB each. These large files
were listed in the output inventory; their read contents were not compared.

The current HTML is about 111 MB. It contains 58 images and 12 tables. Section
and chromosome ordering reflect asynchronous task completion. Fusion images
reach approximately 19,000 pixels wide and chromosome images approximately
28,000 pixels. Scaling these images to page width makes their labels unreadable.

The small-variant tables mix PASS, LowQual, NonSomatic and GERMLINE. Fixed VEP
CSQ column positions make the current presentation parser sensitive to field
order. Structural tables can confuse consequence with event type. The
Stellerator presentation contains single-read records despite a three-read
caption because repeated rows can satisfy a row-count condition.

Coverage joins must use interval identity rather than gene label alone.
Coverage z-scores are QC measurements, not copy-number diagnoses. Mean region
depth cannot establish the percentage of bases above a depth threshold.

The current WhatsHap `--skip-missing-contigs` option can omit reads on contigs
absent from the phasing VCF. Keep both aligned BAMs pending author guidance.
The follow-up must address unmapped and supplementary records, MM/ML/MN and
move tags, and a comparison establishing whether one BAM can be retired.

### COSMIC audit

The downloaded `Cosmic_CompleteTargetedScreensMutant_Normal_v105_GRCh38.vcf.gz`
is an eight-column, BGZF annotation VCF. Its archive checksum matches COSMIC's
published v105 checksum. Complete parsing with bcftools succeeded.

| Property                                          | Result                        |
| ------------------------------------------------- | ----------------------------- |
| Transcript rows                                   | 1,350,468                     |
| Distinct genomic alleles                          | 360,791                       |
| Distinct COSV identifiers                         | 360,803                       |
| Chromosomes                                       | 1-22, X and Y                 |
| Leukemia panel genes represented                  | 74 of 74                      |
| Run-derived target regions with catalogue alleles | 1,073 of 1,134                |
| Review-panel genes with spatial target overlap    | 74 of 74                      |
| Chromosome 17 REF validation                      | 144,982 rows, zero mismatches |
| Supplied index                                    | None                          |

The spatial audit used the run's generated no-padding BED. Full production
validation requires the authoritative unpadded assay BED and exact FASTA.
Chromosome 17 validation used the available GRCh38 primary chromosome fixture.

Confirmed catalogue examples include FLT3 D835Y, NPM1 W288Cfs\*12, IDH1 R132H,
IDH2 R140Q, NRAS Q61R, KRAS G12D, JAK2 V617F and BRAF V600E.

The VCF provides COSV and legacy IDs, gene, transcript, HGVS, canonical status,
sequence change, targeted-screen recurrence and Gene Census tier. All records
have FILTER `.`. This is valid catalogue input and must not be subjected to
the caller PASS rule. `TIER` is a gene-level Cancer Gene Census tier.

There are many transcript rows per allele, twelve alleles with multiple COSV
IDs and seven COSV IDs with conflicting recurrence counts. Preserve all IDs,
deduplicate identical observations and flag inconsistent counts. Do not sum
counts across transcripts, IDs or products. Counts are not leukemia-specific.
Catalogue absence does not establish benignity or sequencing coverage. This
coding resource does not replace fusion or noncoding analysis.

Source VCF SHA-256:
`a82e7cbe9beb652a21e73fa5acbc65d73b96bfa9c1da3b91f4361b9dcb73822c`.

Source archive SHA-256:
`b4c9e91346b5929f80a50d82c4d6f0ca8d9d3f9e51e01a30a1dab658c944e4ca`.

## Interfaces

| Parameter                   | Default                | Behaviour                                      |
| --------------------------- | ---------------------- | ---------------------------------------------- |
| `report_review`             | false                  | Add the offline review HTML                    |
| `review_references`         | unset                  | Versioned JSON resource manifest               |
| `review_panel`              | bundled leukemia panel | Gene and preferred-transcript TSV              |
| `review_fusion_distance`    | 100                    | Maximum bases at each breakpoint for consensus |
| `publish_fastq`             | false                  | Publish generated FASTQ                        |
| `publish_ubam`              | true                   | Publish consolidated all-read sample uBAM      |
| `publish_intermediate_bams` | false                  | Publish diagnostic BAMs, excluding shards      |

Declare parameters in the typed entrypoint, parameter config and JSON schema.
Require valid review references when the report is enabled. An absent optional
resource is displayed as not configured. An explicitly supplied invalid
resource fails preflight.

Publish `reports/<sample>.review.html` and machine-readable review data beneath
`reports/review/<sample>/`. Include manifest, complete small-variant and
structural exports, caller evidence, report JSON and figure assets. The HTML
embeds its data, styles, scripts and displayed figures and works offline.

Provide standalone import, annotation and render commands for completed runs.
Rendering prepared review data must not launch analytical tools.

## Publication implementation

Use the existing publishDir mechanism with process-specific save rules.
Suppress FASTQ in both dedicated reads publishers and the global publisher.
Cover SAMTOOLS_TOFASTQ, SAMTOOLS_TIME_FASTQ, CHOPPER and concatenation producers.

Publish multi-shard simplex uBAM from the existing merge producer. With one
actual shard, publish its completed output under the sample archive name.
Publish multiplex sample uBAM from DORADO_TRIM before quality splitting.
Multiplex archives contain assigned sample reads; unclassified barcode reads
are outside those archives. Prevent the demultiplexer directory publisher
from exposing intermediate BAMs. Avoid an extra full-uBAM staging task.

Preserve explicit alignment and phasing publishers and their indexes. Suppress
other BAM publication by default. Diagnostic opt-in never enables shard
publication. Verify publication after resume and with multiple samples.

## Annotation resources

Use the exact production FASTA/FAI, authoritative unpadded targets, reporting
panel, VEP 115 cache, gnomAD genome and exome VCFs, ClinVar, exact common-pathogenic
rescue VCF, GRCh38 GFF and curated fusion pairs. Reuse validated institutional
resources. Keep enrichment padding separate from assay targets.

The manifest records assembly, release, source URI, checksum and index for each
resource. Prepare large resources once: explicit chromosome renaming, allele
splitting, reference normalisation, sorting and indexing. Reject REF mismatches.
Record source and prepared checksums. Do not download mutable releases at runtime.

Preferred reporting transcripts may differ from the VEP cache release. Record
missing preferences and the matching-gene fallback. Do not upgrade the cache
or image independently.

Use the targeted-screen COSMIC v105 GRCh38 normalised VCF initially. Permit the
matching genome-screen VCF as an additional product. Keep source counts
distinct. Tissue-specific counts require separate TSV data and sample
deduplication; they are outside the first implementation.

## Small variants

Normalise native caller VCFs against the exact FASTA and select target overlap.
Build a disk-backed unique-allele catalogue. Annotate candidate union once with
VEP 115, then perform exact gnomAD, ClinVar, rescue and optional COSMIC joins.
Parse CSQ fields from the VCF header. Preserve native FILTER, genotype, depth,
allele depth, AF and identifiers, including failed counterparts to selected
alleles. Keep full exports separate from concise presentation.

Primary review requires exact caller PASS. Exclude population AF greater than
0.01 unless the exact allele is whitelisted. Preserve valid exome/genome AC/AN
and the existing pooled calculation. Missing frequency remains unknown.
COSMIC presence does not override caller filters or population exclusion.

Separate reviewed ClinVar P/LP, predicted loss of function, additional panel
coding review and supporting evidence. Retain conflicting/uncertain/benign
assertions and review status. Predicted loss of function is not automatically
pathogenic. Explain each selection. Do not introduce a new universal depth or
AF exclusion in the report.

Select preferred panel transcript, MANE Select, MANE Plus Clinical and canonical
protein-coding context deterministically. Preserve alternative transcript and
overlapping-gene context. Compare ClairS-TO and DeepSomatic by exact normalised
PASS allele. Keep Clair3 evidence independent. Never average AF or sum caller
support. Apply scientific selection before pagination; no hidden row cap.

## Structural events and fusions

Adapt native Sniffles, Severus, Stellerator and NASVAR records. Preserve sample,
caller, native IDs, filters, both coordinates, orientation, event type, mate
links, native support/depth/AF, confidence intervals and breakpoint context.
Maintain native coordinates alongside a zero-based interbase representation.
Test all four BND orientations and NASVAR coordinate conversion.

Validate reciprocal mate links before combining records within a caller. Use
the conservative minimum mate support for eligibility while retaining both
values. Do not count repeated transcript/table rows as supporting reads.

Across callers, require matching chromosome pair, compatible event type and
known compatible orientation. Require both breakpoints within the configured
distance. Swap complete endpoint context during canonical ordering. Keep DEL,
DUP and INV distinct from BND/TRA junctions. Unknown orientation does not
automatically merge. Gene-pair agreement alone is insufficient.

Merge only unambiguous pairwise-compatible groups with at most one junction
per caller. Reject transitive chains and leave ambiguous groups linked but
separate. Choose a deterministic native representative and retain all native
coordinates. Curated tool-specific fusion tolerances do not redefine consensus.

Prioritised fusion candidates require at least three supporting reads and the
applicable native PASS state, including both explicit mates. Keep lower-support
and failed records in supporting review. NASVAR retains its native semantics.
Label gene-pair structural evidence accurately; do not infer expression or
reading frame from genomic overlap.

## Report and figures

Define versioned report data for sample, reference, policy, caller status,
classification, findings, native NASVAR, coverage/QC, figures and provenance.
Assemble by explicit sample/artifact type. Distinguish completed-empty, not run,
not configured and failed. Missing required output cannot become empty success.
Stable event IDs prevent gene-pair figure collisions.

Use Python, Jinja2, embedded JavaScript and CSS without Quarto. Pin the released
runtime image and package versions. Implement keyboard navigation, search,
sorting, pagination, evidence expansion and print layout without network access.

Navigation order: overview, classification, fusions, small variants, CNV/other
SVs, NASVAR, coverage/QC, target loci, methods/files. Blood classifiers appear
first. Use model-native cutoffs and score semantics. Show Optic disabled.

Small-variant rows show gene, change, reason, caller agreement, native AF/depth
and annotation. Fusion rows show pair, event type, chromosome pair and separate
caller support. Detailed evidence remains accessible. Rename important genes
to target-locus review and explain selection and assay limits.

Regenerate excessive-width figures with readable geometry, labels and legends.
Use separate breakpoint panels when useful. Target 2,000-3,000 pixel raster
width, preserve original outputs and provide full-screen zoom/pan. Deduplicate
embedded images. Chromosomes follow 1-22, X, Y, mitochondrial then reference
order. Correct bp/kb units. Replace the wide methylation SVG as the primary view
with cards and compact scores. Coverage joins use interval identity and keep
alignment filters distinct. Target HTML below 30 MiB for the examined sample;
warn above 50 MiB without silently dropping findings.

## Compute and delivery

Initial annotation and SV tasks use two tool threads and 8 GiB, rendering two
to four threads and 8 GiB, VEP eight threads and 32 GiB, BAM-backed figures four
threads and 24 GiB. Use existing queues and whole-worker reservations for large
staging. One bounded memory retry reaches 16 GiB or 48 GiB as applicable.
Do not retry deterministic input failures. No new GPU or Terraform resources.
Cache by reference/panel checksums, policy, code and parameters.

Work packages cover interfaces/fixtures, publication, small variants, structural
evidence, QC/classifier adapters, figures, HTML and integration. The integration
owner controls shared schema, parameter, image and workflow changes. Track
progress and follow-ups in Beads. Code provenance belongs with owned adapters;
the report has no runtime dependency on an internal comparison checkout.

## Acceptance

Focused synthetic tests cover exact alleles, normalisation, multiallelic data,
FILTER states, AF boundaries, whitelist exactness, ClinVar conflicts, transcript
fallbacks, COSMIC duplicates/conflicts, BND orientations/mates, chain/ambiguous
consensus and native per-caller evidence. Do not commit patient fixtures.

Workflow tests cover one/multiple shards, simplex/multiplex naming, no published
shards, all-read uBAM, both aligned BAMs/indexes, every FASTQ producer, diagnostic
opt-in, report disabled/enabled, typed parameters and resume. Run the affected
nf-test suite, required lint/hooks and a minimal exercising profile. Avoid
full basecalling/calling reruns for report layout.

Render supplied execution outputs in a fresh evaluation directory. Verify
offline use, readable 1280/1440-pixel layouts, no page-wide horizontal scroll,
keyboard/table/print behaviour, chromosome order, truthful units/captions,
conservative fusion consolidation, complete selected findings and traceability.
Record versions, digests, HTML size, duration and peak memory. Commit and push
authorised changes to the fork after required checks.
