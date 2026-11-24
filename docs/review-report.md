# Optional molecular review report

The `clinical` profile enables the molecular review by default. For other
profiles, enable `--report_review` to produce `reports/<sample>.review.html` alongside the
existing report. The dark-themed Oncoseq review report, aggregated by FrankONTstein, embeds its
data, scripts, styles and displayed
figures. It works offline and requires JavaScript for interactive tables.

The report prioritises target variants, presents caller evidence separately,
consolidates compatible fusion breakpoints and provides coverage and classifier
views. It preserves native results. Caller agreement and COSMIC membership do
not establish somatic origin or pathogenicity.

## Run on AWS Batch

Use the existing institutional preset and Terraform-generated configuration.
Supply the input, output and work paths used for your run:

```bash
nextflow run . \
    -profile clinical,aws \
    -c aws.generated.config \
    -params-file assets/institutional_hg38.params.json \
    --input validation-pod5.csv \
    --outdir "$OUT_URI" \
    -work-dir "$WORK_URI" \
    -resume
```

Use the actual filename of your Terraform-generated configuration in place of
`aws.generated.config`. The review manifest accepts local files and S3 URIs.
Clinical defaults use the institutional `hg38/aux/` prefix under
`s3://plm-cpm-sandbox-research-usw2/ONT/Oncoseq_refdata/` and the bundled
`assets/institutional_hg38.review.json` development mapping. No uploaded
manifest or operator-generated checksums are required during development. Explicit
`--reference_aux_cache` and `--review_references` values override those defaults.
Use `--report_review false` to disable the extra report.

Batch stages the resources using the job role. Reference files and indexes are checked before
basecalling. Supplied checksums are validated, but checksums may be omitted in
the explicitly marked development mapping. Strict manifests also verify the
analytical FASTA checksum. Pin and validate reference identity before launch
(`oncoseq-u7o`); reports identify unpinned development inputs.

VEP remains at release 115. Supply the existing matching cache. The report
does not download annotation resources or models at runtime.

For copy/upload commands using the institutional reference root, see
[Stage institutional references](stage-review-references.md).

## Reference manifest

Required resources are `fasta`, `fai`, `targets`, `gff`, `gnomad_genome`,
`gnomad_exome`, `clinvar`, `clinvar_common_pathogenic`, `cosmic_targeted`
and `cosmic_evidence`. The targeted VCF and prepared SQLite evidence database
are mandatory, including in development. Missing resources, different release
labels or a database prepared against different VCF bytes fail validation.
Other `cosmic_*` products are not accepted by this targeted-screen workflow.

Use the authoritative unpadded assay BED for `targets`. Coverage shown in the
report retains the existing workflow's native intervals. Reporting genes and
preferred transcripts are controlled separately by `--review_panel`.

The clinical development mapping records reference paths and source labels
without requiring SHA-256 checksums. Use the following strict manifest workflow
when pinning references before launch; it is optional during development.
Each pinned reference has an assembly-compatible source, release and SHA-256. VCFs
require a matching adjacent `.tbi` or `.csi` index. Use the exact FASTA and its
adjacent `<fasta>.fai`. Annotation chromosome names must match the FASTA.

Create a local resource mapping after preparing the files:

```json
{
  "fasta": { "path": "/ref/hg38.analysisSet.fa", "release": "GRCh38" },
  "fai": { "path": "/ref/hg38.analysisSet.fa.fai", "release": "GRCh38" },
  "targets": { "path": "/ref/targets.unpadded.bed", "release": "assay-panel-v1" },
  "gff": { "path": "/ref/GRCh38.gff", "release": "GRCh38.p14" },
  "gnomad_genome": { "path": "/ref/gnomad.genomes.vcf.gz", "release": "4.0" },
  "gnomad_exome": { "path": "/ref/gnomad.exomes.vcf.gz", "release": "4.0" },
  "clinvar": { "path": "/ref/clinvar.vcf.gz", "release": "YOUR_FIXED_RELEASE" },
  "clinvar_common_pathogenic": { "path": "/ref/rescue.vcf.gz", "release": "assay-rescue-v1" },
  "cosmic_targeted": { "path": "/ref/cosmic.targeted.v105.GRCh38.vcf.gz", "release": "v105" },
  "cosmic_evidence": { "path": "/ref/cosmic.targeted.v105.GRCh38.evidence.sqlite", "release": "v105" }
}
```

Save this mapping as `review.resource-files.json`. Generate the checksums and
index entries with the pinned runtime:

```bash
export REVIEW_IMAGE='community.wave.seqera.io/library/oncoseq_review_1_0_0@sha256:b143e350e1dae33ad26758fc6f9ccf2392cb8b49881de73960f42525602d5eb5'

docker run --rm --platform linux/amd64 \
    -v "$PWD:$PWD" -v /ref:/ref -w "$PWD" \
    --entrypoint python3 "$REVIEW_IMAGE" \
    bin/review_report.py manifest \
    --resources review.resource-files.json \
    --output review.references.json
```

For S3 deployment, add each file's existing S3 `uri` to the mapping before
generating the manifest. The helper reads the local `path` to calculate the
checksum and records the supplied deployment `uri`. Index URIs use the same
file URI plus `.tbi` or `.csi`. Files with relative URIs resolve beside the
manifest. Keep resource filenames unique within the manifest.

The helper records operator-supplied releases; verify each release against its
source metadata. Pin the generated manifest with the run. Do not replace its
resources in place with new releases.

## Prepare COSMIC

The v105 GRCh38 targeted-screen normalised VCF was inspected successfully.
It uses numeric chromosome names and has transcript rows repeated per allele.
Prepare a separate indexed derivative against the full production FASTA:

```bash
docker run --rm --platform linux/amd64 \
    -v "$PWD:$PWD" -v /ref:/ref -v /data:/data -w "$PWD" \
    --entrypoint python3 "$REVIEW_IMAGE" \
    bin/review_report.py prepare-vcf \
    --source /data/Cosmic_CompleteTargetedScreensMutant_Normal_v105_GRCh38.vcf.gz \
    --fasta /ref/hg38.analysisSet.fa \
    --output /ref/cosmic.targeted.v105.GRCh38.vcf.gz \
    --skip-masked-reference
```

The helper records chromosome renaming, checks reference alleles, normalises,
sorts and indexes the derivative. It writes a preparation audit with source,
FASTA, output and index checksums. Keep the original downloaded file.

The analysis-set FASTA masks some reference regions with `N` bases (including
parts of chrY). `--skip-masked-reference` archives those original VCF records
in `.masked-reference.vcf.gz` and excludes them from the derivative, recording
counts by contig and checksums in the preparation audit. It still fails on REF
mismatches at unmasked bases. Default preparation remains strict without this
flag. A source index is unnecessary because preparation reads sequentially.
The inspected v105 targeted-screen catalogue completed against the full
institutional analysis-set FASTA after excluding four chrY transcript records;
all remaining REF alleles passed strict normalization.

COSMIC records with FILTER `.` are valid annotation input. Recurrence values
remain associated with their catalogue and COSV identifier. Transcript rows
are not added together. Inconsistent values are displayed separately. Gene
Census tier is omitted from report displays. These VCF counts are not
tissue-specific; the evidence database supplies tissue-specific counts separately.

Extract the matching GRCh38 release's `Cosmic_CompleteTargetedScreensMutant_Tsv`
and `Cosmic_Classification_Tsv` archives, retaining the compressed TSVs. Prepare
an indexed evidence database once per reference release (outside the repository):

```bash
python3 bin/review_report.py prepare-cosmic \
    --mutations /ref/Cosmic_CompleteTargetedScreensMutant_v105_GRCh38.tsv.gz \
    --classification /ref/Cosmic_Classification_v105_GRCh38.tsv.gz \
    --vcf /ref/cosmic.targeted.v105.GRCh38.vcf.gz \
    --release v105 \
    --output /ref/cosmic.targeted.v105.GRCh38.evidence.sqlite
```

The supplied v105 files were audited across 1,350,464 VCF transcript records
representing 360,801 distinct COSV IDs: every ID joined to positive evidence,
and every retained phenotype ID resolved. Of 63,485,816 mutation TSV rows,
59,683,220 were negative screens; 158,277 positive rows lacked a COSV ID and
could not join. After deduplication the database contained 3,641,814 evidence
rows. These are reference-catalogue counts, not sample results.

Use Python with the review runtime dependencies, or the pinned container above.
The helper streams the mutation TSV into SQLite; annotation uses indexed COSV
lookups instead of rescanning the TSV for each sample. It records source
checksums, release, assembly, row counts and missing-classification counts inside
the database. Existing outputs are never overwritten. Keep licensed source data
and the production database in reference storage, outside this public repository.

The joins are normalized GRCh38 CHROM/POS/REF/ALT → VCF COSV ID →
`GENOMIC_MUTATION_ID` → `COSMIC_PHENOTYPE_ID`. Only `POSITIVE_SCREEN=y`
rows enter the evidence database. Negative screens and positive records without
a COSV ID cannot establish an allele match. The TSV's HGVS positions are not
used as alternative coordinate matches. Transcript duplicates are retained as
annotations but distinct `COSMIC_SAMPLE_ID` values determine recurrence. Counts
across phenotype groups can overlap; they are not summed into unique totals.
`SAMPLE_NAME` is never imported. Full positive evidence retains COSMIC reference
sample IDs, transcript/HGVS, phenotype IDs, publications and reported somatic status.
Missing TSV or classification joins remain explicit warnings; they never become
blood evidence. Conflicting duplicate phenotype definitions fail preparation.

Policy version 4 puts matches classified as
`PRIMARY_SITE=haematopoietic_and_lymphoid_tissue` in the primary table, including
lymphoid disease. A separate high-priority table retains variants without COSMIC matches that meet
the earlier reviewed ClinVar P/LP or curated predicted loss-of-function criteria.
The collapsed expanded table includes other exact COSMIC matches and remaining
ClinVar/VEP candidates. Tables order by distinct COSMIC PubMed publications,
highest first. Caller PASS and population eligibility
(including the explicit whitelist) apply to both tables and consensus selection.
ClinVar conflicts and benign assertions remain visible even for primary COSMIC
matches. COSMIC transcript/HGVS supplies the representative matched annotation;
VEP transcript consequences remain separate evidence. COSMIC's reported somatic
status describes its catalogue samples, not the current sample.

Annotated VCFs add `ONCOSEQ_REVIEW_COSMIC` and
`ONCOSEQ_REVIEW_COSMIC_BLOOD`. JSON and complete annotation TSV exports retain the
joined evidence, phenotype classifications and join warnings. Native caller FILTER
and read evidence remain unchanged. Deploy the prepared database to the URI in
`cosmic_evidence` before starting a run with the updated reference mapping.

## Output retention

| Output                            | Default       |
| --------------------------------- | ------------- |
| Original aligned BAM and index    | Published     |
| Haplotagged BAM and index         | Published     |
| Consolidated all-read sample uBAM | Published     |
| Individual shard BAMs             | Not published |
| Generated FASTQ                   | Not published |
| Diagnostic intermediate BAMs      | Not published |

Use `--publish_fastq` to publish generated FASTQ. Use `--publish_ubam false` to
omit the uBAM archive. `--publish_intermediate_bams` enables diagnostic BAM
publication without exposing shards.

Simplex uBAM archives contain the consolidated basecalled reads before quality
splitting. Multiplex archives contain trimmed, assigned sample reads before
quality splitting, including quality-failing reads. They exclude unclassified
barcode reads. Both aligned BAMs remain pending guidance on archival retention.

Publication flags do not delete existing outputs or work files. Use a fresh
output prefix when checking the new output inventory. Preserve the work path
and `-resume` when reusing analytical tasks.

## Inspect and re-render

Machine-readable data and complete annotation exports are published beneath
`reports/review/<sample>/`. Original caller outputs remain available. The new
exports contain valid target alleles and original caller records; primary
report selection additionally applies caller PASS and review policy.

Copy the sample's review directory locally. Re-render its prepared JSON and
figures without running callers:

```bash
docker run --rm --platform linux/amd64 \
    -v "$PWD:$PWD" -v /data:/data -w "$PWD" \
    --entrypoint python3 "$REVIEW_IMAGE" \
    bin/review_report.py render \
    --input /data/review/SAMPLE/review.json \
    --figures /data/review/SAMPLE/figures \
    --output /data/SAMPLE.review.html
```

The HTML works independently. Links to full exports need the adjacent results
directory. When copying a results tree, keep the report's `review/<sample>`
directory beside the HTML.

For an existing completed run, inventory its native files with
`review_report.py import --results /data/run --sample SAMPLE --output inputs.json`.
Missing required native files produce an explicit error. The `prepare` and
`complete` commands accept this inventory and a local reference manifest;
`complete` requires VEP 115 annotation of the prepared `catalogue/catalog.vcf`.
Using the pipeline's `--report_review -resume` path performs these steps.

## Selection and verification

Primary small-variant review requires exact PASS. Population AF greater than
0.01 is excluded unless the exact allele is whitelisted. Valid gnomAD genome
and exome AC/AN are pooled under the recorded policy, with both datasets kept
in the evidence. ClinVar assertion conflicts and review status remain visible.
Predicted loss of function, additional panel coding review and supporting
evidence have separate views. Allele fractions and read counts are never
averaged or added across callers.

Fusion consensus requires compatible event type, chromosome pair and known
orientation, with both breakpoints within `--review_fusion_distance` bases
(default 100). Ambiguous matches and transitive chains remain separate. Native
sample aliases are preserved and mapped explicitly. Selected fusion candidates
require at least three supporting reads and applicable native filters.

The implementation is tested with synthetic references and native record
fixtures. Production acceptance requires the staged reference manifest and
representative runs. Do not use the synthetic `test_data/review` references for
analysis. The design and evaluation are in [report-redesign-plan.md](report-redesign-plan.md).

## Completion performance and progress

Review completion uses bounded adaptive exact-allele reference caches and
indexed structural-junction matching. Unknown orientation remains unknown,
including when a VCF omits the optional STRANDS header. Scientific selection,
native evidence and output contracts are preserved. No files are split and no
additional parallel jobs or tool threads are introduced.

Count-only logs identify annotation, export and structural stages; annotation
and export counters appear every 10,000 records. Cache statistics report fetch
counts, hits, oversized-window fallbacks and accounted memory. The cache budget
is 16 MiB per reference, and existing CPU/RAM requests remain unchanged.
[Local benchmarks](benchmarks/review-local.md) show substantial clustered-case
gains while documenting sparse-case overhead and limits on extrapolation.

## Resume and input ordering

Review artifacts are ordered canonically by kind and source URI before task
submission. This prevents asynchronous upstream arrival order from changing
the preparation manifest and its task hash on resume. Adopting the ordering
fix can invalidate a previous noncanonical review preparation once; subsequent
runs with unchanged inputs and parameters should reuse it. Downstream review
tasks rerun when their preparation output location changes.

The regression check reverses synthetic artifact arrivals and confirms that
Nextflow reuses the same REVIEW_PREPARE hash with `-resume`. This ordering fix
does not change resource requests or add parallel jobs. Let an already-running
analysis finish rather than cancelling it solely to install this fix.
