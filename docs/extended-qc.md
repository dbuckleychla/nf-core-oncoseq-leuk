# Extended sequencing QC

Cellular-DNA runs now add regional read-length, coverage, and independent methylation QC. The new metrics appear in both the native report and the optional review report's QC tab. Existing analytical steps and outputs are retained.

Coverage breadth includes ≥1×, ≥5×, ≥10×, ≥20×, ≥30×, ≥50×, ≥100×, ≥150× and ≥200× for each region and MAPQ policy. New runs export all thresholds in JSON/TSV metrics. Existing metrics can also render the new columns from their complete, discrete depth histograms; unavailable values remain unavailable.

QC is enabled by default and uses the exact-target BED at `s3://plm-cpm-sandbox-research-usw2/ONT/Oncoseq_refdata/hg38/targets/DARGAN_OCCRA_ALMA_CHLA_wishlist.GENCODE_v50.hg38.bed`. Override it with `--qc_targets_bed` or a per-sample override manifest when using another panel or assembly. The pipeline's enrichment BED still defines enrichment. Removing padding from that BED is not a substitute for supplying the original target intervals.

```bash
nextflow run . -profile clinical \
    --input synthetic_inputs.csv \
    --qc_targets_bed exact_targets.bed \
    --outdir results
```

Use `--qc_extended false` to retain the previous execution without the new scans or target requirement. Use `--qc_methylation false` to retain read/coverage QC without the independent methylation pileup. This release bypasses extended QC in cfDNA mode.

`--qc_sample_manifest` accepts a TSV with `sample` and `targets_bed` columns; each matching row overrides the global target BED. IDs must match original workflow samples, including assigned multiplex sample IDs. Relative paths are resolved against the manifest directory. Unknown or duplicate samples are rejected.

`--qc_cpg_panels` optionally accepts a TSV with `panel_id`, `assembly`, `model_name`, `model_version`, and `bed`. Each feature BED contains one-base forward-C coordinates of reference CpGs. Features are deduplicated within each panel. No model features are downloaded or inferred automatically.

## Measurements and report scope

- Primary-read counts, alignment rate, exact mean/median/N50, and length histograms for all/mapped/unmapped reads and on/off exact targets and enrichment.
- Independent MAPQ 0/20 coverage, zero-inclusive breadth, observed enrichment, mapping-quality attrition, target-base efficiency, and padding burden. Native mosdepth outputs keep their existing definitions.
- Independent CpG valid-depth distributions, reference-denominated coverage, and combined/5mC/5hmC fractions. Valid depth is counted once across modification rows; a modification type not measured remains unavailable.
- Stage-specific observable preprocessing yield, including multiplex run/unclassified material where available.
- Theoretical low-VAF sampling adequacy and optional classifier-panel CpG coverage.

Metrics describe the full mapping-stage BAM, before downstream filtering, phasing, and classifier subsampling. In native time-series mode the QC source remains the original full sample, not each time window. Imported BAMs describe only retained records and cannot reveal earlier yield.

Theoretical sampling adequacy is a binomial coverage calculation, not assay sensitivity, MRD detection, or a clinical negative-result guarantee. Observed enrichment and short off-target reads do not establish adaptive-sampling acceptance or rejection rates. No sequencer exports are required and no clinical pass/fail thresholds are introduced.

## Reading the QC report

The report opens with alignment rate, on-target median read length, mean exact-target depth at MAPQ ≥20, genome-wide median CpG coverage, and the genome-wide fraction of reference CpGs covered at ≥1×. Median CpG coverage is calculated exactly from the existing combined valid-call depth histogram, including uncovered reference CpGs at zero depth. The ≥1× proportion uses the existing reference fraction at one valid call; both CpG widgets use the primary genome (autosomes and sex chromosomes) and count shared 5mC/5hmC coverage once. Component statuses indicate available evidence, not clinical acceptability.

Read-length distributions use separate on/off-target bars with independent percentage y-axis scales. Compare groups using the labeled percentages; the enrichment comparison retains a shared scale. Select **0–10 kb** (500-bp bins), **0–2 kb** (100-bp bins), or **All lengths** (labeled, unequal length bands). Display bins merge whole original bins without smoothing or estimating counts. Every view uses all reads with known length in its own group as the denominator. Zoomed views explicitly report the reads beyond the displayed range; the full-range view includes the entire tail. On/off comparisons exclude unmapped reads. Enrichment comparisons and the complete read-statistics table, including missing-length counts, expand below the plots.

Coverage and methylation depth use compact heatmaps: each cell gives the percentage of reference bases or CpGs meeting a minimum depth. Darker shading indicates a larger fraction, not pass/fail. Coverage defaults to MAPQ ≥20; MAPQ 0 coverage is expandable. Methylation beta distributions are expandable and display only observed primary-genome CpGs with at least 10 valid calls, separately for 5mC and 5hmC.

Theoretical low-VAF adequacy, classifier feature coverage, and detailed regional metrics remain available in expandable sections. In the review report, additional native-tool summaries and gene coverage are grouped separately because their inputs and filtering can differ. All metric JSON/TSV values and analytical calculations are unchanged by this presentation.

## Outputs and private intermediates

Aggregate metric JSON/TSV files are published under `qc_extended/<sample>/`: `metrics.json`, `summary.tsv`, `read_lengths.tsv`, `coverage.tsv`, `targets.tsv`, `methylation.tsv`, `classifier_cpgs.tsv`, `adequacy.tsv`, and `preprocessing.tsv`. Run summaries are `qc_extended/run_summary.tsv` and `qc_extended/run_preprocessing.tsv`.

Fractions in JSON/TSV use `[0,1]`; report tables display percentages. JSON null and TSV `NA` indicate unavailable or undefined values, with explicit component/reason codes. Histograms are aggregate observations in JSON. The metric schema is [metrics.schema.json](../assets/qc_extended/metrics.schema.json).

The genome-wide CpG bedMethyl is also published as `methylation/<sample>.genomewide.cpg.bedmethyl.gz` when methylation QC runs successfully. This BGZF file contains the independent full-BAM modkit pileup, with separate 5mC and 5hmC calls and combined strands. No target BED restricts its coordinates. It contains observed CpGs across the reference genome; unobserved reference CpGs are included in aggregate coverage denominators rather than emitted as synthetic bedMethyl rows. It is independent of Classy's model-specific preprocessing and is not described as the exact classifier input. Disabled/unavailable methylation emits no bedMethyl file.

The genome-wide bedMethyl is an explicitly requested final methylation deliverable. Classifier-specific locus pileups, indexes, normalized references, logs, report fragments, and other intermediates remain in work storage. Report plots are embedded inline in the existing final reports. Existing native output files are retained. With QC enabled, final native HTML publication is routed through an additional report finalization task; QC-only fragments never enter the published native report source directory.

The additional scans and independent pileup add I/O and compute cost. Native Classy continues to perform its own processing. New BAM readers follow the existing isolated-worker policy.

The independent pileup invokes modkit directly in its pinned Classy slim image, which does not provide `python3`. The workflow routes only alignments with complete or partial modification availability to this process; disabled or unavailable methylation bypasses pileup. Python aggregation runs separately in the pinned QC Python runtime.

## Synthetic validation

The fixture builder and acceptance workflow use synthetic sequences and names only:

```bash
python3 tests/qc_extended/build_fixture.py --outdir /tmp/synthetic-qc-fixture
ONCOSEQ_QC_TEST_PATH=/path/to/python/environment/bin \
nextflow run . -main-script tests/subworkflows/local/extended_qc/acceptance.nf \
    -c tests/subworkflows/local/extended_qc/local.config \
    --qc_fixture_dir /tmp/synthetic-qc-fixture \
    --qc_targets_bed /tmp/synthetic-qc-fixture/targets.bed \
    --outdir /tmp/synthetic-qc-results
python3 tests/qc_extended/check_publication.py /tmp/synthetic-qc-results
```

This local acceptance disables methylation and validates actual read/coverage calculations, sample completion, report content, and publication. The separate production-image methylation smoke test requires a supported Linux amd64 container environment; it must execute actual modkit rather than use a stub. See the detailed [implementation handoff](qc-expansion-handoff.md) for numerical definitions and test requirements.

Python QC tasks use the pinned review runtime's packaged interpreter (`/opt/wave/.pixi/envs/default/bin/python3`). Depth is streamed from the Samtools implementation bundled with pysam, and both versions are recorded. This avoids assuming that a Dorado/preprocessing image or the task PATH provides the required Python package or a separate Samtools executable.
