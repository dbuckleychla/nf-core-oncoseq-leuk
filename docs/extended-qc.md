# Extended sequencing QC

Cellular-DNA runs now add regional read-length, coverage, and independent methylation QC. The new metrics appear in both the native report and the optional review report's QC tab. Existing analytical steps and outputs are retained.

QC is enabled by default and requires an **explicit exact-target BED**, supplied globally or through a per-sample override manifest. The pipeline's enrichment BED still defines enrichment. Removing padding from that BED is not a substitute for supplying the original target intervals.

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

## Outputs and private intermediates

Aggregate metric JSON/TSV files are published under `qc_extended/<sample>/`: `metrics.json`, `summary.tsv`, `read_lengths.tsv`, `coverage.tsv`, `targets.tsv`, `methylation.tsv`, `classifier_cpgs.tsv`, `adequacy.tsv`, and `preprocessing.tsv`. Run summaries are `qc_extended/run_summary.tsv` and `qc_extended/run_preprocessing.tsv`.

Fractions in JSON/TSV use `[0,1]`; report tables display percentages. JSON null and TSV `NA` indicate unavailable or undefined values, with explicit component/reason codes. Histograms are aggregate observations in JSON. The metric schema is [metrics.schema.json](../assets/qc_extended/metrics.schema.json).

Only the final aggregate metrics are published by new QC tasks. QC-only pileups, indexes, normalized references, logs, report fragments, and other intermediates remain in work storage. Report plots are embedded inline in the existing final reports. Existing native output files are retained. With QC enabled, final native HTML publication is routed through an additional report finalization task; QC-only fragments never enter the published native report source directory.

The additional scans and independent pileup add I/O and compute cost. Native Classy continues to perform its own processing. New BAM readers follow the existing isolated-worker policy.

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
