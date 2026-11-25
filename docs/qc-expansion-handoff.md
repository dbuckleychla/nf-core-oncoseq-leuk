# Additive cellular-DNA QC: comparison and implementation handoff

## 1. Purpose, authority, and delivery boundary

This document specifies an additive QC expansion for adaptive-sampling profiling of cellular blood/marrow DNA. It is intended for implementation by agents that need explicit definitions, interfaces, small work packages, and numerical acceptance criteria. It is a design specification, not a task-status tracker; subsequent implementation work belongs in Beads.

**The action that produced this document is documentation only. No pipeline implementation, configuration change, tracker mutation, commit, or push is part of that action.** The implementation instructions below are for a later, separately initiated session.

The user resolved the following requirements during planning:

1. Preserve native nf-core-oncoseq analyses and their current outputs. Add steps instead of replacing or changing native steps.
2. The publication restriction applies to the **QC add-on**. All current pipeline outputs remain publishable according to their existing settings. Only aggregate QC metric files from the add-on may be published; QC-only intermediates must stay in work storage.
3. Ignore cfDNA for this expansion. Prioritize cellular blood/marrow DNA and adaptive sampling.
4. Enable the new QC by default, with an explicit disable option.
5. Require an explicit exact-target BED for enabled cellular-DNA QC. Do not reconstruct it by removing enrichment padding.
6. Do not rely on sequencer output files. There are no MinKNOW, sequencing-summary, Readfish, POD5-metadata, or decision-log requirements for these QC calculations. Existing BAM tags may provide optional technical metadata.
7. Optional externally supplied classifier CpG panels may enable feature-coverage QC.

The execution request subsequently added QC summaries and distributions to both native and review reports. This supersedes the original metric-only presentation boundary: additional report sections may embed aggregate QC and inline plots, while QC-only intermediate files remain unpublished. No additional caller, tumor classifier, clinical cutoff, or biological reinterpretation is included. Native analytical computations remain unchanged.

## 2. Inspected baselines and evidence

### 2.1 Baselines

Inspection date: **25 November 2026**.

| Source                           | Revision                                             | Role                                              |
| -------------------------------- | ---------------------------------------------------- | ------------------------------------------------- |
| This fork                        | `6571a2811c47a82d4c9f2a413c309e8c3f65bd0e`           | Current implementation baseline                   |
| FrankONTstein                    | `46553118352caa218f82762f74c74ce635220177`           | QC behavior and reusable implementation reference |
| Original nf-core-oncoseq release | `v1.0.2`, `a663a10012e467ded972f1b98ded372d5f5ef2fa` | Historical upstream baseline only                 |

The sibling FrankONTstein checkout had an uncommitted `conf/modules.config` modification. Do not port that local configuration change. Use the pinned committed sources below. The current fork was clean at inspection.

[The existing comparison](frankontstein-comparison.md) is useful migration context, but compares an earlier baseline. It must not be treated as a current inventory of everything present or absent in this fork.

### 2.2 Source map

These links are pinned so the handoff remains interpretable after later edits:

| ID  | Source                                                                                                                                                                                     | Evidence                                                                                                 |
| --- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------- |
| F1  | [FrankONTstein QC documentation](https://github.com/dbuckleychla/FrankONTstein/blob/46553118352caa218f82762f74c74ce635220177/documentation/qc.md)                                          | Published outputs, coverage and CpG definitions, bounded execution                                       |
| F2  | [FrankONTstein QC engine](https://github.com/dbuckleychla/FrankONTstein/blob/46553118352caa218f82762f74c74ce635220177/bin/adaptive_qc.py)                                                  | `read_length_groups`, `lengths_summary`, `contig_qc`, `bam_qc`, `methylation_qc`, `parallel_methylation` |
| F3  | [FrankONTstein QC modules](https://github.com/dbuckleychla/FrankONTstein/blob/46553118352caa218f82762f74c74ce635220177/modules/local/qc/main.nf)                                           | `MODKIT_PILEUP`, `INDEX_BEDMETHYL`, `SAMPLE_QC`, `DEMUX_QC`                                              |
| F4  | [FrankONTstein workflow](https://github.com/dbuckleychla/FrankONTstein/blob/46553118352caa218f82762f74c74ce635220177/workflows/adaptive.nf)                                                | Independent pileup and QC attachment; classifier remains separate                                        |
| F5  | [FrankONTstein QC tests](https://github.com/dbuckleychla/FrankONTstein/blob/46553118352caa218f82762f74c74ce635220177/tests/test_qc.py)                                                     | Synthetic numerical oracles and overlap tests                                                            |
| F6  | [FrankONTstein real-container smoke test](https://github.com/dbuckleychla/FrankONTstein/blob/46553118352caa218f82762f74c74ce635220177/tests/run_qc_container_smoke.py)                     | Actual modkit execution rather than stub wiring                                                          |
| F7  | [FrankONTstein input checker](https://github.com/dbuckleychla/FrankONTstein/blob/46553118352caa218f82762f74c74ce635220177/bin/check_bam.py)                                                | Bounded presence inspection, not necessarily whole-input yield                                           |
| F8  | [FrankONTstein ONT metadata extraction](https://github.com/dbuckleychla/FrankONTstein/blob/46553118352caa218f82762f74c74ce635220177/bin/ont_bam_metadata.py)                               | Header provenance and first-1,000-record tag inspection                                                  |
| N1  | [Current main workflow](https://github.com/dbuckleychla/nf-core-oncoseq-leuk/blob/6571a2811c47a82d4c9f2a413c309e8c3f65bd0e/workflows/oncoseq.nf)                                           | Mapping, time-series routing, independent coverage, classifier subsampling                               |
| N2  | [Current mapping workflow](https://github.com/dbuckleychla/nf-core-oncoseq-leuk/blob/6571a2811c47a82d4c9f2a413c309e8c3f65bd0e/subworkflows/local/mapping/mapping.nf)                       | Cramino, route-dependent SeqKit, review-gated alignment QC                                               |
| N3  | [Current adaptive coverage](https://github.com/dbuckleychla/nf-core-oncoseq-leuk/blob/6571a2811c47a82d4c9f2a413c309e8c3f65bd0e/subworkflows/local/adaptive_specific/coverage_separate.nf)  | Background/panel mosdepth tasks and native reporting                                                     |
| N4  | [Current module configuration](https://github.com/dbuckleychla/nf-core-oncoseq-leuk/blob/6571a2811c47a82d4c9f2a413c309e8c3f65bd0e/conf/modules.config)                                     | Flag masks, MAPQ thresholds, broad default publication                                                   |
| N5  | [Current initialization](https://github.com/dbuckleychla/nf-core-oncoseq-leuk/blob/6571a2811c47a82d4c9f2a413c309e8c3f65bd0e/subworkflows/local/utils_nfcore_oncoseq_pipeline/main.nf)      | Per-sample enrichment BEDs and negative-padding construction of `bed_nopad`                              |
| N6  | [Current methylation classification](https://github.com/dbuckleychla/nf-core-oncoseq-leuk/blob/6571a2811c47a82d4c9f2a413c309e8c3f65bd0e/subworkflows/local/methylation_analysis/classy.nf) | Existing Classy analyses and outputs                                                                     |
| N7  | [Current alignment QC documentation](https://github.com/dbuckleychla/nf-core-oncoseq-leuk/blob/6571a2811c47a82d4c9f2a413c309e8c3f65bd0e/docs/alignment-qc.md)                              | Correct primary-read alignment denominator and full-sample scope                                         |
| N8  | [Current review presentation](https://github.com/dbuckleychla/nf-core-oncoseq-leuk/blob/6571a2811c47a82d4c9f2a413c309e8c3f65bd0e/bin/review/presentation.py)                               | Existing Cramino, SeqKit, NASVAR, coverage, and review summaries                                         |

### 2.3 Delta matrix

| QC capability                 | FrankONTstein at inspected revision                                    | Current fork at inspected revision                                                 | Required action                                                       |
| ----------------------------- | ---------------------------------------------------------------------- | ---------------------------------------------------------------------------------- | --------------------------------------------------------------------- |
| Overall length/yield          | Primary-read counts, sequenced bases, exact mean/median/N50, histogram | Cramino and route-dependent SeqKit already provide general length/yield statistics | Add consistent scoped metrics; do not replace existing statistics     |
| Alignment rate                | Mapped primary / all primary, including unmapped denominator           | Correct separate `SAMTOOLS_ALIGNMENT_QC`, enabled with `report_review`             | Make add-on rate independent of reporting; leave existing task alone  |
| On/off-target lengths         | Exact-target and enrichment comparisons separately                     | No equivalent comprehensive regional length output found                           | Restore all groups                                                    |
| Coverage filters              | MAPQ 0/20; primary; includes duplicate/QC-failed flags                 | Native presets include MAPQ 0/60 and masks `1540`/`1796`                           | Add independent comparable calculations                               |
| Exact-target definition       | Independent exact-target and enrichment assets                         | `bed_nopad` is generated by eroding enrichment intervals                           | Require an independent QC target BED                                  |
| Coverage breadth/distribution | Merged region denominators, zeros, original target rows                | Existing mosdepth coverage and plots provide partial coverage                      | Add explicit union and per-target breadth/histograms                  |
| Methylation QC                | Separate CpG pileup and reference-denominated summaries                | Classy processing exists, but no equivalent independent comprehensive QC product   | Add separate modkit QC                                                |
| 5mC/5hmC                      | Separate and combined valid-call fractions                             | Classifier predictions are not those metrics                                       | Restore with availability states                                      |
| CpG coverage                  | Reference CpGs and valid depths ≥1/5/10/20                             | Not comprehensively surfaced                                                       | Restore; add off-enrichment scope                                     |
| Preprocessing yield           | Demux/unclassified counts, prepared-input check, trim metadata         | Pass/fail SeqKit and source artifacts available in some routes                     | Add explicit stage/scope accounting                                   |
| ONT provenance                | Header metadata plus bounded tag sample                                | Partial technical provenance and tag checks                                        | Add selected technical fields and exhaustive alignment-stage counters |
| Identity/other QC             | Tool-specific metrics and reporting                                    | Cramino identity statistics and NASVAR QC already exist where those tools run      | Acknowledge as present; no replacement or reinterpretation            |
| Reports                       | Standalone HTML plus JSON/TSV; bedMethyl published                     | Native reports plus optional clinical review                                       | JSON/TSV plus existing report integration; retain current outputs     |
| Acceptance/rejection          | Explicitly not reconstructed from BAM                                  | No equivalent direct decision accounting                                           | Exclude; observed enrichment is not an acceptance rate                |

Important qualifications:

- F7 inspects at most a bounded number of reads; its `reads`/`bases` can be null when the check did not reach EOF. Do not describe those as full input yield.
- F8's first 1,000 records are not a representative sample. Do not port a sampled tag fraction as a whole-BAM percentage.
- Classy may receive a subsampled BAM in the current workflow. QC of the full mapping BAM cannot establish how many features actually entered that classifier execution.
- Removing padding from previously merged enrichment intervals cannot generally reconstruct original targets. The new target requirement is intentional.
- Do not copy F2's report rendering, raw header export, unconditional output directory publishing, or any ambiguous missing-value behavior.

## 3. Non-negotiable implementation invariants

1. Attach new work to existing outputs. Do not alter native basecalling, alignment, filtering, phasing, callers, coverage presets, Classy inputs, or reporting commands.
2. Do not replace native Classy pileups with the new pileup. Their inputs, filters, and model expectations differ.
3. Preserve native process names, scripts, container choices, `ext.args`, and publication settings. New names use the `QC_EXT_` prefix.
4. Orchestration may gain includes, branches, new parameters, and additional subworkflow emits. Existing process output declarations should not need edits: expose already-emitted artifacts at the subworkflow boundary.
5. Never depend on a published output path. Wire actual process channels, including optional preprocessing statistics.
6. QC defaults to one full mapping-stage sample. Do not append native time-window suffixes or imply that metrics describe each time window.
7. Preserve unmapped records in counts when present. Imported BAMs describe only retained records, not original instrument yield.
8. Absent optional inputs never cause sample loss through an inner join.
9. Publish only the explicit aggregate metric allowlist in section 10. Work storage remains private runtime storage and is still subject to normal retention policy.
10. Use synthetic fixtures and identifiers only in the public repository, tests, documentation, tracker, commits, and PRs. Patient results belong in private runtime storage.

## 4. User-facing configuration and input validation

### 4.1 Parameters

Add declarations to the typed `params` block in `main.nf`, defaults in `conf/params.config`, definitions in `nextflow_schema.json`, and usage documentation. Use underscores consistently; introduce no additional hyphenated aliases.

| Parameter            | Typed declaration | Default | Meaning                                               |
| -------------------- | ----------------- | ------- | ----------------------------------------------------- |
| `qc_extended`        | `Boolean`         | `true`  | Enable the new cellular-DNA QC branch                 |
| `qc_methylation`     | `Boolean`         | `true`  | Enable the branch's independent methylation component |
| `qc_targets_bed`     | `Path?`           | null    | Global exact-target BED                               |
| `qc_sample_manifest` | `Path?`           | null    | Optional per-sample exact-target overrides            |
| `qc_cpg_panels`      | `Path?`           | null    | Optional classifier-feature panel manifest            |

Effective enablement is `qc_extended && !cfdna`. If false, schedule no new QC tasks and do not require or inspect QC reference inputs. `qc_methylation=false` retains alignment/coverage QC and emits methylation status `disabled`. No new cfDNA outputs are produced in this release.

With effective enablement true, every original sample must have an explicit exact-target BED. A missing BED is a configuration error, not permission to use `bed_nopad`. The user can retain the prior pipeline behavior with `--qc_extended false`.

### 4.2 Target override manifest

Use a UTF-8 TSV with exactly these required columns:

```tsv
sample	targets_bed
SYNTHETIC_A	targets_a.bed
SYNTHETIC_B	targets_b.bed
```

Rules:

- Match `sample` against the original sample IDs after demultiplexing assignment and before time-series suffixes. Match strings exactly; do not normalize case or strip substrings.
- A manifest row overrides `qc_targets_bed` for that sample. Samples without a row use the global BED. Reject unresolved samples, duplicate rows, unknown IDs, empty fields, and extra columns.
- Resolve relative paths against the manifest's directory, including URI-relative paths for remote manifests. Absolute paths/URIs retain their meaning.
- Use Nextflow `file`/`path` staging for remote assets; Python helpers only read staged local files.
- Do not extend or reinterpret the native input samplesheet in this release.
- Validate sample-key uniqueness and output-directory safety. Reject path separators, `.`/`..`, and control characters in add-on output identifiers; never silently rewrite a sample ID.

### 4.3 Optional classifier CpG panel manifest

Use a UTF-8 TSV with columns `panel_id`, `assembly`, `model_name`, `model_version`, `bed`:

```tsv
panel_id	assembly	model_name	model_version	bed
SYNTHETIC_PANEL	hg38	SYNTHETIC_MODEL	1	features.bed
```

`panel_id` must be unique and safe for use as a table key. It does not become a patient identifier. `assembly` must match the run's declared assembly. BED rows are one-base, zero-based, half-open intervals naming the forward-strand C of a reference CpG. Validate the actual `CG` sequence in the supplied FASTA. Deduplicate identical coordinates within a panel and report input rows, duplicate rows, and unique features. Features shared between panels are counted independently per panel.

Do not automatically mine model files, infer model versions, lift coordinates, or download panels. An incompatible or malformed supplied panel is a validation error. No panel supplied is `unavailable/no_panel_supplied`, not a QC failure. The manifest is independent of native classifier selection.

### 4.4 General reference validation

- Use the same reference FASTA/FAI as mapping. Validate coordinate-sorted, readable BAM and usable index, and exact BAM/FASTA sequence-dictionary name/length agreement; ordering differences alone need not fail.
- Do not silently drop mismatched contigs, auto-add `chr`, or liftover intervals. A mismatch fails the QC task with an actionable message.
- Accept BED3 or greater, plain text or gzip/BGZF. Ignore blank lines and BED comments/track/browser lines. Require integer coordinates `0 <= start < end <= contig_length`.
- Reject an empty exact-target union and empty enrichment union. Empty off-enrichment space is allowed and yields null ratio metrics.
- Preserve the original target row order and labels. Use `line_number:label` as the stable target key, or `line_number:contig:start-end` when column 4 is absent.
- Normalize and merge overlapping or adjacent intervals for aggregate denominators. Keep original rows separate for per-target summaries.
- Calculate `targets_outside_enrichment_bases` and its target-union fraction; emit a factual warning if nonzero. Do not require target containment or clip the exact targets.
- Record SHA-256 hashes of supplied BED bytes, normalized interval representation, and reference bytes. Compute reference summaries/checksums once per shared reference asset, not once per sample. Do not add a full BAM hashing pass.

## 5. Workflow interfaces and placement

### 5.1 Proposed code layout

Create new components under these locations; names below are the intended interfaces, not suggestions for unrelated refactoring:

| Location                                 | Responsibility                                                                                                      |
| ---------------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| `subworkflows/local/extended_qc/main.nf` | `EXTENDED_QC` orchestration and optional-component completion                                                       |
| `modules/local/extended_qc/main.nf`      | All `QC_EXT_*` processes, explicit outputs, fixed tool commands                                                     |
| `bin/extended_qc.py`                     | CLI entrypoint with subcommands `reference`, `alignment`, `preprocessing`, `methylation`, `assemble`, `run-summary` |
| `bin/qc_extended/`                       | Python modules `contracts.py`, `regions.py`, `alignment.py`, `methylation.py`, `preprocessing.py`, `assemble.py`    |
| `assets/qc_extended/metrics.schema.json` | Versioned output JSON Schema and metric registry contract                                                           |
| `tests/qc_extended/`                     | Synthetic builders, pure-Python tests, numerical oracles, container smoke driver                                    |
| `tests/modules/local/extended_qc/`       | Process execution, publication, and resource tests                                                                  |
| `tests/subworkflows/local/extended_qc/`  | Multiple samples, optional inputs, statuses, and channel tests                                                      |
| `workflows/tests/extended_qc.nf.test`    | Native/add-on integration and scope regression tests                                                                |
| `conf/profiles/test_qc_extended.config`  | Self-contained synthetic pipeline acceptance profile                                                                |

Retain the FrankONTstein MIT attribution for adapted algorithms, using the existing `scripts/LICENSE.FrankONTstein` where appropriate and a source-revision comment in adapted modules. Do not require a sibling checkout at runtime or import its code dynamically.

### 5.2 Attachment point and scope

In `workflows/oncoseq.nf`, add `EXTENDED_QC` as a consumer of `MAPPING.out.bam`, the existing reference channels, and original per-sample enrichment BEDs. Do not consume `ch_bam_for_calling`, haplotagged BAMs, or the Classy subsampling output.

All new sample outputs use `input_scope=full_mapping_stage_sample`. In native time-series mode, one QC sample corresponds to the original full mapping-stage BAM; no per-window recalculation or duplicate output is performed. Document this beside every run summary. A later window-specific extension is out of scope.

### 5.3 Channel contracts

Use an explicit scalar key for joins, not equality of entire metadata maps. Retain `meta` unchanged as data attached to that key.

| Channel                      | Tuple contract                                                | Cardinality                                                |
| ---------------------------- | ------------------------------------------------------------- | ---------------------------------------------------------- |
| `ch_qc_samples`              | `(sample_id, meta, bam, index)`                               | Exactly one per original mapped sample                     |
| `ch_qc_regions`              | `(sample_id, enrichment_bed, exact_targets_bed)`              | Exactly one per enabled sample                             |
| `ch_qc_reference`            | `(reference_key, fasta, fai, assembly)`                       | Current pipeline reference; shared value input             |
| `ch_qc_preprocessing_inputs` | `(scope_type, scope_id, stage, format, files, associations)`  | Zero or more existing artifacts per sample/run             |
| `ch_qc_reference_summary`    | `(reference_key, normalized_regions, reference_summary_json)` | Reusable for identical reference/region/panel combinations |
| `ch_qc_alignment`            | `(sample_id, meta, alignment_json, targets_partial_tsv)`      | Exactly one per sample                                     |
| `ch_qc_mod_status`           | `(sample_id, meta, status_json)`                              | Exactly one per sample, including skip reasons             |
| `ch_qc_methylation`          | `(sample_id, meta, methylation_json)`                         | Exactly one completion record per sample                   |
| `ch_qc_components`           | `(sample_id, component_name, metric_file)`                    | Small metric files only                                    |
| `ch_qc_final`                | `(sample_id, meta, metric_files)`                             | Exactly one final publication per sample                   |

The reference helper may internally key raw reference-derived CpG inventories separately from region-specific summaries so that different panels do not force repeated genome scans. Internal CpG positions/normalized BEDs are private intermediate files and may be cached in work storage.

Validate mandatory key sets and duplicate keys before expensive QC tasks. For optional artifacts, construct an explicit empty/disabled/unavailable completion record from the expected sample registry. Do not inner-join an empty optional channel. Collect only the small reference/sample registry and component-file metadata in the Nextflow driver; never BAM contents or per-base records.

Data flow:

```text
existing mapping BAM/index ────> QC_EXT_ALIGNMENT ────────────────┐
existing FASTA + E + exact T ──> QC_EXT_REFERENCE ───────────────┤
                                    │                         │
alignment modification status ──────┴─> QC_EXT_MODKIT_PILEUP   │
                                           │ private pileup  │
                                           v                 │
                                      QC_EXT_METHYLATION ─────┤
existing preprocessing artifacts ──> QC_EXT_PREPROCESSING ─────┤
disabled/unavailable component records ───────────────────────┤
                                                             v
                                                     QC_EXT_ASSEMBLE
                                                             │
                                                aggregate QC JSON/TSV
                                                             v
                                                     QC_EXT_RUN_SUMMARY
```

### 5.4 Preprocessing taps

Expose already-existing artifacts through additional emits in `BASECALL_SIMPLEX` and `BASECALL_MULTIPLEX`; do not edit Dorado, Samtools, or SeqKit process bodies. The existing `BASECALL_SHARDED.out.ubam` supplies merged prepared-input material. Multiplex demultiplexed directories must be tapped before native filtering discards unclassified outputs. Post-trim BAMs already come from `DORADO_TRIM`.

Carry run/sample/barcode relationships explicitly from the existing demultiplexing join. Do not recover them with `.contains('pass')`, arbitrary `.replace('_pass', '')`, or filename substring matching. Where native output metadata already adds a known stage suffix, use the exact association established at the producing channel boundary and validate it; do not change native naming behavior.

Reuse already-computed SeqKit pass/fail totals whenever their scope and stage match. Add a streaming sidecar count only for stages not already counted. For imported FASTQ/uBAM/aligned BAM, report only observable stages. Do not generate missing preprocessing stages.

## 6. Metric definitions and algorithms

### 6.1 Shared conventions

Fractions are stored on `[0,1]`, not as percentages. Counts and base totals are nonnegative integers. Units and numerator/denominator definitions are part of the metric registry. JSON null and TSV `NA` represent undefined values; never serialize NaN, infinity, or the string `null` as a number.

An observed empty group has count zero, but its mean, median, N50, and normalized histogram fractions are null. A missing component is unavailable, not an observed empty group. Distinguish both in status fields.

Primary means neither secondary (`0x100`) nor supplementary (`0x800`). Duplicate (`0x400`) and QC-failed (`0x200`) primary records remain included in the parity metric family, and their flag counts are reported separately. One primary record is not a proven independent molecule: duplex outputs, split reads, and unmarked duplication can violate that interpretation.

### 6.2 Alignment and read-length metrics

Do one full sequential alignment-record pass for read statistics and modification-tag integrity. This avoids losing or double-counting unmapped reads with assigned coordinates and naturally visits each record once. Indexed contig passes are used separately for depth. Do not combine a full pass with an additional `fetch('*')` count.

Required counters: total alignment records, primary records, mapped primary, unmapped primary, secondary records, supplementary records, records carrying both nonprimary flags, duplicate-flagged primary, QC-failed primary, missing-length primary, hard-clipped primary, and supported modification-tag categories. For total-record reconciliation, use mutually exclusive primary/secondary-only/supplementary-only/both categories; flag-occurrence counters may overlap and must be labeled accordingly.

`alignment_rate = mapped_primary_records / primary_records`, null when there are no primary records. Include both QC-passed and QC-failed records. Never substitute Cramino's `% from total alignments` for this denominator.

Calculate these seven length groups independently:

```text
all, mapped, unmapped,
on_target, off_target,
on_enrichment, off_enrichment
```

For a mapped primary record, classify overlap using only reference blocks containing aligned query bases (`M`, `=`, `X`). Insertions, deletions (`D`), skips (`N`), and soft clips do not establish overlap. Any qualifying overlap assigns the read to the on group; mapped reads without overlap belong to its off group. Unmapped reads belong to neither off group. A read can be on target and on enrichment simultaneously; never sum these independent comparisons.

Use sequence length when present; otherwise infer query length from a valid CIGAR if possible. Include soft clips; exclude unrecoverable hard clips. Store `length_source` counts (`sequence`, `cigar_inferred`, `unavailable`). If neither source exists, preserve the read count but exclude it from length-based denominators and increment `missing_length_records`.

Per group store total records, records with known lengths, missing lengths, sequenced bases with known length, min/max, exact arithmetic mean, exact median, exact N50, and 100-bp histogram counts. Exact statistics use a length-frequency counter, not binned midpoints. Histogram bin `[100*k,100*(k+1))` includes its lower endpoint; the last observed bin may be arbitrarily long. A bin's read fraction uses known-length records in that group. Also expose each bin's sequenced-base sum.

N50 is the largest descending length at which cumulative sequenced bases reach at least half the total. For an even count, median is the mean of the two central lengths. A total of zero known sequenced bases makes N50 null.

Additional descriptive metrics:

- Full mapped-primary MAPQ histogram; separately count MAPQ 255 as unavailable mapping-quality semantics. Threshold coverage below remains numerically filtered for parity and must state that 255 passes that numerical filter.
- Soft-clipped query-base sum and fraction over mapped-primary known query bases; hard-clipped record count, with no recovered length claim.
- `qs` distribution/mean on primary records carrying a finite valid Dorado `qs` tag, with denominator and missing-tag count. Do not silently mix this with inferred base-quality averages.
- For records with complete base qualities, mean error probability and Phred conversion `-10*log10(mean(10^(-Q/10)))`, stored under a separate metric family. Missing/sentinel base qualities are excluded and counted. Do not call arithmetic mean Phred the same quantity.

### 6.3 Coverage engine

Use two streamed `samtools depth` passes per contig at MAPQ thresholds 0 and 20. Base quality is 0. The parity command shape is:

```bash
samtools depth -r CONTIG -q 0 -Q MAPQ -g 0x600 -G 0x900 INPUT.bam
```

Execute as an argument array from Python, not a shell interpolation. Verify the pinned Samtools behavior with duplicate/QC-failed/secondary/supplementary fixtures. `-q` is base quality and `-Q` is mapping quality for **depth**. Do not add `-J`, overlap removal, a maximum-depth truncation, or pileup BAQ behavior.

Stream depth into counters. Avoid a full per-base depth output file. Maintain merged-interval membership and original-target membership with half-open interval events. Accumulate positive-depth histogram counts and recover zeros as region length minus observed covered positions. Check this difference is nonnegative. Test coincident start/end events and overlapping duplicate target rows.

Aggregate coverage uses primary human chromosomes 1–22/X/Y. Support UCSC, bare chromosome, and the RefSeq primary accession families explicitly tested in F5. Do not identify every `NC_*` sequence as a primary chromosome. Keep mitochondrial/alternate/decoy/other contigs separate. Retain all supplied target rows in `targets.tsv` even when outside the primary-genome aggregate.

For each MAPQ and scope (`genome_primary`, `enrichment`, `targets`, `off_enrichment`, plus per-contig scopes), store:

- `reference_bases`: union length in the declared contig universe, including reference Ns.
- `aligned_bases`: sum of depth across those reference positions.
- `mean_depth = aligned_bases/reference_bases`.
- Exact depth histogram including depth zero.
- `breadth_ge_D = positions_with_depth_at_least_D/reference_bases` for `D = 1,5,10,20,30,50,100`.
- `zero_coverage_fraction`, exact depth median, population depth CV, and fraction of bases at least `0.2 * mean_depth` when mean depth is positive. CV is `sqrt(E[d^2]-E[d]^2)/mean`; clamp only floating-point roundoff below zero.
- `mean_depth_mapq20_over_mapq0`, null when the MAPQ0 mean is zero or unavailable.

Per-target metrics use each original interval independently, including zero-covered positions. Aggregate targets use the union, not a sum of original rows. Never average interval means without length weighting.

`observed_enrichment_depth_ratio = mean_depth(enrichment)/mean_depth(off_enrichment)`, separately at MAPQ0 and MAPQ20. Null when the denominator is zero or unavailable. Off-enrichment is the primary-genome coordinate complement of enrichment, not a BAM of reads that never touched enrichment. This distinction matters for long reads that span boundaries.

The original masks remain unchanged: native `MOSDEPTH_NOFILT` uses `-F 1540 -Q 0`, native primary/background use `-F 1796 -Q 0`, and native unique uses `-F 1796 -Q 60`. These results are not interchangeable with this family.

### 6.4 Novel target allocation and short/long-read metrics

For MAPQ0 and MAPQ20, measure aligned query/reference bases within exact-target union and within enrichment-minus-target union, using the same coverage policy. Use primary-genome scopes for these ratios and state that universe in each row.

| Metric ID                             | Formula                                                                              | Interpretation                                                                             |
| ------------------------------------- | ------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------ |
| `target_aligned_base_fraction`        | Target aligned bases / all primary-genome aligned bases                              | Useful coordinate allocation among mapped sequence                                         |
| `target_yield_per_primary_query_base` | Target aligned bases / known-length bases of all primary records, including unmapped | Target yield relative to observed sequencing yield; mixed contig scope explicitly recorded |
| `padding_aligned_base_fraction`       | Aligned bases in enrichment minus targets / enrichment aligned bases                 | Yield allocated to enrichment padding rather than exact targets                            |
| `targets_outside_enrichment_fraction` | Target-union bases outside enrichment / target-union bases                           | Reference-design mismatch, not sequencing performance                                      |

For the query-base yield metric, all primary lengths must be known; otherwise emit unavailable with `incomplete_query_length_denominator`. Do not divide an exhaustive aligned-base numerator by a partial query-base denominator.

For each length group, report read-count and sequenced-base fractions for lengths **≤500 bp**, **≤1,000 bp**, and **≥10,000 bp**, using known-length records/bases and an explicit missing-length count. These nested descriptive thresholds are not mutually exclusive categories. Do not name them `rejected`, `accepted`, `false_accept`, or `false_reject`.

### 6.5 Theoretical low-VAF sampling adequacy

Compute from the MAPQ20 exact-target depth histogram, retaining zero-covered bases. For assumed variant allele fraction `v` and required alternate-read support `k`:

```text
P(X >= k | n,v) = 1 - sum(j=0..k-1) choose(n,j) v^j (1-v)^(n-j)
```

Use VAFs `{0.01,0.02,0.05,0.10}`, supports `{3,5}`, and probability target `0.95`. Find the smallest integer depth meeting the probability criterion, using stable binomial recurrence/log arithmetic. Report the fraction of exact-target union bases with depth at least that threshold, and corresponding per-target fractions. No caller VCF is required.

Numerical oracles:

| Assumed VAF | Depth for ≥3 supports with 95% probability | Depth for ≥5 supports with 95% probability |
| ----------- | ------------------------------------------ | ------------------------------------------ |
| 1%          | 628                                        | 913                                        |
| 2%          | 313                                        | 456                                        |
| 5%          | 124                                        | 181                                        |
| 10%         | 61                                         | 89                                         |

Call this **theoretical sampling adequacy**, never sensitivity, limit of detection, MRD detection, or a clinical negative-result guarantee. It ignores error, strand imbalance, allele-specific mapping, copy number, purity, molecule dependence, and caller filtering. VAF is an assumed allele fraction, not tumor fraction. Count duplicate/duplex caveats remain applicable. Native callers and thresholds are unchanged.

## 7. Independent methylation and classifier-feature QC

### 7.1 Containers and command contract

Use these inspected immutable images for the initial compatibility baseline of **new processes only**:

| New process family                                               | Image                                                                                                        |
| ---------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| Reference, alignment, preprocessing, aggregation, index handling | `ghcr.io/dbuckleychla/frankontstein@sha256:15c315c121d37b5c115d3d965836e60d56756498ca2deebe066a373d41e96785` |
| Independent methylation pileup                                   | `ghcr.io/chusj-pigu/classy-slim@sha256:7ec5fc1e4607be13d30dc65a1b07f2844f8289edbec0ca67a70e9b88ef8beb4c`     |

Verify Python, pysam, Samtools, and modkit availability and actual versions in real Linux containers before merging. These are reuse candidates established by inspected source, not a claim that they were executed during this documentation session. If a required capability is absent, stop that implementation package and document the concrete incompatibility; do not silently substitute a mutable image or alter native Classy.

Start from F3's command, with deterministic threshold-sampling settings made explicit after the capability gate:

```bash
modkit pileup INPUT.bam sample.cpg.bedmethyl.gz \
    --ref reference.fa \
    --cpg --modified-bases 5mC 5hmC --combine-strands --bgzf \
    --threads ALLOCATED_CPUS \
    --num-reads 10042 --filter-percentile 0.1 \
    --log-filepath modkit.log
```

Treat support for all explicit options as part of the container smoke test. F3's original invocation did not explicitly pin sampling options; therefore record this as an intentional reproducibility refinement, and compare threshold/output behavior on the synthetic fixture rather than claiming byte-for-byte historical identity.

Do not add an invented `--min-mapq` option: the inspected current Modkit advanced-usage documentation did not list it. This release has one independent pileup using the pinned tool's documented alignment policy. It does **not** claim MAPQ20 methylation. Store a methylation-specific filter policy separately from the Samtools coverage policies and verify primary/secondary/supplementary/duplicate/QC-failed behavior in a container fixture. Any such tool-native difference must appear in provenance.

The 0.1 setting is an estimated-confidence filtering percentile, not a fixed methylation probability cutoff of 0.1 or 0.9. Capture the actual estimated thresholds when available through a version-tested parser; otherwise report `threshold_values_unavailable` while retaining the explicit settings. Do not scrape and publish arbitrary log text.

The 8,000-record maximum is an explicit pileup-engine cap, not a maximum for Samtools depth. Record it. If source coverage suggests that it may constrain results, emit `potential_methylation_depth_cap`; label the warning as potential because the tool's eligible-record universe can differ. Never present capped CpG depth as an exhaustive count. Higher-depth behavior requires a separately validated change rather than silently enabling tool downsampling.

Index the private BGZF output with tabix through pysam, then aggregate by contig. The pileup, index, log, tool version scratch file, and any helper files have publication disabled. The new pileup never enters a classifier or variant-calling channel.

### 7.2 Modification availability and integrity

During the full primary-record scan, count at least:

- Presence of MM/Mm, ML/Ml, and MN tags, with canonical names and any legacy spelling explicitly recorded.
- Supported cytosine modification codes observed (`m` and `h`) and records with probabilities usable for pileup.
- Missing probabilities, missing MN, explicit/implicit skip modes, empty valid encodings, unsupported cytosine codes, and malformed encodings.
- MN/sequence inconsistencies, invalid modification positions, invalid ML array cardinality, and conflicting duplicate tag spellings.

Use a standards-aware parser, with pysam's decoder behavior tested against the SAM tag specification. Never hand-parse modification coordinates as reference offsets: MM positions are relative to the query's canonical bases and orientation semantics matter. Missing MN is not automatically corruption. Empty MM encodings and absent probability arrays are not automatically malformed; distinguish legal-but-insufficient-for-this-calculation from invalid encoding. Do not invent modification probabilities for missing ML.

Availability is based on actual supported encoding/provenance, not the presence of a zero-valued `h` row produced because the command requested 5hmC. If a BAM only measures 5mC, mark 5hmC unavailable even if modkit emits placeholder rows. Combined modified-C statistics are relative to the supported model/encoding and must not imply that unmeasured 5hmC was tested and absent.

Routing rules:

| Condition                                                                          | Methylation behavior                                                                                                                 |
| ---------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| `qc_methylation=false`                                                             | Do not run pileup; component `disabled/user_disabled`                                                                                |
| No supported usable modification data                                              | Do not run pileup; `unavailable/no_usable_modifications`                                                                             |
| Some untagged or legal non-probabilistic records, with usable records elsewhere    | Run pileup on the original BAM; mark coverage of source data partial and report fractions excluded from usable modification evidence |
| Malformed/inconsistent modification encodings on eligible primary records          | Do not run pileup; `unavailable/invalid_modification_encoding`; retain complete alignment QC and factual integrity counters          |
| Usable tags but no valid CpG calls after filtering                                 | Successful empty result; `complete` with zero observed sites, null beta, and `no_valid_cpg_calls` warning                            |
| Unsupported cytosine modification codes that prevent valid combined interpretation | Withhold affected methylation metrics with `unsupported_cytosine_modification`; do not silently count them as 5mC/5hmC               |
| Corrupt BAM or a tool execution failure                                            | Fail visibly; do not turn it into a normal missing-data result                                                                       |

A malformed bedMethyl produced by a supposedly successful pileup is an execution/contract failure, not a biological empty sample. Raw parser diagnostics and read identifiers stay in private task logs; public fixture assertions and aggregate warnings use counts and stable codes.

### 7.3 bedMethyl aggregation

Adapt F2's streaming grouped-row approach. Require the pinned extended 18-column bedMethyl layout, validate contig/coordinate ordering, and group contiguous rows by `(contig, forward_C_position)`. Require strand-combined records, with one supported row per modification code. Reject duplicate same-code rows and inconsistent shared valid/canonical counts.

Use columns according to the pinned Modkit specification, not the BED score column or percent-modified column:

```text
valid = column 10 (1-based)
modified_for_row = column 12
canonical = column 13
other_modification = column 14
deletion = column 15
filter_fail = column 16
differing_base = column 17
no_call = column 18
```

Validate integer, nonnegative counts and the expected relationship between valid, canonical, modified, and other-modification counts. For supported consistent m/h records, combined modified count is `valid - canonical`; do not sum two copies of `valid`. If both m/h rows are present, their modified counts must be consistent with the shared denominator and other-modification fields. Reject unsupported interpretations instead of manufacturing values.

For each supported modification and scope, accumulate:

- Reference CpG dyads, covered CpG dyads with valid depth >0, and covered/reference fraction.
- Valid-call depth histogram; site counts at valid depth ≥1/5/10/20; corresponding fractions of reference CpGs.
- Total valid and modified calls.
- Site-average beta: `sum(modified_i/valid_i)/number_of_valid_sites`.
- Depth-weighted beta: `sum(modified_i)/sum(valid_i)`.
- 100 beta bins `[0,.01), [.01,.02), ..., [.99,1]`; beta exactly 1 belongs to the final bin.
- Separate beta histograms for all valid sites and those with valid depth ≥10.
- Counts of filtering/no-call/deletion/different-base categories, counted once per dyad when shared across m/h rows. For `filter_fail_fraction`, use `filter_fail/(valid+filter_fail)`; for `no_call_fraction`, use `no_call/(valid+filter_fail+no_call)`. Labels must retain these exact denominators and must not imply a fraction of all aligned reads.

Required scopes: `genome_primary`, `other_contigs`, `enrichment`, `targets`, `off_enrichment`. Target/enrichment aggregates are restricted to primary chromosomes for comparability with coverage; original nonprimary targets remain represented separately. Assign CpGs to a scope using the forward-strand C coordinate, not whether either base of the dyad intersects the BED. This intentionally differs from read-overlap semantics.

Scan reference sequence in bounded chunks with a one-base overlap so a CpG spanning a chunk boundary is counted once. Do not infer the reference denominator from observed pileup rows. Deduplicate CpGs in merged-region scopes, retain zero-covered denominators, and count no reference Ns as CpGs.

### 7.4 Optional classifier-feature summaries

For each supplied panel, intersect its unique forward-C features with the new full-BAM methylation summary. Report counts and panel fractions at valid depths ≥1/5/10/20, zero-observed features, and valid-depth median including zeros. Keep panel model name/version, assembly, and BED checksums in provenance.

Do not treat these counts as a classifier's minimum feature requirement or predict classifier accuracy. Full-BAM feature availability is not the same as feature availability in native Classy's possibly subsampled input. No scores, confidence thresholds, classification calls, or model preprocessing are changed.

If methylation is unavailable/disabled, panel reference size may still be known, but observed feature coverage remains null with the methylation reason code. Do not emit zero observed coverage merely because pileup did not run.

## 8. Preprocessing metrics and safe technical provenance

### 8.1 Stage accounting

Use stage IDs `prepared_input`, `demultiplexed_assigned`, `demultiplexed_unclassified`, `post_trim`, `quality_pass`, `quality_fail`, `mapping_input`, and `mapping_output`. Every row records `scope_type` (`sample` or `run`), workflow scope ID, format, stage, completeness, source kind, count, and unit.

Prepared-input BAM means the existing merged basecalled material, not original POD5 yield. Multiplex prepared input is run-scoped. Count reads/bases in existing unclassified output even though native mapping excludes it. Post-trim counts are sample-scoped. Reuse SeqKit totals for matching quality-pass/fail FASTQ stages. Count unique existing input files once per stage; do not add aggregate and component statistics for the same material.

Read retention ratios may be calculated only between compatible successive scopes/stages with exhaustive counts and explicit denominators. Label trimming base retention separately from read retention. A zero denominator is null. Run-wide counts must not be attributed independently to every sample. Do not infer barcode cross-contamination from unclassified yield.

For imported aligned BAM, only `mapping_output` is guaranteed. For imported uBAM, original input counts can be obtained by an added streaming sidecar, but previous trimming/basecalling yield is unknown. Imported FASTQ can use existing input SeqKit results when available. Missing stages produce explicit unavailable rows, not fabricated counts.

If a source is bounded, record `exhaustive=false`, `records_examined`, and method. Do not scale a sample to total yield or label it an exact total. F7's bounded checker is a semantic reference, not an exhaustive counter to port unchanged.

### 8.2 Provenance allowlist

Publish selected technical values: pipeline revision, tool versions, immutable image references, assembly, reference/region checksums, coordinate/filter definitions, basecalling and modification model names when unambiguously recorded, and aggregate missing/mixed-model flags.

Do not publish raw `@RG`/`@PG`/`@CO` records, command lines, sample/library header text, full source URIs, source filenames, read names, instrument run IDs, flow-cell IDs, access parameters, or arbitrary logs. A model path may contain private information; only retain a validated model basename matching documented model-name syntax, otherwise emit unknown. The workflow sample ID is the intended runtime association; synthetic IDs alone appear in this repository.

Do not conflate a requested model selector with a resolved model identifier. Preserve mixed-model status rather than selecting the first header record. Header settings attach to their original processing stage; basecall trimming settings do not describe a later trimming task.

## 9. Output data model

### 9.1 Versioning and component states

Create a new schema identity `oncoseq_extended_qc`, version `1.0.0`; do not masquerade as FrankONTstein's schema 2. Use JSON Schema Draft 2020-12. Require these top-level fields and reject unexpected top-level properties:

| Field            | Type / constraint                                                                                |
| ---------------- | ------------------------------------------------------------------------------------------------ |
| `schema_name`    | Constant `oncoseq_extended_qc`                                                                   |
| `schema_version` | Constant `1.0.0`                                                                                 |
| `sample_id`      | Validated workflow sample string                                                                 |
| `input_scope`    | Constant `full_mapping_stage_sample`                                                             |
| `status`         | `complete` or `partial` for an assembled sample                                                  |
| `components`     | Object containing `alignment`, `coverage`, `methylation`, `classifier_features`, `preprocessing` |
| `provenance`     | Technical allowlist object described above                                                       |
| `definitions`    | Metric registry entries keyed by stable metric ID                                                |
| `metrics`        | Array of scalar metric records                                                                   |
| `histograms`     | Array of histogram records                                                                       |
| `warnings`       | Array of `{code, component, count, message}`; no raw record data                                 |
| `tables`         | Array of published table basenames and schema versions; relative basenames only                  |

Each component has `status`, `reason`, and `exhaustive`. Allowed component statuses are `complete`, `partial`, `unavailable`, `disabled`, and `not_applicable`; `reason` is null for normal completion, otherwise a stable machine-readable code. Lack of an optional panel does not make overall sample status partial. Expected missing preprocessing history also does not degrade overall status. Overall status is partial when an enabled measurable component is unavailable/partial because of data limitations, such as missing modification evidence or incomplete known query lengths. Counts can be exhaustive even when some derived measurements are unavailable; state this at metric level.

If effective QC is disabled, there is no sample artifact. Do not publish stale-looking `disabled` files from a run that never scheduled the branch. A disabled methylation component within an enabled branch does receive an explicit status.

### 9.2 Scalar record contract

Use the same long-form concept in JSON and TSV. Every scalar JSON record has the following keys; absent dimensions are null:

```json
{
  "metric_id": "coverage.mean_depth",
  "scope": "targets",
  "mapq": 20,
  "target_id": null,
  "panel_id": null,
  "stage": null,
  "modification": null,
  "threshold": null,
  "unit": "aligned_bases_per_reference_base",
  "value": 2.0,
  "numerator": 10,
  "denominator": 5,
  "status": "complete",
  "reason": null
}
```

Counts use integer values. Measurements use finite JSON numbers. `numerator` and `denominator` are null for metrics without a ratio definition; counts do not get a made-up denominator of one. Null values require a non-null reason. A valid observed count of zero has status complete.

The registry must define each emitted metric's description, unit, eligible-record policy, numerator, denominator, coordinate scope, and null conditions. Use stable namespaces: `alignment.*`, `read_length.*`, `read_quality.*`, `coverage.*`, `allocation.*`, `adequacy.*`, `methylation.*`, `classifier_features.*`, and `preprocessing.*`. Once a field ships, changing its denominator requires a new metric ID or schema-major change.

Do not place string provenance or warning text in the numeric metric `value` field. Enumerated nonnumeric context belongs in dimensions/provenance/status fields.

### 9.3 Histogram contract

Each histogram has `metric_id`, the scalar dimensions above, `unit`, `binning`, `denominator`, `status`, `reason`, and `bins`. Each bin has `lower`, `upper`, `upper_inclusive`, `count`, and optional `base_sum`.

- Length histograms use 100-bp bins, `denominator=known_length_records`; `base_sum` stores actual sequenced bases in that bin.
- MAPQ and depth histograms use exact discrete bins, with lower=upper and upper-inclusive true.
- Beta histograms use width 0.01 with the last upper boundary inclusive.
- Empty observed histograms have empty bins and denominator zero. Unavailable histograms have null denominator and an explicit reason.
- Exact length-frequency counters remain internal; do not publish per-read observations.

### 9.4 Table shapes

All TSVs use UTF-8, tabs, a single header row, deterministic ordering, LF endings, no thousands separators, and `NA` for null. Quote/escape fields with a real TSV writer. Reject control characters in user-provided labels that would create misleading rows. A supplied-but-unavailable optional component still gets its table with schema header and explicit unavailable metric rows where a denominator is known. No supplied panels yields a header-only `classifier_cpgs.tsv` plus JSON status.

Common scalar suffix columns are, in order:

```text
metric_id, unit, value, numerator, denominator, status, reason
```

The following dimension prefixes precede that suffix:

| File                  | Dimension prefix / content                                                                                                                           |
| --------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| `summary.tsv`         | `sample_id, input_scope, scope, mapq, target_id, panel_id, stage, modification, threshold`; scalar metrics only, excluding per-target/per-panel rows |
| `read_lengths.tsv`    | `sample_id, scope`; exact group summaries and ≤500/≤1000/≥10000 fractions; detailed bins remain in JSON                                              |
| `coverage.tsv`        | `sample_id, scope, contig, mapq, threshold`; aggregate/per-contig scalar coverage; detailed depth histograms remain in JSON                          |
| `targets.tsv`         | `sample_id, target_id, label, contig, start, end, mapq, threshold`; original target rows, one record per metric                                      |
| `methylation.tsv`     | `sample_id, scope, modification, threshold`; aggregate scalar methylation values; beta/depth distributions remain in JSON                            |
| `classifier_cpgs.tsv` | `sample_id, panel_id, model_name, model_version, threshold`; unique-feature summaries                                                                |
| `adequacy.tsv`        | `sample_id, scope, target_id, mapq, assumed_vaf, min_support, probability_target, required_depth`; theoretical breadth metrics                       |
| `preprocessing.tsv`   | `scope_type, scope_id, stage, source_kind, exhaustive, records_examined`; stage counts/retention                                                     |

Add target coordinate/label metadata to the JSON provenance's `target_intervals` array so every `target_id` in JSON resolves without reading a separate table. This is panel-level reference metadata, not per-read or per-CpG results. Histograms and large per-target tables need not be duplicated in `summary.tsv`.

Row ordering follows reference contig order, then original target line number or sorted panel ID, then MAPQ, metric ID, and numeric threshold. JSON arrays use the same stable ordering. Floating-point equality tests use defined tolerances; do not round stored measurements for presentation.

## 10. Publication and retention contract

### 10.1 Allowlist

Only `QC_EXT_ASSEMBLE` publishes per-sample files, at:

```text
<outdir>/qc_extended/<sample_id>/metrics.json
<outdir>/qc_extended/<sample_id>/summary.tsv
<outdir>/qc_extended/<sample_id>/read_lengths.tsv
<outdir>/qc_extended/<sample_id>/coverage.tsv
<outdir>/qc_extended/<sample_id>/targets.tsv
<outdir>/qc_extended/<sample_id>/methylation.tsv
<outdir>/qc_extended/<sample_id>/classifier_cpgs.tsv
<outdir>/qc_extended/<sample_id>/adequacy.tsv
<outdir>/qc_extended/<sample_id>/preprocessing.tsv
```

`QC_EXT_RUN_SUMMARY` publishes only `qc_extended/run_summary.tsv`, the union of sample scalar summary rows in stable sample order. A dedicated `QC_EXT_RUN_PREPROCESSING` publisher may publish only `qc_extended/runs/<workflow_run_key>/preprocessing.tsv` for multiplex run-scoped data. Workflow run keys come from existing workflow associations, not extraction of raw instrument identifiers.

All other `QC_EXT_*` processes have `publishDir = [enabled: false]`. Final publishers must override the generic default with exact basename allowlists, not `.json`/`.tsv` suffix acceptance. Output declarations list the expected metric files explicitly. Never emit/publish the whole working directory or a directory containing them.

A configuration sketch, adapted to actual module syntax during implementation:

```groovy
withName: 'QC_EXT_.*' {
    publishDir = [enabled: false]
}
withName: 'QC_EXT_ASSEMBLE' {
    publishDir = [
        enabled: true,
        path: { "${params.outdir}/qc_extended/${meta.id}" },
        mode: params.publish_dir_mode,
        saveAs: { filename ->
            filename in [
                'metrics.json', 'summary.tsv', 'read_lengths.tsv',
                'coverage.tsv', 'targets.tsv', 'methylation.tsv',
                'classifier_cpgs.tsv', 'adequacy.tsv', 'preprocessing.tsv'
            ] ? filename : null
        }
    ]
}
```

Test the **merged configuration and actual publication**, including regex/exact-selector precedence. Do not assume that an output declaration alone prevents the repository's broad generic publication rule from exporting intermediates.

### 10.2 Explicit exclusions and native preservation

No new publication of BAM/CRAM/SAM/FASTQ, read-level TSVs, per-CpG bedMethyl, `.tbi`/`.bai`/`.csi`, normalized BED/FASTA, per-base depth tables, HTML, images, logs, SQLite files, or scratch version files. Do not create extra split on/off-target BAMs for the calculations; interval membership can be streamed. If a future calculation requires such a file, it stays private and requires a publication regression test.

Do not change `publish_fastq`, `publish_ubam`, `publish_intermediate_bams`, existing review directories, native VCF outputs, classifier results, or pipeline-info reports. Existing native publication of a format forbidden for **new QC intermediates** is not a violation of this scope. Acceptance tests compare native outputs separately from the new QC subtree.

Work-directory intermediates remain available for `-resume`; do not add eager cleanup that breaks caching. Existing trace/report metadata may mention new task names and private runtime paths as part of normal execution; the add-on creates no separate log publication.

## 11. Resources, performance, errors, and reproducibility

### 11.1 Initial resource settings

| Process                                          | CPUs | Memory | Time | Notes                                                          |
| ------------------------------------------------ | ---- | ------ | ---- | -------------------------------------------------------------- |
| `QC_EXT_REFERENCE`                               | 2    | 8 GB   | 8 h  | Shared reference/CpG/region preparation                        |
| `QC_EXT_ALIGNMENT`                               | 8    | 32 GB  | 24 h | One read-statistics scan plus MAPQ0/20 depth streams           |
| `QC_EXT_PREPROCESSING`                           | 8    | 32 GB  | 12 h | Only material not already represented by compatible statistics |
| `QC_EXT_MODKIT_PILEUP`                           | 8    | 32 GB  | 24 h | Independent full-BAM pileup                                    |
| `QC_EXT_METHYLATION`                             | 8    | 32 GB  | 12 h | Tabix/contig aggregation                                       |
| `QC_EXT_ASSEMBLE`                                | 1    | 4 GB   | 2 h  | Small component files and validation                           |
| `QC_EXT_RUN_SUMMARY`, `QC_EXT_RUN_PREPROCESSING` | 1    | 4 GB   | 2 h  | Aggregate metrics only                                         |

Apply existing `max_cpus`, `max_memory`, `max_time`, resource limits, and bounded retry conventions. Do not let test-profile `errorStrategy='ignore'` silently swallow a required new QC failure; give new QC processes explicit failure/retry behavior with no ignore path. Retry resource/transient failures using the existing repository policy, not malformed-input/contract errors.

Every BAM-reading process and alias must be added to `conf/profiles/aws_bam.config` and covered by `tests/test_resource_coverage.py`. Preserve whole-worker reservations: low workers use 8 vCPUs up to 48 GB, larger memory requests promote to the high queue with the configured 32-vCPU reservation. Do not route QC to general CPU or GPU queues accidentally. Pure aggregation tasks need not reserve a BAM worker when they stage only metric files.

New selectors may be appended to configuration, including new immutable container selectors if needed. Do not broaden a selector so that an existing native process receives a different image, resource allocation, or publication rule.

### 11.2 Bounded execution

- At most `min(4, max(1, floor(task.cpus/2)))` simultaneous depth streams; each consists of a Python reader and Samtools subprocess. No nested unrestricted thread pools. At a one-CPU ceiling, serialize contigs through one stream.
- Modkit uses the allocated CPU count. Methylation aggregation uses at most the allocated count and available contigs.
- Length memory scales with distinct observed lengths, not number of reads. Coverage memory scales with distinct depths and intervals. Report these bounds honestly; they are not mathematically constant memory.
- Exact per-target depth histograms can be expensive for very large panels. Retain per-target sufficient statistics, threshold counts, and a compact exact histogram only where needed for median/uniformity; release contig-local structures after aggregation.
- Avoid unnecessary duplicate full-BAM staging. Share existing immutable input paths; do not materialize on/off-target BAMs or split BAMs per contig. New methylation work follows the alignment integrity result rather than duplicating a preliminary tag scan.
- Reference CpG scanning/checksums are reference-keyed and reusable; panel intersections can be recomputed without rescanning all BAMs.
- A run-level summary waits for small per-sample metrics, not BAMs. Staging filenames must include sample/component directories so repeated `metrics.json` basenames do not collide.
- Record wall time, peak resident memory, input size, and temporary disk footprint during private performance validation. No large clinical input or performance log is committed.

Default-on adds substantial I/O and an independent methylation pileup. Document this cost clearly. Do not eliminate it by changing Classy behavior or sacrificing denominators. If measured runtime violates operational expectations, report the evidence and design a separate validated optimization.

### 11.3 Error contract

| Failure or absence                                                            | Required behavior                                                                   |
| ----------------------------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| Missing exact target, invalid BED, incompatible reference/panel               | Fail validation before expensive QC calculations                                    |
| Duplicate/unknown sample or missing mandatory component                       | Fail with an explicit contract error; never silently drop a sample                  |
| Expected optional input absent                                                | Emit unavailable/disabled component status and finish required metrics              |
| No reads/no mapped reads/no coverage                                          | Valid empty observations with zero counts, null undefined ratios, factual warning   |
| Invalid modification encoding                                                 | Preserve alignment metrics; withhold affected methylation as specified in section 7 |
| Truncated/corrupt BAM, failed subprocess, malformed tool output               | Task failure; no successful final sample publication                                |
| JSON/schema inconsistency, negative denominator correction, impossible counts | Task failure                                                                        |
| Partial preprocessing history                                                 | Explicit unavailable stages; no invented input yield                                |

Do not write broad `try/except` fallbacks returning zero, use `errorStrategy='ignore'`, or silently omit missing files. The final assembler requires all mandatory completion records before emitting any publishable sample files. Atomic file writes inside that task avoid accepting a truncated metric file after a crash.

Factual warnings include `no_mapped_primary_reads`, `no_target_coverage`, `no_valid_cpg_calls`, `targets_outside_enrichment`, `mapq20_target_dropout`, `mixed_basecalling_models`, and `potential_methylation_depth_cap`. Do not port FrankONTstein's illustrative 95%-at-10x warning as a clinical threshold. No pass/fail designation is based on an unvalidated biological cutoff.

### 11.4 Determinism and compatibility

Sort rows and dimensions explicitly. Avoid timestamps, absolute work paths, and random UUIDs in metric payloads; execution timestamps already belong to native provenance. Pin containers and record actual versions. Record threshold-sampling seed/settings and tested limitations of tool determinism. Test serial versus parallel aggregation equality with exact integer equality and a relative/absolute tolerance of `1e-12` for derived floating-point values.

Do not promise that an image pin alone produces identical threshold estimates across all hardware/thread settings. Acceptance requires repeated identical-input runs in the supported Linux environment to agree on counts and within documented tolerances on derived values. Unexpected differences must be explained before release.

Native resume compatibility is established by trace/task-hash comparison, not a statement that any workflow edit is automatically cache-neutral. Purely added QC inputs must not become inputs of native processes. If an existing native task's cache key changes, identify the exact cause and remove avoidable coupling.

## 12. Ordered implementation packages for later agents

These packages are dependency ordered specifications. Create/update Beads issues during the later implementation session; do not use this document as a mutable checklist. A later agent must read sections 1–11 before editing code and inspect the current repository revision for intervening changes.

### Package A: contracts, fixtures, and reference validation

**Depends on:** none.

**Edits:** new Python package/CLI, output schema, fixture builder, reference helper, new tests. Parameter wiring is deferred to package F so default-on behavior is not released before working metrics exist.

Implement the long-form scalar/histogram serializers and null handling first. Add a metric registry with the denominator definitions above, reference/BED validation, interval unions, stable original target IDs, primary-contig selection, CpG denominator scanning, and panel validation. Keep functions independent of Nextflow.

Create `tests/qc_extended/build_fixture.py --outdir DIR`. Generate synthetic FASTA, index, sorted BAM/index, exact-target and enrichment BEDs, tagged/untagged/bad-tag variants, panel manifests, and tiny valid bedMethyl fixtures in that directory. No downloaded data or existing clinical output is used.

**Completion evidence:** unit tests validate schemas, paths, interval/CpG boundaries, status behavior, and fixture expectations. Source attribution is present. A real-container capability probe records required tools/options without claiming a stub proves compatibility.

### Package B: alignment, read-length, and coverage restoration

**Depends on:** A.

**Edits:** `alignment.py`, regions/contract helpers, `extended_qc.py alignment`, tests; new process wrapper may be introduced without top-level wiring.

Implement one-pass primary-record accounting, exact length groups, mapped blocks, quality/clipping counts, modification-integrity completion status, MAPQ0/20 depth streams, and per-target coverage. Use F2 algorithms selectively; do not import reporting/provenance code wholesale. Keep the original target intervals alongside union indexes.

**Completion evidence:** F5-derived numerical fixture, independent brute-force tiny-reference depth oracle, flag/CIGAR cases, all-primary count reconciliation, one-versus-multiple worker equivalence, missing-length behavior, and bounded-memory review. Prove no per-read or per-base output file is written.

### Package C: independent methylation QC

**Depends on:** A and the integrity output from B.

**Edits:** new pileup/index wrapper, `methylation.py`, explicit tool-version parser, tests, optional-panel aggregation.

Implement the new private modkit pileup, index, CpG aggregation, missing-modification semantics, reference denominators, feature coverage, and completion records. Use the supplied reference, not an image-bundled assembly. Retain separate native Classy work.

**Completion evidence:** pure aggregation oracles plus a real pinned-container modkit run using synthetic MM/ML data. Verify supported options, strand combination, m-only/m+h availability, tool flag behavior, threshold metadata, empty output, and no published pileup/index/log. Resolve any source-versus-actual-tool discrepancy before proceeding.

### Package D: preprocessing and technical provenance

**Depends on:** A.

**Edits:** preprocessing Python code, additive emits in simplex/multiplex subworkflows, new sidecar modules and tests. Do not modify existing tool process bodies.

Expose matching pass/fail statistics; count prepared/demux/unclassified/post-trim material only where necessary. Use explicit sample/run/barcode associations. Implement the provenance allowlist, mixed-model flags, and partial-history semantics. Keep run-scoped summaries separate from sample rows.

**Completion evidence:** simplex and multiplex synthetic stage tables, exact run/sample accounting, retained unclassified counts, no duplicated input files, a sample name containing `_pass`/`_fail` that remains intact, and fixtures proving raw headers/read names/source URIs do not enter aggregate payloads.

### Package E: novel descriptive metrics

**Depends on:** B; classifier-feature summaries also require C.

**Edits:** allocation and adequacy calculations in the Python metric engine, stable registry entries, tests.

Implement target-base efficiency, padding burden, nested short/long length fractions, MAPQ attrition, theoretical binomial adequacy, and optional classifier features. Use existing aggregated counters rather than adding a fresh BAM scan for each metric. Do not implement the deferred research extensions in section 14.

**Completion evidence:** conservation checks for coordinate partitions, all eight binomial depth oracles, target unions with overlap, unavailable partial query-base denominators, and feature-panel deduplication. No metric is named clinical sensitivity, tumor fraction, or adaptive rejection rate.

### Package F: Nextflow integration, publication, resources

**Depends on:** A–E.

**Edits:** new subworkflow/modules; additive wiring in `workflows/oncoseq.nf`; typed/default/schema parameters; new selectors in module/container/resource configuration; tests. Existing selectors retain their values.

Implement sample registry validation and explicit completion records. Add all publication overrides before enabling the workflow. Route BAM readers through isolated AWS workers and clamp worker concurrency. Emit tool versions through the existing version mechanism without separately publishing scratch YAML/logs.

Add `test_qc_extended` to `nextflow.config` and its synthetic profile. Supply explicit synthetic exact targets in the QC-enabled acceptance profile. Existing test profiles without a verified exact-target asset should set `qc_extended=false` explicitly to preserve their native-path regression role; never relabel the eroded enrichment BED as exact to make a test pass. Add focused tests proving the actual global default remains true.

For fast full-entrypoint wiring validation, use the existing native stub capability where available and the new synthetic profile; mark that evidence as wiring only. Real metrics are separately required through the real-container QC subworkflow acceptance harness. Do not invent new caller-skip flags or modify native tools to make a tiny synthetic full analysis succeed.

**Completion evidence:** all expected sample completions, independent report/classifier enablement, native time-series full-sample scope, cfDNA bypass, exact allowlist publication, actual merged-profile resource routing, and resume/native-command comparisons.

### Package G: acceptance, documentation, handoff

**Depends on:** F.

**Edits:** user usage/output documentation and acceptance tests; fixes remain confined to the additive feature unless separately authorized.

Document required exact targets, intentional default-on behavior, disable options, outputs, status codes, full-BAM scope, compute cost, optional panel requirements, and interpretive limitations. Add a concise link from README to the usage documentation. Preserve all existing report documentation.

Run affected tests and at least one `nextflow run` profile exercising the new path, plus the real QC container/subworkflow harness. Record exact commands, versions, results, and unavailable infrastructure. A remote/private full-size performance run is a release validation item, not permission to copy clinical artifacts into the repo.

**Completion evidence:** every acceptance category below has results or a precisely identified external blocker. Do not label missing runtime tests as passed. Future commits/pushes must follow the repository's public-content/Gitleaks/outgoing-history requirements; the documentation-only action that created this file performs none of those mutations.

## 13. Tests, synthetic oracles, and validation commands

### 13.1 Core parity fixture

Recreate F5's synthetic 100-bp `chr1` reference as `CG` repeated 50 times. Enrichment intervals are `[0,10)` and `[5,15)`; exact target is `[5,10)`.

Write these records, coordinate-sort, and index:

| Synthetic record     | Start    | Query length | MAPQ | Flag             | CIGAR      |
| -------------------- | -------- | ------------ | ---- | ---------------- | ---------- |
| `SYNTHETIC_READ_A`   | 0        | 10           | 60   | primary mapped   | `10M`      |
| `SYNTHETIC_READ_B`   | 10       | 10           | 0    | primary mapped   | `10M`      |
| `SYNTHETIC_SUPP`     | 10       | 10           | 60   | supplementary    | `10M`      |
| `SYNTHETIC_READ_C`   | 20       | 20           | 60   | primary mapped   | `10M5D10M` |
| `SYNTHETIC_UNMAPPED` | unplaced | 5            | 0    | unmapped primary | none       |

Expected values:

| Metric                                      | Expected    |
| ------------------------------------------- | ----------- |
| Primary / mapped / unmapped records         | 4 / 3 / 1   |
| Alignment rate                              | 0.75        |
| All-primary known query bases               | 45          |
| All-primary median / N50                    | 10 / 10     |
| On-target / off-target mapped reads         | 1 / 2       |
| On-enrichment / off-enrichment mapped reads | 2 / 1       |
| On-enrichment fraction among mapped reads   | 2/3         |
| Enrichment union length                     | 15          |
| Genome aligned bases MAPQ0 / MAPQ20         | 40 / 30     |
| Enrichment aligned bases MAPQ0 / MAPQ20     | 15 / 10     |
| Target aligned bases MAPQ0 / MAPQ20         | 5 / 5       |
| Genome breadth ≥1 MAPQ0 / MAPQ20            | 0.4 / 0.3   |
| Target mean depth MAPQ0 / MAPQ20            | 1 / 1       |
| Query-base target yield MAPQ0 / MAPQ20      | 5/45 / 5/45 |
| Reference primary-genome CpGs               | 50          |

Use an independent tiny Python brute-force CIGAR/reference oracle in tests, not a second invocation of the implementation under test. Samtools output is an additional cross-check. The supplementary alignment must not inflate primary coverage.

For a CpG at coordinate 0, create matching m/h rows with valid depth 10, canonical count 5, modified m=3 and h=2, and reciprocal other-mod counts. Expected one dyad, valid depth 10 once, combined beta 0.5, 5mC beta 0.3, 5hmC beta 0.2, primary reference denominator 50, and target observed CpG sites zero because coordinate 0 is outside `[5,10)`.

The example scalar record in section 9 is illustrative schema data and is not the oracle for this fixture.

### 13.2 Required behavioral cases

| Category              | Cases and assertions                                                                                                                                      |
| --------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Interval unions       | Overlap, adjacency, duplicate rows, duplicate labels, retained original line IDs, zero-covered targets; no union double count                             |
| Interval validation   | Negative/reversed/out-of-bounds coordinates, unknown contigs, wrong reference, empty files; fail with defined error                                       |
| Overlap semantics     | `10M10D10M` and `10M10N10M` over a target in the gap remain off target; insertion/soft-clip-only overlap does not count                                   |
| Read length           | Even median, N50 base weighting, hard clipping, absent sequence with/without inferable CIGAR, very long read, empty group                                 |
| Alignment flags       | Secondary, supplementary, both flags, duplicate/QC-failed primary, placed and unplaced unmapped; exact reconciliation                                     |
| Mapping quality       | 0,19,20,60,255; numerical threshold behavior and unknown-MAPQ count                                                                                       |
| Coverage zeros        | Empty BAM, no mapped reads, no target coverage, complete enrichment universe, zero off-enrichment denominator                                             |
| Contig universes      | UCSC, bare, RefSeq primary, mitochondrial/alternate contigs, nonprimary target retained in per-target output                                              |
| Quality metrics       | Valid/missing/nonfinite `qs`, complete/missing base qualities, known error-probability conversion                                                         |
| Modification tags     | m-only, m+h, tag-free, partial tag coverage, valid empty MM, missing optional MN, missing probability arrays, malformed deltas/ML length, incompatible MN |
| bedMethyl             | Duplicate same-code dyad row, inconsistent denominators, negative/impossible counts, unsorted rows, empty valid output, unknown modification code         |
| CpG reference         | Chunk boundary, CpG at interval edge, overlapping target intervals, reference Ns, zero reference CpGs                                                     |
| Methylation fractions | Site-average differs from depth-weighted example; shared m/h valid count once; absent h remains unavailable                                               |
| Feature panels        | Duplicate features, shared features across panels, wrong assembly, non-CpG site, absent panel, unavailable methylation                                    |
| Sampling adequacy     | All eight thresholds in section 6.5, n below k, zero depth, VAF/support boundaries, overlapping targets not double-counted                                |
| Preprocessing         | Run/sample scopes, unclassified yield, matching reused SeqKit stats, partial imported history, suffix-containing sample IDs                               |
| Optional branches     | Two samples with different modification availability; both finish with complete mandatory metrics                                                         |
| Failure propagation   | Nonzero Samtools/modkit exit, corrupt BAM, missing mandatory metric, malformed JSON; final publisher does not claim completion                            |
| Serialization         | JSON Schema validation, no NaN/infinity, TSV `NA`, stable order, exact units/denominators, no raw record/header leakage                                   |

Add metamorphic checks: reordering BED rows changes original-row IDs but not union metrics; adding secondary/supplementary records changes their counters but not primary metrics; doubling independent synthetic primary records doubles coverage and record counts but preserves length-group fractions; changing only exact targets changes target metrics but not alignment totals; worker-count changes preserve integer results.

### 13.3 Publication and native-preservation tests

Use actual Nextflow publication behavior, following the style of `tests/modules/local/review/publication.nf.test`. Create a synthetic task whose work directory also contains decoy BAM, bedMethyl, index, log, HTML, PNG, and unapproved JSON/TSV files. Assert all intended metric basenames are published and every decoy is absent. Test merged generic plus specific configuration and `copy` and `symlink` publication modes.

Assert intermediate files remain available in task work directories and are usable by downstream QC tasks. Do not equate unpublished with absent from work storage.

Test two samples publishing identical basenames under different sample directories, plus two run-scoped preprocessing tables. Assert no overwrites, basename collisions, or orphaned samples. Header-only optional tables are still valid allowed metric files.

For native preservation, compare QC-off and QC-on synthetic runs with identical native inputs/settings:

- Existing native process command scripts and relevant resolved configuration are identical.
- Existing outputs are retained. Compare semantic contents for files with unavoidable timestamps; compare bytes where stable.
- Extra files are confined to the new QC subtree, apart from ordinary native trace/version metadata describing added tasks.
- A QC-on `-resume` run after QC-off reuses native tasks where their hashes are unchanged and adds only necessary QC work.
- Existing standalone alignment QC remains gated by `report_review`; the new metrics exist when review is off.
- Time-series outputs retain their native naming/routing, while new QC has one full-sample entry per original sample.
- cfDNA runs bypass the add-on and do not newly require exact targets.

### 13.4 Test harness and command templates

Create the named fixture builder, container smoke driver, subworkflow acceptance entrypoint, and profile as part of implementation. These commands describe the **future** test interfaces; they are not claimed to run before those files exist.

Use a Python environment with pysam and jsonschema and the pinned production-compatible Samtools for unit tests. Prefer the pinned QC container for tool-dependent calculations. Do not install arbitrary packages into an existing user's Python environment.

```bash
# Pure metric/schema tests; tool-dependent cases must not silently skip in CI.
python3 -m unittest discover -s tests/qc_extended -p 'test_*.py' -v

# Real Linux modkit plus Python/Samtools smoke test, all inputs synthetic.
python3 tests/qc_extended/run_container_smoke.py --threads 4

# New process and channel contracts, including actual publication.
nf-test test tests/modules/local/extended_qc --profile docker --verbose
nf-test test tests/subworkflows/local/extended_qc --profile docker --verbose
nf-test test workflows/tests/extended_qc.nf.test --profile docker --verbose

# Existing configuration/resource contracts and alignment-QC regression.
python3 -m unittest discover -s tests -p 'test_resource_coverage.py' -v
python3 -m unittest discover -s tests -p 'test_typed_configuration.py' -v
python3 -m unittest discover -s tests/review -p 'test_alignment_qc.py' -v

# Inspect merged configuration without launching tasks.
NXF_OFFLINE=true nextflow config -profile test_qc_extended,docker

# Pipeline lint after implementation.
nf-core pipelines lint .
```

The future `run_container_smoke.py` must use temporary directories, readonly source mounts, synthetic data, and explicit immutable images. It must assert actual modkit-produced counts and not just file existence. Container networking can be disabled after images are available. It must exit nonzero when tools/options are absent rather than reporting a skipped success.

Provide a real subworkflow acceptance entrypoint at `tests/subworkflows/local/extended_qc/acceptance.nf`, plus a fixture builder producing `inputs.json` in the generated fixture directory. `inputs.json` contains only synthetic sample/reference/region paths and settings required by that harness. The test profile resolves fixture paths from the environment variable `ONCOSEQ_QC_TEST_FIXTURE_DIR`; do not add a production parameter for the test fixture directory. Validate that the variable names a complete synthetic fixture directory before launch.

```bash
qc_acceptance_root=$(mktemp -d "${TMPDIR:-/tmp}/oncoseq-qc-acceptance.XXXXXX")
python3 tests/qc_extended/build_fixture.py --outdir "$qc_acceptance_root/fixture"

# Real additive subworkflow: actual coverage and modkit, no native caller stubs.
nextflow run tests/subworkflows/local/extended_qc/acceptance.nf \
    -profile docker \
    -params-file "$qc_acceptance_root/fixture/inputs.json" \
    -work-dir "$qc_acceptance_root/work-real" \
    --outdir "$qc_acceptance_root/results-real"

# Full pipeline entrypoint wiring through a synthetic test profile.
# Stub evidence is additional to, not a replacement for, the real test above.
ONCOSEQ_QC_TEST_FIXTURE_DIR="$qc_acceptance_root/fixture" \
nextflow run . -profile test_qc_extended,docker -stub-run \
    -work-dir "$qc_acceptance_root/work-wiring" \
    --outdir "$qc_acceptance_root/results-wiring"
```

The pipeline has typed parameters. Keep parameter validation enabled in the test profile and declare production QC parameters in all three required locations. The environment-variable test-fixture interface above avoids adding an undeclared test-only command-line parameter.

Full production-route biological validation still requires realistic private inputs and reference assets. The tiny fixture proves mathematics and wiring, not clinical equivalence or throughput. Record runtime/tool versions and any test failures separately from pre-existing suite failures.

### 13.5 Release acceptance

Acceptance requires all of the following evidence in the later implementation handoff:

1. Source/parameter diff review shows additive scope and preserves native process behavior.
2. Every enabled synthetic sample produces valid required metrics and explicit optional-component status.
3. Numerical fixtures and independent overlap/depth/CpG oracles pass.
4. Actual pinned-container methylation execution passes; stub-only evidence is insufficient.
5. Actual Nextflow publication contains only allowed new metrics and no QC-only intermediates.
6. Existing output retention, time-series scope, cfDNA bypass, and resume behavior are tested.
7. New BAM readers have isolated-worker coverage, and measured resource use is recorded.
8. Tests, lint, usage documentation, and PHI/credential review are complete, with any infrastructure blockers stated honestly.

No automatic clinical pass/fail threshold or decision-rate claim is a release condition or output.

## 14. Additional useful QC proposals and deferred work

“Novel” here means useful beyond the inspected workflow, not a claim of scientific novelty.

| Proposal                                        | Why useful for cellular liquid-tumor adaptive sampling                                                    | Status / prerequisite                                                                                                                                         | Interpretation limit                                                                                     |
| ----------------------------------------------- | --------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------- |
| Target-base efficiency and padding burden       | Separates long-read target overlap from actual bases allocated to intended loci                           | Implement now, section 6.4                                                                                                                                    | Not enrichment versus a nonadaptive control                                                              |
| Short/long off-enrichment burden                | Reveals changes in observed background length composition                                                 | Implement now                                                                                                                                                 | Cannot label individual reads as accepted/rejected                                                       |
| MAPQ attrition and target dropout               | Highlights loci where coverage depends on ambiguous alignments, relevant to paralogous/homologous regions | Implement now                                                                                                                                                 | Additional depth may not fix ambiguity                                                                   |
| Theoretical low-VAF support adequacy            | Relates target depth distribution to assumed subclonal allele sampling                                    | Implement now, section 6.5                                                                                                                                    | Not assay sensitivity, MRD, or tumor fraction                                                            |
| Modification evidence attrition                 | Exposes loss of tags/probabilities or callable CpGs despite apparently adequate sequence depth            | Implement now                                                                                                                                                 | Global methylation changes may be biological; no universal cutoff                                        |
| Classifier-feature coverage                     | Measures whether unique relevant CpGs are observed rather than repeatedly oversampling a small subset     | Implement when a matching panel is supplied                                                                                                                   | Does not predict score validity or native classifier input coverage                                      |
| Off-enrichment CNA background diagnostics       | Background coverage may support CNAs, but target halos, GC, and mappability can distort it                | Defer: define versioned bin/mappability/exclusion resources and a reference cohort; begin with depth/zero-bin/GC summaries and robust adjacent-bin dispersion | Biological CNAs and sex chromosomes also change dispersion; do not call raw MAD a calibrated noise score |
| Target-boundary decay                           | Can reveal padding/fragment-length interactions and uneven usable coverage near edges                     | Defer: strand-aware boundary windows, overlapping-target policy, and synthetic/null validation                                                                | Coverage shape is not a direct rejection-delay measurement                                               |
| Long-read spanning adequacy                     | Potentially useful for fusion/intron/ITD interrogation beyond mean depth                                  | Defer: explicit locus/anchor manifest and contiguous aligned-query support definitions                                                                        | Uninterrupted reference spans do not prove structural-variant detection sensitivity                      |
| Temporal yield/CpG accumulation                 | Can support operational decisions about whether more sequencing adds distinct useful information          | Defer: use only preserved BAM `st`/`du` tags with completeness checks; no sequencer-export dependency                                                         | Read-start timestamps do not by themselves provide full pore occupancy or decision latency               |
| Identity/contamination checks                   | Useful for multiplex and repeated leukemia samples                                                        | Defer: validated SNP panel, ploidy/LOH/tumor-aware model, consented reference samples                                                                         | Low purity and clonal copy changes confound generic germline assumptions                                 |
| Methylation mixture or tumor-fraction inference | Could contextualize classification in mixed blood/marrow populations                                      | Defer: independently validated cell-type reference and disease context                                                                                        | Global beta/5hmC shifts alone are not tumor-fraction estimates                                           |

Do not implement sequencer decision adapters, estimated unblock latency, acceptance/rejection rates, false accept/reject rates, or cfDNA fragmentomics as part of this handoff. Those exceed the agreed scope. Optional BAM end-reason tags also must not be promoted to an experimentally validated decision metric without separate validation.

## 15. External technical references and evidence limits

Use pinned tool versions for behavior. Online manuals are explanatory references and may describe newer releases than the chosen image.

- [Samtools depth manual](https://www.htslib.org/doc/samtools-depth.html): flag inclusion/exclusion, swapped `-q`/`-Q` meanings relative to mpileup, deletion handling, zero output, and absence of the old depth cap. The inspected page described Samtools 1.25; verify the image's actual version.
- [SAM optional-tag specification](https://samtools.github.io/hts-specs/SAMtags.pdf): MM/ML/MN encoding, coordinate/strand semantics, canonical-state modes, and legal optional fields.
- [Modkit pileup and bedMethyl documentation](https://nanoporetech.github.io/modkit/intro_pileup.html): strand combination, modified-base selection, reference CpGs, and extended bedMethyl semantics.
- [Modkit confidence filtering](https://nanoporetech.github.io/modkit/filtering.html): estimated filtering thresholds versus fixed per-base/per-modification cutoffs.
- [Modkit advanced usage](https://nanoporetech.github.io/modkit/advanced_usage.html): explicit command capability verification, sampling options, and depth caps.
- [Dorado SAM specification](https://software-docs.nanoporetech.com/dorado/latest/basecaller/sam_spec/): meanings of recorded model/read-group metadata and `qs`, `st`, `du`, `dx`, and related tags. Their presence is optional for imported data.
- [Payne et al., Readfish](https://doi.org/10.1038/s41587-020-00746-x): adaptive sequencing and targeted cancer-panel context, including a leukemia fusion example. This supports the assay context, not a BAM-only decision-rate estimator.
- [Scheinin et al., QDNAseq](https://doi.org/10.1101/gr.175141.114), [open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4248318/): GC/mappability/problem-region considerations motivating deferred CNA background diagnostics. It does not validate a new adaptive-sampling QC cutoff.

During the planning/documentation session, repository source and public technical documentation were inspected. The current shell had Samtools and a Docker executable, but no discoverable Nextflow or nf-test executable and no pysam/jsonschema in its default Python. No production container, pipeline, methylation calculation, or biological equivalence test was executed. The numeric binomial depth thresholds were independently calculated with Python's standard library.

The immediate document-only completion criterion is a single added Markdown file with these contracts and no other repository changes. Later implementation acceptance is governed by section 13, not by the existence of this document.

## 16. Execution update: report integration

The user explicitly requested that QC also appear in reports. `QC_EXT_REPORT` generates private aggregate QC content. `QC_EXT_ATTACH_REPORT` adds that content to the completed native HTML, with QC navigation, so native Quarto output directories do not acquire QC-only fragments. Native HTML publication is routed through this additive finalization only when QC is enabled; final report filenames and existing report-directory outputs are retained. Review reporting consumes `qc:extended` metric artifacts and embeds the same content in the QC tab. Both display regional lengths, coverage breadth, methylation distributions, sampling adequacy, feature availability, component statuses, and the full-mapping-stage scope. Standalone QC intermediate HTML, pileups, indexes, and logs remain unpublished.

The implementation uses common long-form TSV columns across the metric tables and publishes multiplex run stage counts in `qc_extended/run_preprocessing.tsv`. These are aggregate metric-file interface refinements of sections 9 and 10. The synthetic acceptance profile is a dedicated QC entrypoint selected with `-main-script tests/subworkflows/local/extended_qc/acceptance.nf`; it validates real metric execution without altering native callers.

## 17. Locked-image command validation

Actual Linux execution showed that this modkit rejects `--seed` with `--num-reads`, and explicit `--max-depth` requires the downsampling `--high-depth` path. The production command therefore uses indexed count-based threshold sampling (`--num-reads 10042 --filter-percentile 0.1`), without a seed, explicit depth cap, or downsampling. Provenance records those settings accurately. This supersedes the seed/depth-cap assumptions in section 7; source depth near the documented optimized-path high-depth range produces a factual warning, not an automatic sampling-policy change.
