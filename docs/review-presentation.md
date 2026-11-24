The optional review report retains Oncoseq as the primary workflow and uses the aggregator to organize native evidence. Presentation changes do not change variant calls, selection thresholds, or resource requests.

Classification cards show the top class, numeric score, emitted diagnostic cutoff, and a red cutoff marker. Input methylation thresholds remain in model evidence; they are not substituted for diagnostic cutoffs. Missing cutoffs are explicitly marked. Nanomix reports mixture proportions rather than diagnostic classifications. Marlin, ALMA, Lamprey, TUCAN, Nanomix and Optic are visible by default; other classifiers appear in a collapsed section. Missing results remain visible as unavailable, and Optic remains disabled pending model availability.

Small variants primarily use mandatory COSMIC targeted-screen annotations. Exact
normalized matches with positive blood/lymphoid tissue evidence appear in the
primary table. The collapsed expanded view includes all other COSMIC matches and
additional ClinVar/VEP candidates. Caller PASS and population eligibility apply
to both. Each allele appears once, with caller-specific percentages and read
counts. The evidence popup shows distinct COSMIC sample counts, tissue, histology
and subtypes, publications and reported somatic status, alongside ClinVar
assertions and VEP consequences. The COSMIC evidence column precedes ClinVar and
shows distinct catalogue sample counts, tissue category, gene tier and source
somatic-status categories. These use labelled colour badges; ClinVar significance
has separate formatting. Catalogue sample counts are deduplicated across COSV IDs
and callers and do not imply distinct patients. IDs remain in the evidence dialog.
Caller support shows native variant/total read counts; unavailable counts use a
dash and are never inferred from variant percentage. A single overview disclaimer
explains the limits of COSMIC evidence, without repeating it in variant reasons. Missing joins are explicit. Catalogue membership
and reported somatic status do not classify the current sample. Caller controls
filter selected rows without pooling evidence. Complete annotation exports retain
all target alleles and the joined evidence.

Copy-number views show three separate plots: NASVAR chromosome-arm copy estimates, GC-corrected read counts per bin, and B-allele frequency (BAF). The coverage plot is raw read-count coverage, not a log ratio or inferred copy number; its axis starts at zero and has count ticks. The displayed coverage range excludes extreme outliers, marked amber, while original counts remain unchanged. BAF uses alternate-allele reads divided by total reads on a labelled 0–1 axis, requiring depth of at least 20. Density bins preserve both sides of 0.5 rather than folding or averaging the fractions. Chromosome widths follow observed coverage-bin extents. All available chromosome figures appear in separate rows, with an unavailable notice for missing chromosomes. Numeric tables remain collapsed. Target-gene figures show coverage and rearrangements rather than independent calls.

QC includes [primary-read alignment rate](alignment-qc.md) and summarizes available read metrics and interval-weighted coverage, with up to 20 lowest-coverage genes. Overlapping reported intervals may contribute more than once; no assay-specific pass/fail threshold is inferred. Basic plots show gene coverage, target coverage distribution, read counts, N50 and mean read quality when the measurements are available. Read-length histograms require read-level data or precomputed distributions; summary N50 and mean length cannot reconstruct them. The current local outputs contain neither BAM/FASTQ reads nor histogram data, so the report marks these distributions unavailable. QC is the penultimate tab. The full interval table remains collapsed. Figure dialogs provide a full-window toggle instead of slider zoom.

To re-render an existing local report without rerunning analysis:

```bash
python3 bin/review_report.py render \
  --input /path/to/results/reports/review/SAMPLE/review.json \
  --figures /path/to/results/reports/review/SAMPLE/figures \
  --results-root /path/to/results \
  --output /path/to/results/reports/SAMPLE.review.updated.html
```

Use a Python environment containing the pinned review image dependencies (Jinja2, Pillow and pysam). The optional results root supplies existing NASVAR coverage/allele counts and read QC for older outputs. Updated pipeline outputs stage these artifacts directly into review data. Rendering writes small derived SVGs under the private review directory; it does not upload results or modify native analytical files.

Adding raw NASVAR and QC artifacts changes review inputs. The first resumed run after this update may rerun review preparation, VEP and completion; upstream calling and mapping remain eligible for their existing caches. Local re-rendering is preferable when only the layout needs updating.
