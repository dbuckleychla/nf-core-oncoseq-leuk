The review QC tab reports read-level alignment rate:

```text
mapped primary reads / (mapped primary reads + unmapped primary reads) × 100
```

Secondary and supplementary alignment records are excluded. QC-passed and QC-failed primary records are included. The counts come from `samtools flagstat -O json` on the mapping-stage BAM, before downstream report-specific filtering. For time-series reports, the metric describes the whole mapping-stage sample, not just the displayed time window. For imported aligned BAMs, it describes the reads retained in the supplied file.

The clinical review workflow runs a separate `SAMTOOLS_ALIGNMENT_QC` task using the existing pinned Samtools image and BAM-worker queue policy. Alignment commands and resource requests are unchanged, so existing alignments remain eligible for `-resume`. No extra scan is requested when the optional review report is disabled. The review data retains the raw counts in `qc.alignment.json`.

The old Cramino 1.1.0 `% from total alignments` metric is not substituted for this rate: its denominator includes supplementary records, while its numerator counts mapped primary reads.

An existing HTML report can receive already-computed primary-read counts without regenerating figures:

```bash
samtools flagstat -O json /path/to/aligned.bam > alignment_qc.json
python3 bin/review_report.py patch-alignment-qc \
  --input /path/to/sample.review.html \
  --output /path/to/sample.review.with-alignment-qc.html \
  --alignment-qc alignment_qc.json
```

Use the existing review Python environment. The updater changes only the embedded QC payload and preserves all other report text and figures. It can also accept normalized primary-count JSON with explicit source/scope metadata.
