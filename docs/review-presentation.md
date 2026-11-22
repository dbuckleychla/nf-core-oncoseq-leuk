The optional review report retains Oncoseq as the primary workflow and uses the aggregator to organize native evidence. Presentation changes do not change variant calls, selection thresholds, or resource requests.

Classification cards show the top class, numeric score, emitted diagnostic cutoff, and a red cutoff marker. Input methylation thresholds remain in model evidence; they are not substituted for diagnostic cutoffs. Missing cutoffs are explicitly marked. Nanomix reports mixture proportions rather than diagnostic classifications. Missing classifier results remain visible as unavailable, and Optic remains disabled pending model availability.

Small variants show gene/context, change, caller-specific AF and depth, ClinVar assertions and review status, and optional COSMIC matches. ClinVar records retain reported conditions and conflicts, with links when a variation ID is emitted. A diagram explains aggregation and filtering. Caller support is not pooled, and catalogue membership is not a pathogenicity determination.

Copy-number views prioritize re-rendered NASVAR chromosome-arm estimates and a genome coverage/allele-fraction overview, followed by selectable chromosome figures. The overview uses native GC-adjusted binned coverage signal, not base depth or inferred copy number. Folded allele fractions require total depth of at least 20 and are averaged in 1 Mb display bins. The coverage display uses a 1st–99th percentile range extended to zero, with clipped outliers marked amber; native values are unchanged. Chromosome widths follow observed coverage-bin extents. Native numeric tables remain available in collapsed sections. Focal target-locus figures provide CN/SV context rather than independent calls.

QC summarizes available read metrics and interval-weighted coverage, with up to 20 lowest-coverage genes. Overlapping reported intervals may contribute more than once; no assay-specific pass/fail threshold is inferred. The full interval table remains collapsed. Figure dialogs provide a full-window toggle instead of slider zoom.

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
