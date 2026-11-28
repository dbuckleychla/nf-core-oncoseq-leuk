# Classy blood classification

The additive `CLASSY_BLOOD_WORKFLOW` runs MARLIN, Lamprey, ALMA3.1 and OPTiC using the private Classy blood runtime. The existing native Classy execution remains in place. Both executions consume the same classification BAM selected by oncoseq, including its existing adaptive-coverage/subsampling rules. This addition does not change basecalling, mapping, native variant callers or classification input selection.

## Defaults and interpretation

| Parameter                              | Default | Effect                                                                       |
| -------------------------------------- | ------- | ---------------------------------------------------------------------------- |
| `classy_blood`                         | `true`  | Run the new classifier and use its results in the native and review reports. |
| `classy_blood_enable_alma3`            | `true`  | Enable ALMA3.1's research-only native approximation.                         |
| `classy_blood_alma3_modification_mode` | `5mc`   | ALMA3 signal input; alternatively `5mc-plus-5hmc`.                           |

ALMA3.1 is explicitly labeled **native approximation** in the report. Numerical parity with the reference implementation is not established. Reports include the required attribution, **“Analysis performed using ALMA3.”** The source image retains its model licenses; its availability does not grant permission to redistribute private model artifacts.

The ALMA3 signal mode is independent of the Dorado model specification. A model that emits both 5mC and 5hmC does not automatically select the combined mode. MARLIN, Lamprey and OPTiC run when ALMA3 is disabled. Brain/solid classifiers absent from this image are excluded explicitly; they remain available in the legacy execution.

The report retains native classification scores, cutoffs, hierarchy labels and accepted-call decisions. It uses a native `accepted_call` when emitted, including decisions that cannot be represented by a single scalar cutoff. Missing model results and zero matched CpG observations are shown as unavailable. Report adaptation does not manufacture predictions or turn inadequate coverage into successful calls.

For cfDNA runs (`cfdna=true`), the new Classy step lowers its minimum read length to 20 bp; otherwise it retains the image's native 500 bp default. Read-quality and other preprocessing settings remain the image's native defaults. ALMA3 requires at least 5,000 matched release CpGs. Model/input adequacy errors can fail the native classifier process and should be investigated rather than bypassed.

## Files and reports

New final outputs publish to `reports/tumor_classifiers/`:

- `<sample>_combined_classification.json`: unmodified native combined JSON.
- `<sample>_combined_classification_combined_blood.svg`: primary blood summary.
- `<sample>_combined_classification_combined_task_pies.svg`: additional combined summary.
- `<sample>_{marlin,lamprey,alma3,optic}.tsv`: concise native score tables, or explicit unavailable/disabled rows.
- `<sample>_classification.metadata.json`: executed image, ALMA3 settings, model versions, target coverage and a plot allowlist.

Legacy Classy JSON, plots and parsed tables publish to `reports/classy_legacy/<sample>/`. The legacy execution still runs. Its outputs do not supply the native methylation section or review classification tab while `classy_blood=true`.

The native report gets a new methylation section containing the blood SVG and score tables. The review report gets the four blood-model cards, primary blood SVG and an additional-plot dropdown. Model releases and the executed image are included in provenance. The metadata allowlist prevents stale legacy brain/solid SVGs in a resumed output directory from appearing in the new report.

The independent QC branch additionally publishes [genome-wide CpG bedMethyl](extended-qc.md) as an explicit final methylation deliverable; it is separate from classifier-specific locus pileups.

Parser inputs, augmented review JSON, per-model inference outputs, generated pileups and native section working files are not separately published by the new processes. Final report generation follows the existing report publication configuration. The native combined JSON remains unchanged; adapter annotations are added only to a separate private work-directory copy.

## Private ECR delivery

The process and clinical lock pin the **Linux/amd64 image manifest**, not the build attestation or OCI index:

```text
288525090463.dkr.ecr.us-west-2.amazonaws.com/oncoseq-leuk/classy-blood-nextflow@sha256:77aed6f284d94eda23f27a75080919a7b6888accb75cbbe73d7fd08c72121561
```

The [wrapper Dockerfile](../assets/containers/classy_blood_nextflow/Dockerfile) inherits the verified GHCR source manifest and clears its application entrypoint so AWS Batch can launch Nextflow's Bash task script. No classifier binary or model files are changed. The source manifest is `sha256:cd1f5ea7a1ecf140f7ff5da013018f5f5f22c257bf8ad861899d2882d6dc3e91` (Classy 0.22.0, modkit 0.6.3).

Workers use their AWS IAM ECR pull permissions; no GitHub PAT is required on workers for this route. The private repository must be reachable from the Batch compute environment and the upload must finish before tasks can pull it. An image tag is registered only after its layers and manifest have been uploaded. Build output alone does not prove the ECR tag exists.

Verify upload completion with the standard AWS credential chain:

```bash
aws ecr describe-images \
    --profile plm-sandbox --region us-west-2 \
    --repository-name oncoseq-leuk/classy-blood-nextflow \
    --image-ids imageTag=sha-8713032 \
    --query 'imageDetails[0].imageDigest' --output text
```

A build with attestations can report an OCI index digest here; the pipeline deliberately pins its amd64 child. After successful upload, resume the existing run with the same work directory and `-resume`. To omit ALMA3, add `--classy_blood_enable_alma3 false`. To revert report selection to the legacy classifier and omit the new execution, add `--classy_blood false`; its files still publish under `classy_legacy`.

The new BAM reader uses the existing isolated high-memory AWS queue, reserves the full high-worker CPU allocation and caps pileup processing at eight threads in the clinical profile. Its clinical baseline memory request is 64 GB, with existing retry/cap behavior.

## Validation boundary

Synthetic adapter tests verify native decision preservation, unavailable/disabled results, ALMA3 labeling, metadata/model provenance and stale-plot exclusion. Local Nextflow workflow tests execute the new processes with an explicitly synthetic CLI stand-in to verify arguments, channel contracts, ALMA3 defaults/opt-out, cfDNA configuration and publication filtering. Local report tests render final review HTML with synthetic blood results and verify the new model identities, attribution and model provenance.

The real wrapper was inspected and its executables exercised locally. Bounded synthetic pileup trials reached MARLIN and ALMA3 inference, but a 50% methylation fixture was rejected by Lamprey preprocessing; these trials are not validation of prediction accuracy or a complete four-model modified-BAM run. A subsequent adequate synthetic pileup trial completed all four models successfully in the real wrapper and confirmed their native combined JSON contract. This checks executable/schema compatibility, not predictive accuracy or the BAM preprocessing path. Production modified-BAM execution and AWS pulling remain deployment checks after the private upload completes. No GitHub or cloud pipeline smoke run was launched for this integration.
