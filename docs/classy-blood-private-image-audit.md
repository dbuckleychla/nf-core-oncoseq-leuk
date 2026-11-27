# Private Classy blood image: compatibility and AWS delivery

Audit date: 27 November 2026. This is an inspection and integration design; production classifier pins and AWS infrastructure have not been changed.

## Recommendation

Use **private ECR pull-through caching of GHCR**, following the two ECR resources in `../cpm_msq_pipeline_research/terraform/AWS/ecr.tf`. The existing sandbox AWS GHCR credential can read this image. ECR can therefore fetch the source image from GHCR while Batch workers authenticate to ECR using AWS IAM.

Two independent compatibility changes remain necessary:

1. The source image's `ENTRYPOINT ["classy"]` is incompatible with the normal Nextflow/AWS Batch Bash launcher. Prefer an upstream image with no application entrypoint. Otherwise use the tested, configuration-only [launcher wrapper](../assets/containers/classy_blood_nextflow/Dockerfile) in a private ECR repository.
2. Add a blood-specific workflow/module and output contract. This image cannot replace the current all-classifier container globally: required brain/solid outputs and the legacy ALMA model are absent, while ALMA3.1 and OPTiC need explicit report integration.

Direct worker-to-GHCR pulling is also supported by the existing Terraform registry-auth hook, but ECR caching avoids putting the GitHub credential on every Batch worker and avoids a GHCR dependency on every cold worker pull.

## Image inspected

| Item | Recorded value |
| --- | --- |
| Source tag | `ghcr.io/bwbioinfo/classy-blood:sha-8713032` |
| Source OCI index | `sha256:e593b00b44b6ad50b84045d3e2dac67f12780317c4a02f9f415983f4a20a28df` |
| Linux/amd64 manifest | `sha256:cd1f5ea7a1ecf140f7ff5da013018f5f5f22c257bf8ad861899d2882d6dc3e91` |
| Declared source revision | `8713032d4eba385ca041f1b8ec207282d6814ba9` |
| Runtime | Linux/amd64; executed locally under ARM-host Docker emulation |
| Classy executable/version | `/usr/local/bin/classy`; `classy 0.22.0` |
| Modkit version | `modkit 0.6.3` |
| Working directory | `/opt/classy` |
| Entrypoint / command | `["classy"]` / `["--help"]` |
| Local unpacked image size | 12,522,390,907 bytes, approximately 11.66 GiB |

The index also contains an unknown-platform attestation manifest. Use the **Linux/amd64 child digest** for platform-specific locks; the index digest is not the same identifier.

The image includes Bash. `python3` and `samtools` were not found on its default PATH. The existing combined classification command only needs Classy and modkit; Python-based QC must retain its separate runtime.

### Packaged classifier assets

The directory inspection found:

- `/opt/classy/models/MARLIN`, including its HDF5 model, feature RData and hg19/hg38/T2T probe files.
- `/opt/classy/models/Lamprey`, including probe and runtime assets.
- `/opt/classy/models/alma3/release-3.1.0`, including release, ensemble, taxonomy, CpG manifest and threshold metadata.
- `/opt/classy/models/optic`, including `optic.safetensors`, lineage and probe files.
- `/opt/classy/models/crossNN`, including runtime, static and probe assets.
- `/opt/classy/models/liftover_chains`.

There were no MPACT, Sturgeon or TUCAN model directories, and no legacy `/opt/classy/models/alma` directory. Asset presence establishes packaging, not successful inference for every model.

## Local execution evidence

The current invocation's main arguments are accepted by `classy combined`: input, output, sample, reference, genome, motif, `--use-pileup` and `--pileup-threads`. Supported genome values remain `hg19`, `hg38` and `t2t`.

A bounded local combined run used a synthetic precomputed bedMethyl file, retained MARLIN only and excluded the other routes. It completed successfully under 2 CPUs/4 GiB and emitted:

- Combined JSON with `combined.classification_scores`, `combined.classifier_runs`, `combined.charts` and `sample_label`.
- A `marlin` report carrying `model_label`, `model_version`, `model_version_source`, `inference` and `cutoff`.
- The familiar `*_combined_blood.svg` and `*_combined_task_pies.svg` files, plus a new blood hierarchy SVG and route-specific figures.

There were **zero matched probe observations**. This validates the launch and serialization path; it does not establish classification accuracy, sparse-sample adequacy, all-model inference, or equivalence to the old image. No clinical sample was used.

The supported exclusion names were confirmed by CLI validation: `marlin`, `crossnn`, `tucan`, `sturgeon`, `alma`, `alma3`, `lamprey`, `mpact`, `nanomix`, `optic`.

### Launcher incompatibility and tested fix

Running the source image with a normal Bash command produced `unrecognized subcommand '/bin/bash'`: Docker appended the command to its Classy entrypoint. The inspected Nextflow AWS Batch launcher submits a Bash command; its `AwsContainerOptionsMapper` does not map Docker `--entrypoint` into Batch container properties. A local Docker `--entrypoint` override therefore does not establish that the same configuration works on Batch.

The wrapper is:

```dockerfile
ARG SOURCE_IMAGE=ghcr.io/bwbioinfo/classy-blood@sha256:cd1f5ea7a1ecf140f7ff5da013018f5f5f22c257bf8ad861899d2882d6dc3e91
FROM ${SOURCE_IMAGE}
ENTRYPOINT []
CMD ["/bin/bash"]
```

It was built locally as `local/classy-blood-nextflow:8713032`. A Bash launch returned Classy 0.22.0 and modkit 0.6.3. Docker inspection confirmed that its filesystem layer list is identical to the source image; the change is launcher configuration only. It has not been pushed to any registry. Its eventual private ECR manifest digest must be recorded after publication, rather than treating a local image/config ID as a registry manifest digest.

## Fit with this workflow

### Why a global container swap will fail

`modules/local/clinical_threads/main.nf:CLASSY_COMBINED` requires blood, brain, solid, combined and nanomix SVGs. `modules/local/classifier_process/main.nf:PARSE_JSON_COMBINED` requires TSVs for legacy ALMA, both crossNN classifiers, MARLIN, Lamprey, MPACT, both Sturgeon classifiers, TUCAN and Nanomix.

The parser writes tables only for emitted `classifier_run_id` values. A valid blood-only JSON can therefore produce a successful parser exit followed by a Nextflow missing-output failure. Do not manufacture empty results for missing classifier families.

ALMA3.1 is not a rename of the legacy ALMA model. The inspected help explicitly describes its Rust/native implementation as an approximation whose decisions can differ from upstream, and retires ALMA3.0. Its default modification mode is `5mc`; `5mc-plus-5hmc` is an explicit alternative. Record this choice separately from the Dorado modified-base model.

The current review adapter enumerates legacy `alma` and always marks OPTiC disabled. It must gain explicit ALMA3/OPTiC support. OPTiC's native accepted-call decision must be retained; the existing R parser already documents that it is not simply `score >= cutoff`.

### Recommended implementation boundary

Introduce an opt-in local `CLASSY_BLOOD` module/subworkflow. Initially run it additively against the same classification-stage BAM as the existing workflow. Keep the old classifier path available for comparison and for brain/solid samples.

The initial blood command can retain MARLIN, Lamprey and OPTiC while explicitly disabling ALMA3 pending its validation:

```bash
classy combined \
  --input SYNTHETIC.bam \
  --output SYNTHETIC_classy_blood.json \
  --sample SYNTHETIC \
  --reference reference.fa \
  --genome hg38 \
  --use-pileup --motif CpG:CG --pileup-threads 8 \
  --exclude crossnn,tucan,sturgeon,alma,mpact,nanomix \
  --enable-alma3 false
```

This is a proposed blood-route invocation, not a validated all-blood result. Enable ALMA3 only after checking its native approximation, ONT modification mode and expected calls. Confirm OPTiC participation explicitly in the emitted combined metadata.

The module should have a blood-only output contract and separate filenames, with absent models represented as unavailable in reports. Reuse the existing SVG discovery/cropping/full-view presentation. Extend the provenance tab with Classy/modkit versions, model release metadata, source GHCR digest, executed wrapper/ECR digest and route/modification settings.

Preserve each model's native coverage and adequacy semantics. The zero-observation MARLIN probe still emitted probabilities: a successful process exit and a populated score vector must not by themselves mark a sample as an adequately supported classification. Record matched-probe counts and the native availability/accepted-call fields.

Route this new BAM reader to the isolated high-BAM queue and start with the existing Classy 8-thread/64-GiB budget; measure memory before reducing it. Add the new process to resource coverage checks and image locks. Clinical `withName` container pins override the process body's `params.classy_container` fallback, so `--classy_container` alone is not a complete clinical promotion mechanism.

### Model license requirements found

The bundled ALMA3 license authorizes internal hospital/laboratory use and evaluation. It requires preservation of its license/copyright in copies, attribution in shared products/services, and the statement **“Analysis performed using ALMA3.”** in reports shared with third parties. Commercial diagnostic/model services to third parties require a separate agreement. Its shared-product attribution is **“Powered by ALMA3.”**

The bundled OPTiC license is MIT. Preserve the bundled licenses in any wrapper or private cache. The configuration-only wrapper retains these files automatically.

## AWS credential verification

The MSQ source configuration identifies a GHCR Secrets Manager credential in the sandbox account. The tested us-west-2 credential is the existing `ecr-pullthroughcache/github-3` secret referenced by that Terraform configuration.

Using the sandbox AWS profile, the credential was retrieved into process memory, then tested against GHCR. No token, password, Authorization header or secret payload was printed or saved.

| Request | Result |
| --- | --- |
| Repository-scoped GHCR pull token for `bwbioinfo/classy-blood` | Success |
| Manifest for `sha-8713032` | HTTP 200; matched the local OCI index digest |
| Manifest for the Linux/amd64 child digest | HTTP 200 |
| Image-layer access via authenticated HEAD | HTTP 200 |

This establishes that the **current AWS-stored PAT has pull access to this specific image**. It does not establish that an ECR cache rule has been deployed or that workers can launch the image. The tested sandbox us-west-2 registry had no active GHCR `github` pull-through rule at audit time. The default AWS credential context also differed from the explicitly selected sandbox profile; select account/profile and region explicitly when deploying.

The existing completed-run metadata records Batch in the same sandbox account and us-west-2 region as the tested secret. Select `AWS_PROFILE=plm-sandbox` for the corresponding desktop deployment commands; the desktop's default credential context points to a different account.

## ECR pull-through cache design

Reuse the **pattern**, rather than applying or copying the whole MSQ stack. Its `aws_ecr_pull_through_cache_rule` and repository-creation template are appropriate. Its timestamp-triggered invalidation provisioners delete cache repositories/images on apply/destroy and should not be brought into this fork. Its broad worker `SecretsManagerReadWrite` attachment is also unnecessary for ECR-managed upstream authentication.

Use a fork-owned prefix such as `oncoseq-ghcr` in this repository's separate Terraform stack:

```hcl
resource "aws_ecr_pull_through_cache_rule" "classy_ghcr" {
  ecr_repository_prefix = "oncoseq-ghcr"
  upstream_registry_url = "ghcr.io"
  credential_arn        = var.ghcr_cache_secret_arn
}

resource "aws_ecr_repository_creation_template" "classy_ghcr" {
  prefix               = aws_ecr_pull_through_cache_rule.classy_ghcr.ecr_repository_prefix
  applied_for          = ["PULL_THROUGH_CACHE"]
  image_tag_mutability = "MUTABLE"
}
```

Supply the existing credential ARN, not a Terraform data source that reads its SecretString. ECR expects a secret in the cache's account and region, with the `ecr-pullthroughcache/` name prefix and the `username`/`accessToken` credential shape. Ensure the ECR pull-through service role can read it, and can decrypt a customer-managed KMS key if applicable.

The cached source reference would be:

```text
<account>.dkr.ecr.us-west-2.amazonaws.com/oncoseq-ghcr/bwbioinfo/classy-blood:sha-8713032
```

Warm the cache by pulling the source **tag** with a deployment identity, then inspect the cached image and pin the verified Linux/amd64 digest. Do not assume the first pull by an upstream digest will populate an empty cache, or that a cached index digest is its amd64 child. ECR refreshes tag content opportunistically, at most once per 24 hours; a workflow digest pin prevents silent promotion when a source tag moves.

The deployment/warm-up identity needs normal ECR pull access plus scoped `ecr:CreateRepository` and `ecr:BatchImportUpstreamImage` for the cache prefix. Prewarming allows worker roles to use ordinary ECR pull permissions for the pinned cached artifacts. ECR handles upstream GHCR credentials; workers do not need the GitHub PAT.

For the launcher fix, prefer a corrected upstream image and cache it directly. If that is unavailable, build the supplied wrapper from the verified cached source and publish it to a **private** ECR repository such as `oncoseq-leuk/classy-blood-nextflow`. Pin the resulting wrapper digest for the new blood module. The upstream image still comes from GHCR through ECR; only the launcher configuration is maintained locally.

## Alternative: direct GHCR pulls

This repository already implements worker-only registry authentication in `terraform/aws/templates/configure_registry_auth.py`. Setting `ghcr_secret_arn` grants narrowly scoped secret access to worker instance roles and installs an ECS startup hook that configures `ECS_ENGINE_AUTH_DATA` before the agent starts. The existing credential is compatible with the hook.

This route preserves direct source pulls but requires new/refreshed workers to pick up the launch-template user data and requires GHCR connectivity/authentication on cold pulls. Headnode `docker login` does not authenticate Batch workers. It also leaves the launcher and blood-only output-contract changes unresolved.

## Evidence limits and rollout boundary

The current audit verified architecture, CLI, packaged paths, a MARLIN-only serialization run, launcher-only wrapper behavior, and AWS PAT access. Full MARLIN/Lamprey/OPTiC/ALMA3 inference on a representative modified BAM, old/new model comparisons, report parsing and AWS Batch execution remain separate validation work.

No Terraform apply, ECR rule creation, cache import, remote image push, production pin change or clinical sample run was performed in this audit.
