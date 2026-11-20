# Institutional AWS Batch deployment

The pipeline baseline is parent repository `chusj-pigu/nf-core-oncoseq` release `v1.0.2`, commit `a663a10012e467ded972f1b98ded372d5f5ef2fa`. The adaptation PR targets this fork's `main` branch.

This fork has a separate Terraform stack under `terraform/aws`. It derives from FrankONTstein commit `46553118352caa218f82762f74c74ce635220177` and uses its institutional queue, storage, registry-authentication and IAM conventions. The stack uses its own Terraform state, resource names and S3 prefix. Its default name and prefix are `oncoseq-leuk`. The copied infrastructure and helper code retain the FrankONTstein MIT licence.

Terraform 1.13.2 and AWS provider 6.13.0 are pinned. Supply an existing S3 bucket and deployment region. The example uses the institutional research bucket, six explicit CPU instance sizes, and independent 512-vCPU CPU/basecalling/general GPU capacities. Both GPU fleets use g6e sizes. Optional g7e basecalling candidates require Dorado/AMI acceptance testing; the general fleet retains g6e for caller-image compatibility. The default optional VPC CIDR is `10.83.0.0/16`.

The example pins AWS ECS-optimised AL2023 x86-64 release `2023.0.20261113` in `us-west-2`: CPU `ami-0f9003a6aee571ca8` and NVIDIA GPU `ami-0c050c3f1924995a9`. These IDs were resolved on 20 November 2026 from the public SSM parameters `/aws/service/ecs/optimized-ami/amazon-linux-2023/recommended` (version 165) and `/aws/service/ecs/optimized-ami/amazon-linux-2023/gpu/recommended` (version 94). EC2 confirmed both are public, available, x86-64 images owned by AWS account `591542846629`. Both contain ECS agent 1.108.3 and Docker 25.0.18. IDs are region-specific; select and pin equivalent images when changing region. Worker boot, CUDA inference and Dorado compatibility still require deployment acceptance tests.

```bash
cp terraform/aws/terraform.tfvars.example terraform/aws/terraform.tfvars
# Fill deployment values and choose the intended state backend/workspace.
terraform -chdir=terraform/aws init
terraform -chdir=terraform/aws plan -out=deployment.tfplan
# Review the plan before applying it to the selected AWS account.
terraform -chdir=terraform/aws apply deployment.tfplan
```

Use a distinct backend state key or an independently managed local state for this stack. An existing VPC can be selected with `vpc_id` and `subnet_ids`. Work and results IAM writes are scoped to the configured prefix. The coordinator requires the exported `coordinator_policy_arn`. Supply readable input/reference bucket ARNs and relevant KMS keys in tfvars.

All 15 GHCR image digests in the current lock were verified anonymously readable on 20 November 2026. Public pulls use the default `ghcr_secret_arn = null`; workers receive no registry-secret policy or authentication hook. If private images are added, supply an existing same-region secret ARN with `username`/`accessToken` or `username`/`password` JSON. Workers retrieve it before ECS starts. Supply `ghcr_secret_kms_key_arn` only when that secret uses a customer-managed key.

## Launching

Install Nextflow 26.04.6 and initialise the recorded module commits. The environment helper reads this stack's outputs and emits shell-quoted `ONCOSEQ_*` variables. Use a distinct run ID for each analysis.

```bash
git submodule update --init --recursive
python3 scripts/check_release.py
python3 scripts/aws_batch_env.py --terraform-dir terraform/aws > /tmp/oncoseq-aws.env
source /tmp/oncoseq-aws.env
NXF_VER=26.04.6 nextflow run . -profile aws,clinical \
  --run_id validation-001 \
  -params-file validation.params.json
```

The helper also accepts a saved full `terraform output -json` export with `--terraform-outputs`. `ONCOSEQ_WORK_PREFIX` and `ONCOSEQ_RESULTS_PREFIX` combine with `--run_id`. Explicit CLI or params-file settings retain Nextflow's normal precedence. Runs require S3 work and result locations, a region, task role and CPU queue.

`aws` enables Docker, pins `nf-amazon@3.9.2`, and routes Dorado to `aws_basecall_queue` with one accelerator. Other GPU-labelled processes use `aws_gpu_queue`. The two queues have separate compute environments and independently configurable instance families and maximum vCPUs. Clair3, ClairS-TO and DeepSomatic use CPU execution in this candidate. General resource ceilings are 32 CPUs and 224 GB. The basecalling fleet uses `basecall_instance_types` and `basecall_max_vcpus`. The general fleet uses `gpu_instance_types` and `gpu_max_vcpus`. Choose GPU instance sizes that can satisfy Dorado's requested CPU and host memory resources.

When dedicated BAM queues are enabled, larger tasks reserve a 32-vCPU worker. Selected smaller BAM tasks reserve an 8-vCPU worker while their memory request is at most 48 GB. Retries exceeding that threshold move to the high queue. These reservations also set the thread counts used by the existing wrappers. Confirm resource use with representative clinical datasets.

## Reproducibility

`clinical` pins Nextflow 26.04.6 and 36 Linux x86-64 image digests covering all 87 module processes. `nf-schema` remains pinned at 2.7.2. `images.lock.json` records the source tags, selected digests and process assignments. The clinical selectors take precedence over module image declarations. `scripts/check_release.py` checks coverage, config consistency and submodule revisions. Update the lock and its generated config together when changing an image.

The nf-core institutional config reference, PharmCAT 3.2.0 positions, pbsv tandem repeats and Severus PoN URLs use fixed commits. Clinical execution requires an explicit reference cache and a VEP cache matching the chosen assembly and VEP version. Basecalling requires complete versioned model names and pre-staged model directories. Clinical runs fail when these cached assets are absent.

The reference cache contains the selected FASTA and index, auxiliary files under their original basenames, and `dorado_models/` when basecalling. A direct FASTA path also permits auxiliary lookup in its parent directory. Use `--vep_cache` for the VEP cache root. Ordinary runs that download VEP publish it to that explicit root or `outdir/reference_cache`.

Version the actual reference bundle, model files, panels and run parameters with checksums and immutable S3 object versions as part of the validated release. The current candidate checks asset presence and version names. Content verification for those operator-supplied assets, host-package reproducibility, and clinical concordance remain open work. AWS bootstrap currently installs OS packages with `dnf`, even with a pinned AMI. The stack has been tested with a mocked provider and has not been deployed from this fork.

For SUP v5.2.0, use the compatibility evidence and GPU acceptance commands in the accompanying guide. [Dorado SUP v5.2.0](dorado-sup520.md)

## Verification

The candidate passes 20 mocked Terraform cases, 14 Python helper/registry tests and 10 focused nf-test cases covering input validation, GPU routing, BAM reservations and aliased image pins. Registry tests cover anonymous and authenticated worker startup. Nextflow lint reports no errors in the checked scripts. The pipeline compiles through its normal `--version` entry point.

```bash
terraform -chdir=terraform/aws init -backend=false
terraform -chdir=terraform/aws validate
terraform -chdir=terraform/aws test
python3 scripts/check_release.py
python3 -m unittest discover -s tests -p 'test_aws_*.py'
nf-test test tests/config --profile aws,clinical
nf-test test tests/subworkflows/local/utils_nfcore_oncoseq_pipeline/validation.nf.test --profile ''
```

The full `test,docker` run requires the missing `test_data/basecalling/samplesheets/samplesheet.csv` and its data. Broad nf-core lint also reports template, schema and vendored-module differences, including expectations that conflict with this fork's pinned institutional configuration. These results are tracked with fixture and CI restoration in `oncoseq-y3t`. GPU execution, reference-content verification and clinical concordance remain pending.
