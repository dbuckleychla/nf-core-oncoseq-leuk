# Institutional AWS Batch deployment

The pipeline baseline is parent repository `chusj-pigu/nf-core-oncoseq` release `v1.0.2`, commit `a663a10012e467ded972f1b98ded372d5f5ef2fa`. The adaptation PR targets this fork's `main` branch.

This fork has a separate Terraform stack under `terraform/aws`. It derives from FrankONTstein commit `46553118352caa218f82762f74c74ce635220177` and uses its institutional queue, storage, registry-authentication and IAM conventions. The stack uses its own Terraform state, resource names and S3 prefix. Its default name and prefix are `oncoseq-leuk`. The copied infrastructure and helper code retain the FrankONTstein MIT licence.

Terraform 1.16.4 and AWS provider 6.13.0 are pinned. Supply an existing S3 bucket and deployment region. The example uses the institutional research bucket, six explicit CPU instance sizes, and independent 512-vCPU CPU/basecalling/general GPU capacities. Basecalling uses 48-vCPU g7e.12xlarge/g6.12xlarge candidates; g7e requires Dorado/AMI acceptance testing. The general GPU fleet uses g6 sizes. The default optional VPC CIDR is `10.83.0.0/16`.

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

Clinical AWS execution requires the dedicated BAM high queue. Larger tasks reserve a whole 32-vCPU worker. Smaller BAM/read tasks reserve a whole 8-vCPU worker while their memory request is at most 48 GB. Retries exceeding that threshold move to the high queue. Local wrappers retain separate tool thread caps; imported nf-core modules retain their interfaces. `max_cpus` must be at least 32, reservations must match 32/8, and `max_memory` must not exceed 224 GB. Multi-shard Dorado jobs request 8 CPUs, 30 GB and one GPU, allowing two shards to share a two-GPU worker; single-shard runs retain the larger reservation. See [step budgets and failure handling](resource-allocation.md).

## Reproducibility

`clinical` pins Nextflow 26.04.6 and 36 Linux x86-64 image digests covering all 87 module processes. `nf-schema` remains pinned at 2.7.2. `images.lock.json` records the source tags, selected digests and process assignments. The clinical selectors take precedence over module image declarations. `scripts/check_release.py` checks coverage, config consistency and submodule revisions. Update the lock and its generated config together when changing an image.

The nf-core institutional config reference, PharmCAT 3.2.0 positions, pbsv tandem repeats and Severus PoN URLs use fixed commits. Clinical execution requires an explicit reference cache and a VEP cache matching the chosen assembly and VEP version. Basecalling requires complete versioned model names and pre-staged model directories. Clinical runs fail when these cached assets are absent.

The reference cache contains the selected FASTA and index, auxiliary files under their original basenames, and `dorado_models/` when basecalling. A direct FASTA path also permits auxiliary lookup in its parent directory. Existing S3 layouts can provide separate `dorado_model_cache` and `reference_aux_cache` directories. Use `--vep_cache` for the VEP cache root. Ordinary runs that download VEP publish it to that explicit root or `outdir/reference_cache`. The [institutional reference mapping](institutional-references.md) preserves compatible FrankONTstein paths and selects the public VEP 115 cache.

Version the actual reference bundle, model files, panels and run parameters with checksums and immutable S3 object versions as part of the validated release. The current candidate checks asset presence and version names. Content verification for those operator-supplied assets, host-package reproducibility, and clinical concordance remain open work. AWS bootstrap currently installs OS packages with `dnf`, even with a pinned AMI. The stack has been tested with a mocked provider and has not been deployed from this fork.

For SUP v5.2.0, use the compatibility evidence and GPU acceptance commands in the accompanying guide. [Dorado SUP v5.2.0](dorado-sup520.md)

## Verification

The candidate passes 20 mocked Terraform cases, 15 Python helper/registry/routing tests and 35 focused nf-test cases covering input/cache validation, GPU routing, BAM reservations, retry promotion/ceilings, huge-input scaling, aliased image pins and real wrapper argument capture with mock callers. These tests do not run inference. Registry tests cover anonymous and authenticated worker startup. Nextflow lint reports no errors in the checked scripts. The pipeline compiles through its normal `--version` entry point.

```bash
terraform -chdir=terraform/aws init -backend=false
terraform -chdir=terraform/aws validate
terraform -chdir=terraform/aws test
python3 scripts/check_release.py
python3 -m unittest discover -s tests -p 'test_*.py'
nf-test test tests/config --profile aws,clinical
nf-test test tests/modules/local/clinical_threads --profile aws,clinical
nf-test test tests/subworkflows/local/utils_nfcore_oncoseq_pipeline/validation.nf.test --profile ''
nf-test test tests/subworkflows/local/utils_nfcore_oncoseq_pipeline/references.nf.test --profile ''
```

The full `test,docker` run requires the missing `test_data/basecalling/samplesheets/samplesheet.csv` and its data. Broad nf-core lint also reports template, schema and vendored-module differences, including expectations that conflict with this fork's pinned institutional configuration. These results are tracked with fixture and CI restoration in `oncoseq-y3t`. GPU execution, reference-content verification and clinical concordance remain pending.

### AWS CLI staging fails with missing zlib

If `.command.err` reports that the mounted AWS CLI cannot load `libz.so.1`, the failure occurs during S3 staging, before the analysis command runs. Nextflow mounts the worker's CLI installation into the tool container, but the container may not include that library. Downloaded CLI installations now include a copy of the worker's zlib inside the mounted tree and a wrapper that sets the CLI library search path. Existing-AMI CLI installations remain operator-managed and must provide their own container-compatible dependencies.

After updating this checkout, run `terraform plan` and `terraform apply` in `terraform/aws` to publish the corrected worker bootstrap. Existing EC2 workers do not rerun user data; drain and replace affected workers before resuming the workflow. Re-export the Terraform environment if needed, then rerun the original Nextflow command with `-resume`.

For an existing worker, an administrator can apply the equivalent repair below **on the EC2 worker**, outside task containers. This applies to this stack's downloaded CLI path:

```bash
sudo bash <<'REPAIR'
set -euo pipefail
dnf install -y zlib
cli_root=/opt/aws-cli/v2/current
test -x "$cli_root/dist/aws"
install -d "$cli_root/runtime-libs"
cp -L /lib64/libz.so.1 "$cli_root/runtime-libs/libz.so.1"
rm -f "$cli_root/bin/aws"
cat > "$cli_root/bin/aws" <<'WRAPPER'
#!/bin/sh
set -eu
cli_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
export LD_LIBRARY_PATH="$cli_root/runtime-libs${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
exec "$cli_root/dist/aws" "$@"
WRAPPER
chmod 0755 "$cli_root/bin/aws"
"$cli_root/bin/aws" --version
REPAIR
```

Verify from the exact failing image, with the CLI mounted as Nextflow mounts it:

```bash
docker run --rm \
  -v /opt/aws-cli:/opt/aws-cli:ro \
  --entrypoint /opt/aws-cli/v2/current/bin/aws \
  quay.io/biocontainers/cramino@sha256:a0ed72cbb52240e6cb0085bea16f5356b3d790621281dd11cecc2f5173e5f791 \
  --version
```

A host-only `aws --version` check cannot establish container compatibility. Run the container check on each worker type or representative image before clinical acceptance.
