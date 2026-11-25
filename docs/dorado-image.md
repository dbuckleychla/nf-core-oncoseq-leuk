# Build and publish Dorado 2.2.0

The tracked [Dockerfile](../containers/dorado-2.2.0/Dockerfile) installs ONT's
Linux x64 Dorado 2.2.0 release in a digest-pinned CUDA 12.8.1/cuDNN Ubuntu 24.04
runtime. The image name is `ghcr.io/dbuckleychla/dorado-2.2.0`, with tag `latest`.
The pipeline will use its Linux amd64 SHA-256 digest after acceptance, rather
than the mutable tag. The current pipeline image remains pinned to Dorado 1.3.0
until that promotion is validated.

BuildKit checks the official archive's SHA-256 before extraction:
`c2376701a48bc6227a9cb5ad1d0b9ada60b5f3c8760a657066131fdf0d9545c4`.
The 3,473,625,059-byte archive was streamed from ONT's HTTPS release endpoint
and hashed on 25 November 2026. No archive or model files are committed here.
The complete extracted release retains its bundled libraries and licensing files.
A separate download stage keeps the archive out of the final image layers; the
runtime contains the pinned CUDA base and extracted Dorado distribution. The
Dockerfile frontend is digest-pinned as well.
[Official release](https://github.com/nanoporetech/dorado/releases/tag/v2.2.0)
[Versioned model catalogue](https://github.com/nanoporetech/dorado/blob/v2.2.0/dorado/models/models.cpp)

## Local build

Use a current Docker Desktop/Buildx installation. On Apple Silicon, build the
Linux amd64 image explicitly; this image is for the x86-64 Batch workers.
The Dockerfile rejects other target architectures. Use the external volume for
large artifacts when mounted. Docker Desktop's disk-image location controls
image and build-cache storage; setting shell `TMPDIR` alone does not relocate
Docker layers. Allow space for the 3.5 GB archive, its extracted libraries, the
CUDA base image and build cache.

Run from the repository root:

```bash
DORADO_IMAGE='ghcr.io/dbuckleychla/dorado-2.2.0:latest'
docker buildx build --platform linux/amd64 --load \
    --build-arg VCS_REF="$(git rev-parse HEAD)" \
    --tag "$DORADO_IMAGE" containers/dorado-2.2.0
```

This narrow build context contains only the Dockerfile and `.dockerignore`.
Software is downloaded and verified during the build; models and clinical data
are not packaged. The NVIDIA entrypoint remains shell-compatible for Nextflow.
There is no automated image-publishing workflow.

## Verify and push as the owner

```bash
docker image inspect "$DORADO_IMAGE" \
    --format '{{.Os}}/{{.Architecture}} {{index .Config.Labels "org.opencontainers.image.version"}}'
docker run --rm --platform linux/amd64 "$DORADO_IMAGE" dorado --version
docker run --rm --platform linux/amd64 "$DORADO_IMAGE" dorado basecaller --help
docker run --rm --platform linux/amd64 "$DORADO_IMAGE" dorado demux --help
docker run --rm --platform linux/amd64 "$DORADO_IMAGE" dorado trim --help
docker run --rm --platform linux/amd64 "$DORADO_IMAGE" dorado download --list-structured
```

Require Linux/amd64 and executable Dorado 2.2.0. Confirm the existing
`dna_r10.4.1_e8.2_400bps_sup@v5.2.0` and matching
`dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mCG_5hmCG@v2` models. This is a software
upgrade; it does not select the new v6.x basecalling models or change downstream
callers. If Apple Silicon emulation cannot execute the binary, retain that
failure and run these checks on the x86-64 GPU worker before promotion.

Log in interactively; enter your PAT at the password prompt. GHCR uses a classic
PAT with package read/write scopes. Do not put its value in commands, repository
files, build arguments or GitHub Actions secrets for this image.

```bash
docker login ghcr.io --username dbuckleychla
docker push "$DORADO_IMAGE"
docker buildx imagetools inspect "$DORADO_IMAGE"
```

In GitHub, open your **Packages → dorado-2.2.0 → Package settings**:

- Set visibility to **Public** for anonymous pulls.
- Disable **Inherit access from source repository** if enabled.
- In **Manage access**, retain only `dbuckleychla` with Write/Admin access.
- In **Manage Actions access**, grant no repository workflow Write/Admin access.

The OCI source label connects the package to this repository for provenance;
separate package permissions prevent future repository collaborators from
inheriting publishing access. Anyone can pull a public package; publishing still
requires an authorized account token. Verify these settings before promotion.
[GHCR access permissions](https://docs.github.com/en/packages/learn-github-packages/about-permissions-for-github-packages)

Record the published **Linux amd64 manifest digest**, not the local image ID.
For a manifest list, select its Linux/amd64 entry. For a single Linux/amd64
manifest, use that manifest's registry digest. The runtime reference will be:

```text
ghcr.io/dbuckleychla/dorado-2.2.0@sha256:<published-linux-amd64-digest>
```

## GPU acceptance before promotion

The [acceptance workflow](../tests/acceptance/dorado/main.nf) exercises the actual
sharded Dorado wrapper with one public-test POD5 and at most 32 reads. It verifies
Dorado 2.2.0, nonempty valid BAM output, the requested model identity and MM/ML
modification tags. Use public ONT test input with matching R10.4.1 E8.2, 400 bps,
5 kHz chemistry. Keep test input, complete model directories and generated
artifacts in external/runtime storage, outside the repository.

On a Linux NVIDIA machine with Docker and NVIDIA Container Toolkit:

```bash
nextflow run tests/acceptance/dorado/main.nf -profile docker \
    --image "$DORADO_DIGEST" \
    --pod5 "$TEST_POD5" --model "$SUP_MODEL" --mod_model "$CPG_MODEL" \
    --outdir "$SMOKE_OUT" -work-dir "$SMOKE_WORK"
```

Set `DORADO_DIGEST` to the published digest reference above. Set `SUP_MODEL` and
`CPG_MODEL` to the complete versioned model directories named in the verification
section. `TEST_POD5`, `SMOKE_OUT` and `SMOKE_WORK` are operator-supplied runtime
paths. The workflow uses the fixed synthetic label `SYNTHETIC_DORADO_SMOKE`.

To test the existing AWS Batch deployment, use its queue names and registry/CLI
configuration without changing infrastructure:

```bash
AWS_PROFILE=plm-sandbox nextflow run tests/acceptance/dorado/main.nf \
    -profile aws -c aws.generated.config \
    --aws_region us-west-2 \
    --aws_basecall_queue "$BASECALL_QUEUE" --aws_queue "$CPU_QUEUE" \
    --image "$DORADO_DIGEST" \
    --pod5 "$TEST_POD5_URI" --model "$SUP_MODEL_URI" --mod_model "$CPG_MODEL_URI" \
    --outdir "$SMOKE_OUT_URI" -work-dir "$SMOKE_WORK_URI"
```

Obtain queue names and the generated configuration from the current Terraform
deployment. Use S3 paths for Batch inputs, output and work directories. The
basecaller requests one GPU, 8 CPUs and 30 GB, with no retries; validation runs
on the CPU queue. Public image pulls require no registry-write credential.

Retain `validated.tsv`, the BAM, `versions.yml`, the basecall task log containing
`nvidia-smi`, input/model checksums, published image digest and actual Batch
instance type/AMI identity. Check the GPU's driver against the container's CUDA
requirements. ONT additionally documents a minimum driver of 570.26 for its
supported Linux Blackwell cards; the actual runtime constraints also apply.
A failed check must not trigger model changes, automatic infrastructure changes
or promotion. [Dorado 2.2.0 platforms](https://github.com/nanoporetech/dorado/blob/v2.2.0/README.md#platforms)

After acceptance passes, update the Dorado entry in `images.lock.json`, regenerate
`conf/containers_clinical_amd64.config`, and change local Dorado defaults to the
same digest. This covers basecalling (including shards), demultiplexing, trimming
and model-download commands. Run affected configuration/sharding tests and
`python3 scripts/check_release.py`. Retain the existing 1.3.0 digest for rollback
until the candidate is accepted. No default image pin changes merely because a
new `latest` tag exists.
