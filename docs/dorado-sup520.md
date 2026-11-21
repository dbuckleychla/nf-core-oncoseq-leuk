# Dorado SUP v5.2.0 compatibility

The clinical candidate pins `ghcr.io/dbuckleychla/frankontstein@sha256:15c315c121d37b5c115d3d965836e60d56756498ca2deebe066a373d41e96785` for Linux x86-64. Its OCI configuration and build history identify Dorado **1.3.0**, installed from ONT's versioned Linux archive. The image uses the CUDA 12.3.2 runtime. The replacement image configuration was inspected on 21 November 2026 with
`skopeo inspect --config`. It adds preprocessing dependencies to the same
Dorado base release; it is not a Dorado version upgrade. The digest was pulled,
but a local runtime check failed after Docker exhausted scratch storage. Run
the GPU smoke test below on the intended worker before production acceptance.
The parent image-update request is tracked as `oncoseq-4f7`.

Dorado 1.3.0 explicitly lists `dna_r10.4.1_e8.2_400bps_sup@v5.2.0` for R10.4.1 E8.2, 400 bps, 5 kHz DNA data. It also lists matching `5mCG_5hmCG@v1` and `@v2` CpG models, and `5mC_5hmC@v1` and `@v2` all-context models. The software release supports the requested basecalling model. [Dorado model catalogue at the release commit](https://github.com/nanoporetech/dorado/blob/6ea40018974f30d436122f17e422bdc55d8d7945/dorado/models/models.cpp#L308)

GPU execution remains unverified in this fork. Use the selected Batch NVIDIA AMI and instance type for the acceptance test. The local development host has no NVIDIA GPU, and the original base image alone contains approximately 5.4 GiB of compressed layers. Dorado's documented Linux driver minimum is 525.105, and the image's NVIDIA runtime constraints also apply. [Dorado platform requirements](https://github.com/nanoporetech/dorado/blob/6ea40018974f30d436122f17e422bdc55d8d7945/README.md#platforms)

## GPU smoke test

On a Linux NVIDIA worker with Docker and NVIDIA Container Toolkit, place a representative POD5 at `smoke/input.pod5` and the complete versioned model directories under `models/`. Use data from the chemistry and sampling rate above. The commands call at most 32 reads and retain version, driver, header, and log evidence. The CpG modification model below is an explicit test choice.

```bash
mkdir -p smoke
DORADO_IMAGE='ghcr.io/dbuckleychla/frankontstein@sha256:15c315c121d37b5c115d3d965836e60d56756498ca2deebe066a373d41e96785'
docker run --rm --gpus all --entrypoint bash \
  -v "$PWD/models:/models:ro" -v "$PWD/smoke:/smoke" -w /smoke \
  "$DORADO_IMAGE" -c '
    set -euo pipefail
    dorado --version > dorado-version.txt 2>&1
    nvidia-smi > nvidia-smi.txt
    dorado download --list-structured > available-models.json
    grep -Fq "dna_r10.4.1_e8.2_400bps_sup@v5.2.0" available-models.json
    dorado basecaller \
      /models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0 \
      /smoke/input.pod5 \
      --device cuda:all --max-reads 32 \
      --modified-bases-models /models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mCG_5hmCG@v2 \
      > smoke.bam 2> basecall.log
  '

SAMTOOLS_IMAGE='community.wave.seqera.io/library/htslib_samtools@sha256:a55ddea590e567a91df592300a960aa534cfc1bd16e7623e3938ec21f4f3df15'
docker run --rm --entrypoint bash \
  -v "$PWD/smoke:/smoke" -w /smoke "$SAMTOOLS_IMAGE" -c '
    set -euo pipefail
    samtools quickcheck -u smoke.bam
    samtools view -H smoke.bam > header.sam
    grep -Fq "dna_r10.4.1_e8.2_400bps_sup@v5.2.0" header.sam
    samtools view -c smoke.bam > read-count.txt
    samtools view -c -d MM smoke.bam > modified-read-count.txt
    samtools view -c -d ML smoke.bam > probability-read-count.txt
    test "$(cat read-count.txt)" -gt 0
    test "$(cat modified-read-count.txt)" -gt 0
    test "$(cat probability-read-count.txt)" -gt 0
  '
```

A passing test demonstrates model loading, CUDA inference, valid nonempty BAM output, the requested model identity, and modification-tag emission. Retain the image digest, model checksums, AMI ID, GPU type, input checksum and generated evidence together. Repeat through `DORADO_BASECALL` on Batch with the same inputs to exercise Nextflow staging and queue routing. A representative full dataset then measures throughput, memory and result concordance.

## Downstream model pairing

Oncoseq currently chooses `r1041_e82_400bps_sup_v500` for every SUP input. It also omits Dorado `--emit-moves`. FrankONTstein accepts external `r1041_e82_400bps_sup_v520_with_mv` weights, enables move emission and Clair3 dwell-time features, and checks the model pairing. Establishing equivalent SUP v5.2.0 processing through Clair3 requires that additional work and validation. It is tracked as `oncoseq-abz`. [Current caller selection](../subworkflows/local/variant_calling/snp_calling.nf) [Current Dorado wrapper](../modules/local/dorado/main.nf) [FrankONTstein basecalling contract](https://github.com/dbuckleychla/FrankONTstein/blob/46553118352caa218f82762f74c74ce635220177/documentation/basecalling.md)
