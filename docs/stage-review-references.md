# Stage institutional references for Oncoseq

The copy block requires only AWS CLI and performs direct S3-to-S3 copies.
It downloads no reference files locally. Use your usual AWS profile; these
private bucket operations use authenticated access. The later COSMIC and
manifest preparation steps use existing local files, Python 3 and Docker.

The source is the existing institutional reference prefix. The destination is
`s3://plm-cpm-sandbox-research-usw2/ONT/Oncoseq_refdata/`. No AWS uploads were
performed during development.

## Copy existing resources

```bash
set -euo pipefail
export AWS_DEFAULT_REGION=us-west-2
export SOURCE_REF_ROOT='s3://plm-cpm-sandbox-research-usw2/ONT/FrankONTstein_refdata'
export ONCOSEQ_REF_ROOT='s3://plm-cpm-sandbox-research-usw2/ONT/Oncoseq_refdata'
aws s3 cp "$SOURCE_REF_ROOT/hg38/fasta/hg38.analysisSet.fa" \
    "$ONCOSEQ_REF_ROOT/hg38/fasta/hg38.analysisSet.fa" --only-show-errors
aws s3 cp "$SOURCE_REF_ROOT/hg38/fasta/hg38.analysisSet.fa.fai" \
    "$ONCOSEQ_REF_ROOT/hg38/fasta/hg38.analysisSet.fa.fai" --only-show-errors
aws s3 cp "$SOURCE_REF_ROOT/hg38/targets/DARGAN_OCCRA_ALMA_CHLA_wishlist.GENCODE_v50.hg38.bed" \
    "$ONCOSEQ_REF_ROOT/hg38/targets/DARGAN_OCCRA_ALMA_CHLA_wishlist.GENCODE_v50.hg38.bed" --only-show-errors
aws s3 cp "$SOURCE_REF_ROOT/hg38/targets/DARGAN_OCCRA_ALMA_CHLA_wishlist.GENCODE_v50.10000bp-buffer.hg38.bed" \
    "$ONCOSEQ_REF_ROOT/hg38/targets/DARGAN_OCCRA_ALMA_CHLA_wishlist.GENCODE_v50.10000bp-buffer.hg38.bed" --only-show-errors
aws s3 cp "$SOURCE_REF_ROOT/hg38/nasvar/GCF_000001405.40_GRCh38.p14_genomic.gff" \
    "$ONCOSEQ_REF_ROOT/hg38/nasvar/GCF_000001405.40_GRCh38.p14_genomic.gff" --only-show-errors
aws s3 cp "$SOURCE_REF_ROOT/hg38/gnomad/genomes/gnomad_genome_4.0_lite.vcf.gz" \
    "$ONCOSEQ_REF_ROOT/hg38/gnomad/genomes/gnomad_genome_4.0_lite.vcf.gz" --only-show-errors
aws s3 cp "$SOURCE_REF_ROOT/hg38/gnomad/genomes/gnomad_genome_4.0_lite.vcf.gz.tbi" \
    "$ONCOSEQ_REF_ROOT/hg38/gnomad/genomes/gnomad_genome_4.0_lite.vcf.gz.tbi" --only-show-errors
aws s3 cp "$SOURCE_REF_ROOT/hg38/gnomad/exomes/gnomad_exome_4.0_lite.vcf.gz" \
    "$ONCOSEQ_REF_ROOT/hg38/gnomad/exomes/gnomad_exome_4.0_lite.vcf.gz" --only-show-errors
aws s3 cp "$SOURCE_REF_ROOT/hg38/gnomad/exomes/gnomad_exome_4.0_lite.vcf.gz.tbi" \
    "$ONCOSEQ_REF_ROOT/hg38/gnomad/exomes/gnomad_exome_4.0_lite.vcf.gz.tbi" --only-show-errors
aws s3 cp "$SOURCE_REF_ROOT/hg38/clinvar/clinvar.hg38.vcf.gz" \
    "$ONCOSEQ_REF_ROOT/hg38/clinvar/clinvar.hg38.vcf.gz" --only-show-errors
aws s3 cp "$SOURCE_REF_ROOT/hg38/clinvar/clinvar.hg38.vcf.gz.tbi" \
    "$ONCOSEQ_REF_ROOT/hg38/clinvar/clinvar.hg38.vcf.gz.tbi" --only-show-errors
aws s3 cp "$SOURCE_REF_ROOT/hg38/hgnc/clinvar_pathogenic_common.vcf.gz" \
    "$ONCOSEQ_REF_ROOT/hg38/hgnc/clinvar_pathogenic_common.vcf.gz" --only-show-errors
aws s3 cp "$SOURCE_REF_ROOT/hg38/hgnc/clinvar_pathogenic_common.vcf.gz.tbi" \
    "$ONCOSEQ_REF_ROOT/hg38/hgnc/clinvar_pathogenic_common.vcf.gz.tbi" --only-show-errors

aws s3 cp "$SOURCE_REF_ROOT/dorado-models/" \
    "$ONCOSEQ_REF_ROOT/dorado-models/" --recursive --only-show-errors
# Copy only VEP 115.
aws s3 cp "$SOURCE_REF_ROOT/hg38/vep-cache/homo_sapiens/115_GRCh38/" \
    "$ONCOSEQ_REF_ROOT/hg38/vep-cache/homo_sapiens/115_GRCh38/" --recursive --only-show-errors
```

## Prepare auxiliary assets if not already staged

This separate step downloads the small Oncoseq-specific auxiliary assets from
the pinned upstream URLs. Skip it if your destination auxiliary cache is ready.

```bash
export LOCAL_REFS="$PWD/oncoseq-ref-staging"
# Oncoseq-specific auxiliary assets, from the pinned URLs in this repository.
mkdir -p "$LOCAL_REFS/hg38/aux"
while IFS= read -r url; do
    curl --fail --location --retry 3 "$url" \
        -o "$LOCAL_REFS/hg38/aux/${url##*/}"
done < assets/reference_aux_hg38.urls
(cd "$LOCAL_REFS/hg38/aux" && sha256sum *.bgz *.csi *.bed *.gz > SHA256SUMS)
aws s3 sync "$LOCAL_REFS/hg38/aux/" "$ONCOSEQ_REF_ROOT/hg38/aux/" --only-show-errors
```

The explicit VCF names above match the institutional reference mapping. Missing
objects stop the script. The ClinVar input is the chromosome-compatible
`clinvar.hg38.vcf.gz` derivative, rather than the original numeric-contig file.
The review uses the unpadded target BED; the analytical preset uses the buffered
BED. Source catalogue subsets remain subsets after copying.

## Prepare and upload COSMIC

Supply the local path to your licensed targeted-screen VCF. Set `LOCAL_REFS`
to your existing local reference tree, with the analytical FASTA at
`hg38/fasta/hg38.analysisSet.fa`. The S3 copy block does not populate this tree.
Preparation uses that same full FASTA and creates a separate indexed derivative.

```bash
export REVIEW_IMAGE='community.wave.seqera.io/library/oncoseq_review_1_0_0@sha256:b143e350e1dae33ad26758fc6f9ccf2392cb8b49881de73960f42525602d5eb5'
export COSMIC_SOURCE='/Volumes/DNB_External_1/temp_SSD_data/COSMIC/Cosmic_CompleteTargetedScreensMutant_Normal_v105_GRCh38.vcf.gz'
export COSMIC_DIR="$(dirname "$COSMIC_SOURCE")"
mkdir -p "$LOCAL_REFS/hg38/cosmic"

docker run --rm --platform linux/amd64 \
    -v "$PWD:$PWD" -v "$LOCAL_REFS:/ref" -v "$COSMIC_DIR:/cosmic-source:ro" \
    -w "$PWD" --entrypoint python3 "$REVIEW_IMAGE" \
    bin/review_report.py prepare-vcf \
    --source "/cosmic-source/$(basename "$COSMIC_SOURCE")" \
    --fasta /ref/hg38/fasta/hg38.analysisSet.fa \
    --output /ref/hg38/cosmic/cosmic.targeted.v105.GRCh38.vcf.gz \
    --skip-masked-reference

# Includes the VCF, index, and preparation audit. Keep the source archive locally.
aws s3 sync "$LOCAL_REFS/hg38/cosmic/" "$ONCOSEQ_REF_ROOT/hg38/cosmic/" --only-show-errors
```

The analysis-set FASTA masks some reference regions with `N` bases (including
parts of chrY). `--skip-masked-reference` archives those original VCF records
in `.masked-reference.vcf.gz` and excludes them from the derivative, recording
counts by contig and checksums in the preparation audit. It still fails on REF
mismatches at unmasked bases. Default preparation remains strict without this
flag. A source index is unnecessary because preparation reads sequentially.

The matching genome-screen product is optional. Prepare it separately and add
another `cosmic_*` manifest entry if you obtain it; never substitute another
assembly or merge the products' recurrence counts.

## Optional prelaunch checksum manifest

Skip this section for development runs: the clinical profile uses the bundled
`assets/institutional_hg38.review.json` mapping and needs no uploaded manifest.
Reference pinning before launch is tracked as `oncoseq-u7o`.

This optional step hashes existing local reference files beneath `LOCAL_REFS`, using
the same relative paths as the S3 objects. It does not download them. If you
already have a checksum-pinned manifest for these exact bytes, update its URIs
to the destination prefix instead of regenerating hashes.

Set the two release labels from your downloaded source metadata before running
this block. Checksums pin the actual copied content; labels must identify the
ClinVar release and the exact rescue-list preparation, not guessed versions.

```bash
export CLINVAR_RELEASE='REPLACE_WITH_SOURCE_RELEASE'
export RESCUE_RELEASE='REPLACE_WITH_RESCUE_LIST_RELEASE'

python3 - <<'PY'
import json, os
from pathlib import Path
root = Path(os.environ['LOCAL_REFS'])
remote = os.environ['ONCOSEQ_REF_ROOT'].rstrip('/')
clinvar = os.environ['CLINVAR_RELEASE']
rescue = os.environ['RESCUE_RELEASE']
if 'REPLACE' in clinvar or 'REPLACE' in rescue:
    raise SystemExit('Set actual CLINVAR_RELEASE and RESCUE_RELEASE first')
files = {
    'fasta': ('hg38/fasta/hg38.analysisSet.fa', 'GRCh38-analysisSet'),
    'fai': ('hg38/fasta/hg38.analysisSet.fa.fai', 'GRCh38-analysisSet'),
    'targets': ('hg38/targets/DARGAN_OCCRA_ALMA_CHLA_wishlist.GENCODE_v50.hg38.bed', 'GENCODE-v50-assay-targets'),
    'gff': ('hg38/nasvar/GCF_000001405.40_GRCh38.p14_genomic.gff', 'GRCh38.p14'),
    'gnomad_genome': ('hg38/gnomad/genomes/gnomad_genome_4.0_lite.vcf.gz', '4.0-institutional-lite'),
    'gnomad_exome': ('hg38/gnomad/exomes/gnomad_exome_4.0_lite.vcf.gz', '4.0-institutional-lite'),
    'clinvar': ('hg38/clinvar/clinvar.hg38.vcf.gz', clinvar),
    'clinvar_common_pathogenic': ('hg38/hgnc/clinvar_pathogenic_common.vcf.gz', rescue),
    'cosmic_targeted': ('hg38/cosmic/cosmic.targeted.v105.GRCh38.vcf.gz', 'v105'),
}
mapping = {name: {'path': str(root/key), 'uri': remote+'/'+key, 'release': release}
           for name, (key, release) in files.items()}
Path('review.resource-files.json').write_text(json.dumps(mapping, indent=2)+'\n')
PY

docker run --rm --platform linux/amd64 \
    -v "$PWD:$PWD" -v "$LOCAL_REFS:$LOCAL_REFS:ro" -w "$PWD" \
    --entrypoint python3 "$REVIEW_IMAGE" \
    bin/review_report.py manifest \
    --resources review.resource-files.json --output review.references.json

aws s3 cp review.references.json \
    "$ONCOSEQ_REF_ROOT/hg38/review/review.references.json" --only-show-errors
aws s3 ls "$ONCOSEQ_REF_ROOT/hg38/vep-cache/homo_sapiens/115_GRCh38/info.txt"
```

The clinical profile reads the bundled unpinned development mapping by default.
Select the strict manifest explicitly with `--review_references` once prepared.
Keep a copy with validation results.
The review preflight checks resource hashes, indexes and assembly compatibility
before basecalling. Do not overwrite these references with newer catalogue
releases after approving the manifest.

## Run with the staged resources

The institutional preset now points to the Oncoseq reference root. After all
copies finish, the clinical profile supplies the auxiliary cache, review
enablement and bundled development reference mapping automatically. Use the preset:

```bash
-profile aws,clinical \
-params-file assets/institutional_hg38.params.json
```

See [the complete AWS run example](review-report.md#run-on-aws-batch). The Batch
job role must be able to read the destination prefix.
