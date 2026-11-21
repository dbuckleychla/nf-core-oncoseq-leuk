# Reusing institutional S3 reference data

[Copy and upload the required references](stage-review-references.md) before
using the preset, which now points to the dedicated Oncoseq reference root.
The optional molecular review additionally consumes the indexed annotation
VCFs through a checksum-pinned manifest; see [review setup](review-report.md).

`assets/institutional_hg38.params.json` maps compatible inputs from FrankONTstein's `assets/aws_batch_references.yaml`. The source inventory records objects beneath `s3://plm-cpm-sandbox-research-usw2/ONT/Oncoseq_refdata/`. This session cannot list that prefix with its current AWS identity, so current object presence and contents remain unverified.

| Existing asset                                                   | Reuse in Oncoseq                                                                                                                                                 |
| ---------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `hg38/fasta/hg38.analysisSet.fa` and adjacent `.fai`             | Set `ref_cache` to the exact FASTA. Validate contig names/lengths against caller and annotation resources.                                                       |
| GENCODE v50 enrichment BED with 10000-bp buffer                  | Set `bed` and `padding = 10000`; the upstream workflow derives targets by removing padding. This may differ from FrankONTstein's explicit unpadded target BED.   |
| `dorado-models/` with SUP v5.2.0 and matching CpG v2 directories | Set `dorado_model_cache` to reuse the existing prefix without moving it beneath the FASTA directory.                                                             |
| VEP cache `homo_sapiens/115_GRCh38`                              | The preset uses the institutional `s3://plm-cpm-sandbox-research-usw2/ONT/Oncoseq_refdata/hg38/vep-cache/` cache for the pinned VEP image.                       |
| NASVAR repeats, MAF sites and GFF                                | Oncoseq already uses its image-bundled versions. Keep those defaults until external-content parity is established.                                               |
| Delly mappability and ichorCNA assets                            | Reuse Oncoseq's bundled resources after checking the exact FASTA/bin-size pairing. The legacy ichor 1-Mb assets cannot substitute for Oncoseq's 500-kb defaults. |
| External Lamprey/ALMA3/OPTiC models                              | Classy already bundles these classifier families. Validate versions and preprocessing before substituting external directories.                                  |
| Clair3 `r1041_e82_400bps_sup_v520_with_mv`                       | The current wrapper selects bundled v500 weights and omits move emission. External model/move support remains tracked in `oncoseq-abz`.                          |
| gnomAD/ClinVar VCFs                                              | Useful for later annotation parity, but the current VEP/filter workflow does not consume FrankONTstein's external annotation catalog directly.                   |

Clinical execution also requires the pinned auxiliary files used by Oncoseq: PharmCAT positions and CSI, the assembly-specific pbsv tandem-repeat BED, and the Severus PoN. Their source basenames must exist in `reference_aux_cache`, or beside the FASTA when that parameter is absent. The FrankONTstein checked-in S3 reference mapping does not identify compatible staged copies of those files. Use the fixed source URLs recorded by this fork to populate the cache and record checksums/object versions.

The institutional VEP 115 cache root is `s3://plm-cpm-sandbox-research-usw2/ONT/Oncoseq_refdata/hg38/vep-cache/`. The cache is mirrored from the public `s3://annotation-cache/vep_cache/115_GRCh38/` source, whose `homo_sapiens/115_GRCh38/info.txt` was verified anonymously on 20 November 2026. The metadata identifies GRCh38.p14 and GENCODE 49. The preset uses the mirror root above, which must contain `homo_sapiens/115_GRCh38/info.txt`. Reads use the Batch job role and its existing institutional reference-bucket permissions; anonymous access is unnecessary. A new cache build is unnecessary. Complete cache content checksums and clinical acceptance remain required. The preset leaves `input` and `reference_aux_cache` for the operator to supply.

`assets/reference_aux_hg38.urls` lists the four fixed source URLs needed for the auxiliary cache:

```bash
set -euo pipefail
mkdir -p aux
while IFS= read -r url; do
  curl --fail --location --retry 3 "$url" -o "aux/${url##*/}"
done < assets/reference_aux_hg38.urls
sha256sum aux/*.bgz aux/*.csi aux/*.bed aux/*.gz > aux/SHA256SUMS
aws s3 sync aux/ s3://plm-cpm-sandbox-research-usw2/ONT/Oncoseq_refdata/hg38/aux/
```

```bash
NXF_VER=26.04.6 nextflow run . -profile aws,clinical \
  -params-file assets/institutional_hg38.params.json \
  --input validation.samplesheet.csv \
  --reference_aux_cache s3://plm-cpm-sandbox-research-usw2/ONT/Oncoseq_refdata/hg38/aux/ \
  --run_id validation-001
```

Load the Terraform environment first. Queue/job-role parameters, S3 read permissions, and immutable release content checks remain required. The preset uses a SUP v5.2.0 basecaller and does not establish matching Clair3 model support or biological concordance.

ClairS-TO now receives `--indel_min_af 0.05`. DeepSomatic receives `--make_examples_extra_args=vsc_min_fraction_indels=0.05`. These candidate-generation defaults match FrankONTstein and can be selected explicitly with `clairsto_indel_min_af` and `deepsomatic_indel_min_af`. They change indel sensitivity and require representative validation; they are separate from report/display filters.
