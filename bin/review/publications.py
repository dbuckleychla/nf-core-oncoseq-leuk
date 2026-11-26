"""Curated offline bibliography from upstream OncoSeq and DOI metadata.

Reviewed 2026-11-26. Publication links describe workflow capabilities, not
which optional analyses executed for a particular sample. See CITATIONS.md.
"""
from copy import deepcopy


PUBLICATIONS = {'schema_version': 1,
 'reviewed': '2026-11-26',
 'sources': ['https://github.com/chusj-pigu/nf-core-oncoseq/blob/main/CITATIONS.md',
             'https://github.com/chusj-pigu/nf-core-oncoseq/blob/main/README.md',
             'https://api.crossref.org'],
 'entries': [{'id': 'nf-core-oncoseq',
              'name': 'nf-core/oncoseq',
              'group': 'Workflows and infrastructure',
              'citation': 'Geoffrion, Nicholas et al. (2026). Accelerated and comprehensive pediatric cancer '
                          'profiling using whole-genome sequencing with adaptive sampling. Genome Medicine. '
                          'https://doi.org/10.1186/s13073-026-01815-7',
              'software_url': 'https://github.com/chusj-pigu/nf-core-oncoseq',
              'doi': '10.1186/s13073-026-01815-7',
              'kind': 'publication'},
             {'id': 'nf-core',
              'name': 'nf-core',
              'group': 'Workflows and infrastructure',
              'citation': 'Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso '
                          'P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. '
                          'Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: '
                          '32055031.',
              'software_url': 'https://pubmed.ncbi.nlm.nih.gov/32055031/',
              'doi': '10.1038/s41587-020-0439-x',
              'kind': 'publication'},
             {'id': 'nextflow',
              'name': 'Nextflow',
              'group': 'Workflows and infrastructure',
              'citation': 'Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables '
                          'reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: '
                          '10.1038/nbt.3820. PubMed PMID: 28398311.',
              'software_url': 'https://pubmed.ncbi.nlm.nih.gov/28398311/',
              'doi': '10.1038/nbt.3820',
              'kind': 'publication'},
             {'id': 'frankontstein',
              'name': 'FrankONTstein',
              'group': 'Workflows and infrastructure',
              'citation': 'Software and documentation reference.',
              'software_url': 'https://github.com/dbuckleychla/FrankONTstein',
              'doi': None,
              'kind': 'software'},
             {'id': 'chopper',
              'name': 'Chopper',
              'group': 'Read processing and QC',
              'citation': 'De Coster, W., & Rademakers, R. (2023). NanoPack2: Population-scale evaluation of '
                          'long-read sequencing data. Bioinformatics, 39(5), btad311. '
                          'https://doi.org/10.1093/bioinformatics/btad311',
              'software_url': 'https://github.com/wdecoster/chopper/',
              'doi': '10.1093/bioinformatics/btad311',
              'kind': 'publication'},
             {'id': 'cramino',
              'name': 'Cramino',
              'group': 'Read processing and QC',
              'citation': 'De Coster, W., & Rademakers, R. (2023). NanoPack2: Population-scale evaluation of '
                          'long-read sequencing data. Bioinformatics, 39(5), btad311. '
                          'https://doi.org/10.1093/bioinformatics/btad311',
              'software_url': 'https://github.com/wdecoster/cramino',
              'doi': '10.1093/bioinformatics/btad311',
              'kind': 'publication'},
             {'id': 'dorado',
              'name': 'Dorado',
              'group': 'Read processing and QC',
              'citation': 'Software and documentation reference.',
              'software_url': 'https://dorado-docs.readthedocs.io/en/latest/',
              'doi': None,
              'kind': 'software'},
             {'id': 'minimap2',
              'name': 'minimap2',
              'group': 'Read processing and QC',
              'citation': 'Li, H. (2021). New strategies to improve minimap2 alignment accuracy. Bioinformatics, '
                          '37(23), 4572–4574. https://doi.org/10.1093/bioinformatics/btab705',
              'software_url': 'https://github.com/lh3/minimap2',
              'doi': '10.1093/bioinformatics/btab705',
              'kind': 'publication'},
             {'id': 'modkit',
              'name': 'modkit',
              'group': 'Read processing and QC',
              'citation': 'Software and documentation reference.',
              'software_url': 'https://nanoporetech.github.io/modkit/',
              'doi': None,
              'kind': 'software'},
             {'id': 'mosdepth',
              'name': 'mosdepth',
              'group': 'Read processing and QC',
              'citation': 'Pedersen, B. S., & Quinlan, A. R. (2018). Mosdepth: Quick coverage calculation for '
                          'genomes and exomes. Bioinformatics (Oxford, England), 34(5), 867–868. '
                          'https://doi.org/10.1093/bioinformatics/btx699',
              'software_url': 'https://github.com/brentp/mosdepth',
              'doi': '10.1093/bioinformatics/btx699',
              'kind': 'publication'},
             {'id': 'ontime',
              'name': 'ontime',
              'group': 'Read processing and QC',
              'citation': 'Hall, M. (2024). mbhall88/ontime: v0.3.1 (0.3.1). Zenodo. '
                          'https://doi.org/10.5281/zenodo.10652351',
              'software_url': 'https://github.com/mbhall88/ontime',
              'doi': '10.5281/zenodo.10652351',
              'kind': 'software'},
             {'id': 'pysam',
              'name': 'pysam',
              'group': 'Read processing and QC',
              'citation': 'Software and documentation reference.',
              'software_url': 'https://pysam.readthedocs.io/en/latest/',
              'doi': None,
              'kind': 'software'},
             {'id': 'seqkit',
              'name': 'SeqKit',
              'group': 'Read processing and QC',
              'citation': 'Shen, W., Le, S., Li, Y., Hu, F. (2016). SeqKit: A Cross-Platform and Ultrafast '
                          'Toolkit for FASTA/Q File Manipulation. PLoS ONE, 11(10), e0163962. '
                          'https://doi.org/10.1371/journal.pone.0163962',
              'software_url': 'https://bioinf.shenwei.me/seqkit/',
              'doi': '10.1371/journal.pone.0163962',
              'kind': 'publication'},
             {'id': 'alma',
              'name': 'ALMA',
              'group': 'Methylation classifiers',
              'citation': 'Marchi, F., Shastri, V. M., Marrero, R. J., Nguyen, N. H. K., Öttl, A., Schade, A.-K., '
                          'Landwehr, M., Krali, O., Nordlund, J., Ghavami, M., Sckaff, F., Mansinghka, V. K., '
                          'Cao, X., Slayton, W., Starostik, P., Cogle, C. R., Ribeiro, R. C., Rubnitz, J. E., '
                          'Klco, J., … Lamba, J. K. (2025). Epigenomic diagnosis and prognosis of Acute Myeloid '
                          'Leukemia. Nature Communications, 16(1), 6961. '
                          'https://doi.org/10.1038/s41467-025-62005-4',
              'software_url': 'https://github.com/f-marchi/ALMA-classifier',
              'doi': '10.1038/s41467-025-62005-4',
              'kind': 'publication'},
             {'id': 'classy',
              'name': 'Classy',
              'group': 'Methylation classifiers',
              'citation': 'Software and documentation reference.',
              'software_url': 'https://github.com/chusj-pigu/wf-modules/tree/classy',
              'doi': None,
              'kind': 'software'},
             {'id': 'crossnn',
              'name': 'crossNN',
              'group': 'Methylation classifiers',
              'citation': 'Yuan, D., Jugas, R., Pokorna, P., Sterba, J., Slaby, O., Schmid, S., Siewert, C., '
                          'Osberg, B., Capper, D., Halldorsson, S., Vik-Mo, E. O., Zeiner, P. S., Weber, K. J., '
                          'Harter, P. N., Thomas, C., Albers, A., Rechsteiner, M., Reimann, R., Appelt, A., … '
                          'Euskirchen, P. (2025). crossNN is an explainable framework for cross-platform DNA '
                          'methylation-based classification of tumors. Nature Cancer, 6(7), 1283–1294. '
                          'https://doi.org/10.1038/s43018-025-00976-5',
              'software_url': 'https://gitlab.com/euskirchen-lab/crossNN',
              'doi': '10.1038/s43018-025-00976-5',
              'kind': 'publication'},
             {'id': 'marlin',
              'name': 'MARLIN',
              'group': 'Methylation classifiers',
              'citation': 'Steinicke, T. L., Benfatto, S., Capilla-Guerra, M. R., Monteleone, A. B., Young, J. '
                          'H., Shankar, S., Michaels, P. D., Tsai, H. K., Good, J. D., Kreso, A., van Galen, P., '
                          'Schliemann, C., Chen, E. C., Griffin, G. K., & Hovestadt, V. (2025). Rapid epigenomic '
                          'classification of acute leukemia. Nature Genetics, 57(10), 2456–2467. '
                          'https://doi.org/10.1038/s41588-025-02321-z',
              'software_url': 'https://github.com/hovestadt/MARLIN',
              'doi': '10.1038/s41588-025-02321-z',
              'kind': 'publication'},
             {'id': 'mpact',
              'name': 'MPACT',
              'group': 'Methylation classifiers',
              'citation': 'Smith, K. S., Fischer, T. T., Han, K., Kostecka, A., Lin, H., Senfter, D., Soliman, '
                          'T., Stepien, N., Volz, S., Schwarz, N., Wedig, T., Madlener, S., Haberler, C., Dhanda, '
                          'S. K., Upadhyaya, S. A., Blackburn, P. R., Schmook, M. T., de Bont, J., Haapasalo, H., '
                          '… Northcott, P. A. (2026). M-PACT leverages cell-free DNA methylomes to achieve robust '
                          'classification of pediatric brain tumors. Nature Cancer, 7(4), 667–683. '
                          'https://doi.org/10.1038/s43018-026-01115-4',
              'software_url': 'https://github.com/kylessmith/MethylVerse',
              'doi': '10.1038/s43018-026-01115-4',
              'kind': 'publication'},
             {'id': 'nanomix',
              'name': 'Nanomix',
              'group': 'Methylation classifiers',
              'citation': 'Software and documentation reference.',
              'software_url': 'https://github.com/Jonbroad15/nanomix',
              'doi': None,
              'kind': 'software'},
             {'id': 'optic',
              'name': 'OPTiC',
              'group': 'Methylation classifiers',
              'citation': 'Software and documentation reference.',
              'software_url': 'https://github.com/jwanglab/optic',
              'doi': None,
              'kind': 'software'},
             {'id': 'sturgeon',
              'name': 'Sturgeon',
              'group': 'Methylation classifiers',
              'citation': 'Vermeulen, C. et al. (2023). Ultra-fast deep-learned CNS tumour classification during '
                          'surgery. Nature. https://doi.org/10.1038/s41586-023-06615-2',
              'software_url': 'https://doi.org/10.1038/s41586-023-06615-2',
              'doi': '10.1038/s41586-023-06615-2',
              'kind': 'publication'},
             {'id': 'bcftools',
              'name': 'BCFtools',
              'group': 'Variant calling and phasing',
              'citation': 'Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., '
                          'Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years '
                          'of SAMtools and BCFtools. GigaScience, 10(2), giab008. '
                          'https://doi.org/10.1093/gigascience/giab008',
              'software_url': 'https://samtools.github.io/bcftools/bcftools.html',
              'doi': '10.1093/gigascience/giab008',
              'kind': 'publication'},
             {'id': 'clair3',
              'name': 'Clair3',
              'group': 'Variant calling and phasing',
              'citation': 'Zheng, Z., Li, S., Su, J., Leung, A. W.-S., Lam, T.-W., & Luo, R. (2022). Symphonizing '
                          'pileup and full-alignment for deep learning-based long-read variant calling. Nature '
                          'Computational Science, 2(12), 797–803. https://doi.org/10.1038/s43588-022-00387-x',
              'software_url': 'https://github.com/HKU-BAL/Clair3',
              'doi': '10.1038/s43588-022-00387-x',
              'kind': 'publication'},
             {'id': 'clairs-to',
              'name': 'ClairS-TO',
              'group': 'Variant calling and phasing',
              'citation': 'Chen, L., Zheng, Z., Su, J., Yu, X., Wong, A. O. K., Zhang, J., Lee, Y.-L., & Luo, R. '
                          '(2025). ClairS-TO: A deep-learning method for long-read tumor-only somatic small '
                          'variant calling. bioRxiv, 2025.03.10.642523. https://doi.org/10.1101/2025.03.10.642523',
              'software_url': 'https://github.com/HKU-BAL/ClairS-TO',
              'doi': '10.1101/2025.03.10.642523',
              'kind': 'preprint'},
             {'id': 'deepsomatic',
              'name': 'DeepSomatic',
              'group': 'Variant calling and phasing',
              'citation': 'Park, Jimin et al. (2025). Accurate somatic small variant discovery for multiple '
                          'sequencing technologies with DeepSomatic. Nature Biotechnology. '
                          'https://doi.org/10.1038/s41587-025-02839-x',
              'software_url': 'https://github.com/google/deepsomatic',
              'doi': '10.1038/s41587-025-02839-x',
              'kind': 'publication'},
             {'id': 'samtools',
              'name': 'SAMtools',
              'group': 'Variant calling and phasing',
              'citation': 'Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., '
                          'Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years '
                          'of SAMtools and BCFtools. GigaScience, 10(2), giab008. '
                          'https://doi.org/10.1093/gigascience/giab008',
              'software_url': 'https://github.com/samtools/samtools',
              'doi': '10.1093/gigascience/giab008',
              'kind': 'publication'},
             {'id': 'severus',
              'name': 'Severus',
              'group': 'Variant calling and phasing',
              'citation': 'Keskus, A. G., Bryant, A., Ahmad, T., Yoo, B., Aganezov, S., Goretsky, A., Donmez, A., '
                          'Lansdon, L. A., Rodriguez, I., Park, J., Liu, Y., Cui, X., Gardner, J., McNulty, B., '
                          'Sacco, S., Shetty, J., Zhao, Y., Tran, B., Narzisi, G., … Kolmogorov, M. (2026). '
                          'Severus detects somatic structural variation and complex rearrangements in cancer '
                          'genomes using long-read sequencing. Nature Biotechnology, 44(2), 247–257. '
                          'https://doi.org/10.1038/s41587-025-02618-8',
              'software_url': 'https://github.com/KolmogorovLab/Severus',
              'doi': '10.1038/s41587-025-02618-8',
              'kind': 'publication'},
             {'id': 'sniffles',
              'name': 'Sniffles',
              'group': 'Variant calling and phasing',
              'citation': 'Smolka, M., Paulin, L. F., Grochowski, C. M., Horner, D. W., Mahmoud, M., Behera, S., '
                          'Kalef-Ezra, E., Gandhi, M., Hong, K., Pehlivan, D., Scholz, S. W., Carvalho, C. M. B., '
                          'Proukakis, C., & Sedlazeck, F. J. (2024). Detection of mosaic and population-level '
                          'structural variants with Sniffles2. Nature Biotechnology, 42(10), 1571–1580. '
                          'https://doi.org/10.1038/s41587-023-02024-y',
              'software_url': 'https://github.com/fritzsedlazeck/Sniffles',
              'doi': '10.1038/s41587-023-02024-y',
              'kind': 'publication'},
             {'id': 'stellerator',
              'name': 'Stellerator',
              'group': 'Variant calling and phasing',
              'citation': 'Software and documentation reference.',
              'software_url': 'https://github.com/chusj-pigu/stellerator',
              'doi': None,
              'kind': 'software'},
             {'id': 'whatshap',
              'name': 'WhatsHap',
              'group': 'Variant calling and phasing',
              'citation': 'Martin, M., Patterson, M., Garg, S., O Fischer, S., Pisanti, N., Klau, G. W., '
                          'Schöenhuth, A., & Marschall, T. (2016). WhatsHap: Fast and accurate read-based '
                          'phasing. bioRxiv, 085050. https://doi.org/10.1101/085050',
              'software_url': 'https://whatshap.readthedocs.io/en/latest/index.html',
              'doi': '10.1101/085050',
              'kind': 'preprint'},
             {'id': 'ichorcna',
              'name': 'ichorCNA',
              'group': 'Copy-number analysis',
              'citation': 'Adalsteinsson, V. A., Ha, G., Freeman, S. S., Choudhury, A. D., Stover, D. G., '
                          'Parsons, H. A., Gydush, G., Reed, S. C., Rotem, D., Rhoades, J., Loginov, D., Livitz, '
                          'D., Rosebrock, D., Leshchiner, I., Kim, J., Stewart, C., Rosenberg, M., Francis, J. '
                          'M., Zhang, C.-Z., … Meyerson, M. (2017). Scalable whole-exome sequencing of cell-free '
                          'DNA reveals high concordance with metastatic tumors. Nature Communications, 8(1), '
                          '1324. https://doi.org/10.1038/s41467-017-00965-y',
              'software_url': 'https://github.com/broadinstitute/ichorCNA',
              'doi': '10.1038/s41467-017-00965-y',
              'kind': 'publication'},
             {'id': 'nasvar',
              'name': 'NASVAR',
              'group': 'Copy-number analysis',
              'citation': 'Software and documentation reference.',
              'software_url': 'https://github.com/jwanglab/nasvar',
              'doi': None,
              'kind': 'software'},
             {'id': 'qdnaseq',
              'name': 'QDNAseq',
              'group': 'Copy-number analysis',
              'citation': 'Scheinin, I., Sie, D., Bengtsson, H., van de Wiel, M. A., Olshen, A. B., van Thuijl, '
                          'H. F., van Essen, H. F., Eijk, P. P., Rustenburg, F., Meijer, G. A., Reijneveld, J. '
                          'C., Wesseling, P., Pinkel, D., Albertson, D. G., & Ylstra, B. (2014). DNA copy number '
                          'analysis of fresh and formalin-fixed specimens by shallow whole-genome sequencing with '
                          'identification and exclusion of problematic regions  in the genome assembly. Genome '
                          'Research, 24(12), 2022–2032. https://doi.org/10.1101/gr.175141.114',
              'software_url': 'https://www.bioconductor.org/packages/release/bioc/html/QDNAseq.html',
              'doi': '10.1101/gr.175141.114',
              'kind': 'preprint'},
             {'id': 'bedtools',
              'name': 'BEDTools',
              'group': 'Annotation and pharmacogenomics',
              'citation': 'Quinlan, A. R., & Hall, I. M. (2010). BEDTools: A flexible suite of utilities for '
                          'comparing genomic features. Bioinformatics, 26(6), 841–842. '
                          'https://doi.org/10.1093/bioinformatics/btq033',
              'software_url': 'https://bedtools.readthedocs.io/en/stable/',
              'doi': '10.1093/bioinformatics/btq033',
              'kind': 'publication'},
             {'id': 'ensembl-vep',
              'name': 'Ensembl VEP',
              'group': 'Annotation and pharmacogenomics',
              'citation': 'McLaren, W., Gil, L., Hunt, S. E., Riat, H. S., Ritchie, G. R. S., Thormann, A., '
                          'Flicek, P., & Cunningham, F. (2016). The Ensembl Variant Effect Predictor. Genome '
                          'Biology, 17(1), 122. https://doi.org/10.1186/s13059-016-0974-4',
              'software_url': 'https://www.ensembl.org/info/docs/tools/vep/index.html',
              'doi': '10.1186/s13059-016-0974-4',
              'kind': 'publication'},
             {'id': 'pharmcat',
              'name': 'PharmCAT',
              'group': 'Annotation and pharmacogenomics',
              'citation': 'Sangkuhl, K., Whirl-Carrillo, M., Whaley, R. M., Woon, M., Lavertu, A., Altman, R. B., '
                          'Carter, L., Verma, A., Ritchie, M. D., & Klein, T. E. (2020). Pharmacogenomics '
                          'Clinical Annotation Tool (PharmCAT). Clinical Pharmacology & Therapeutics, 107(1), '
                          '203–210. https://doi.org/10.1002/cpt.1568',
              'software_url': 'https://pharmcat.clinpgx.org/using/',
              'doi': '10.1002/cpt.1568',
              'kind': 'publication'},
             {'id': 'snpeff',
              'name': 'SnpEff',
              'group': 'Annotation and pharmacogenomics',
              'citation': 'Cingolani, Pablo et al. (2012). A program for annotating and predicting the effects of '
                          'single nucleotide polymorphisms, SnpEff. Fly. https://doi.org/10.4161/fly.19695',
              'software_url': 'https://pcingola.github.io/SnpEff/',
              'doi': '10.4161/fly.19695',
              'kind': 'publication'},
             {'id': 'clinvar',
              'name': 'ClinVar',
              'group': 'Reference resources',
              'citation': 'Landrum, Melissa J et al. (2017). ClinVar: improving access to variant interpretations '
                          'and supporting evidence. Nucleic Acids Research. https://doi.org/10.1093/nar/gkx1153',
              'software_url': 'https://www.ncbi.nlm.nih.gov/clinvar/',
              'doi': '10.1093/nar/gkx1153',
              'kind': 'publication'},
             {'id': 'cosmic',
              'name': 'COSMIC',
              'group': 'Reference resources',
              'citation': 'Sondka, Zbyslaw et al. (2023). COSMIC:\xa0a curated database of somatic variants and '
                          'clinical data for cancer. Nucleic Acids Research. https://doi.org/10.1093/nar/gkad986',
              'software_url': 'https://cancer.sanger.ac.uk/cosmic',
              'doi': '10.1093/nar/gkad986',
              'kind': 'publication'},
             {'id': 'gencode',
              'name': 'GENCODE',
              'group': 'Reference resources',
              'citation': 'Frankish, Adam et al. (2018). GENCODE reference annotation for the human and mouse '
                          'genomes. Nucleic Acids Research. https://doi.org/10.1093/nar/gky955',
              'software_url': 'https://www.gencodegenes.org',
              'doi': '10.1093/nar/gky955',
              'kind': 'publication'},
             {'id': 'gnomad',
              'name': 'gnomAD',
              'group': 'Reference resources',
              'citation': 'Chen, Siwei et al. (2023). A genomic mutational constraint map using variation in '
                          '76,156 human genomes. Nature. https://doi.org/10.1038/s41586-023-06045-0',
              'software_url': 'https://gnomad.broadinstitute.org',
              'doi': '10.1038/s41586-023-06045-0',
              'kind': 'publication'},
             {'id': 'figeno',
              'name': 'figeno',
              'group': 'Visualization and reporting',
              'citation': 'Sollier, E., Heilmann, J., Gerhauser, C., Scherer, M., Plass, C., & Lutsik, P. (2024). '
                          'Figeno: Multi-region genomic figures with long-read support. Bioinformatics, 40(6), '
                          'btae354. https://doi.org/10.1093/bioinformatics/btae354',
              'software_url': 'https://figeno.readthedocs.io/en/latest/',
              'doi': '10.1093/bioinformatics/btae354',
              'kind': 'publication'},
             {'id': 'imagemagick',
              'name': 'ImageMagick',
              'group': 'Visualization and reporting',
              'citation': 'ImageMagick Studio LLC. (2024). ImageMagick. Retrieved from https://imagemagick.org',
              'software_url': 'https://imagemagick.org/#gsc.tab=0',
              'doi': None,
              'kind': 'software'},
             {'id': 'quarto',
              'name': 'Quarto',
              'group': 'Visualization and reporting',
              'citation': 'Software and documentation reference.',
              'software_url': 'https://quarto.org/',
              'doi': None,
              'kind': 'software'},
             {'id': 'anaconda',
              'name': 'Anaconda',
              'group': 'Supporting software',
              'citation': 'Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. '
                          'Web.',
              'software_url': 'https://anaconda.com',
              'doi': None,
              'kind': 'software'},
             {'id': 'bioconda',
              'name': 'Bioconda',
              'group': 'Supporting software',
              'citation': 'Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster '
                          'J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for '
                          'the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: '
                          '10.1038/s41592-018-0046-7. PubMed PMID: 29967506.',
              'software_url': 'https://pubmed.ncbi.nlm.nih.gov/29967506/',
              'doi': '10.1038/s41592-018-0046-7',
              'kind': 'publication'},
             {'id': 'biocontainers',
              'name': 'BioContainers',
              'group': 'Supporting software',
              'citation': 'da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel '
                          'M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez '
                          'RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and '
                          'community-driven framework for software standardization. Bioinformatics. 2017 Aug '
                          '15;33(16):2580-2582. doi: 10.1093/bioinformatics/btx192. PubMed PMID: 28379341; PubMed '
                          'Central PMCID: PMC5870671.',
              'software_url': 'https://pubmed.ncbi.nlm.nih.gov/28379341/',
              'doi': '10.1093/bioinformatics/btx192',
              'kind': 'publication'},
             {'id': 'docker',
              'name': 'Docker',
              'group': 'Supporting software',
              'citation': 'Merkel, D. (2014). Docker: lightweight linux containers for consistent development and '
                          'deployment. Linux Journal, 2014(239), 2. doi: 10.5555/2600239.2600241.',
              'software_url': 'https://dl.acm.org/doi/10.5555/2600239.2600241',
              'doi': '10.5555/2600239.2600241',
              'kind': 'publication'},
             {'id': 'jsonlite',
              'name': 'jsonlite',
              'group': 'Supporting software',
              'citation': 'Ooms, J. (2014). The jsonlite Package: A Practical and Consistent Mapping Between JSON '
                          'Data and R Objects. arXiv:1403.2805 [stat.CO]. https://arxiv.org/abs/1403.2805.',
              'software_url': 'https://cran.r-project.org/web/packages/jsonlite/index.html',
              'doi': None,
              'kind': 'software'},
             {'id': 'optparse',
              'name': 'optparse',
              'group': 'Supporting software',
              'citation': 'Davis, T. L. (2024). optparse: Command Line Option Parser. '
                          'https://CRAN.R-project.org/package=optparse',
              'software_url': 'https://CRAN.R-project.org/package=optparse',
              'doi': None,
              'kind': 'software'},
             {'id': 'r',
              'name': 'R',
              'group': 'Supporting software',
              'citation': 'R Core Team. (2021). R: A Language and Environment for Statistical Computing. R '
                          'Foundation for Statistical Computing. https://www.R-project.org/',
              'software_url': 'https://www.r-project.org/about.html',
              'doi': None,
              'kind': 'software'},
             {'id': 'singularity',
              'name': 'Singularity',
              'group': 'Supporting software',
              'citation': 'Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of '
                          'compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. '
                          'eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.',
              'software_url': 'https://pubmed.ncbi.nlm.nih.gov/28494014/',
              'doi': '10.1371/journal.pone.0177459',
              'kind': 'publication'},
             {'id': 'tidyverse',
              'name': 'tidyverse',
              'group': 'Supporting software',
              'citation': 'Wickham, H., Averick, M., Bryan, J., Chang, W., McGowan, L. D., François, R., '
                          'Grolemund, G., Hayes, A., Henry, L., Hester, J., Kuhn, M., Pedersen, T. L., Miller, '
                          'E., Bache, S. M., Müller, K., Ooms, J., Robinson, D., Seidel, D. P., Spinu, V., … '
                          'Yutani, H. (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), '
                          '1686. https://doi.org/10.21105/joss.01686',
              'software_url': 'https://tidyverse.tidyverse.org/articles/paper.html',
              'doi': '10.21105/joss.01686',
              'kind': 'publication'}]}


def catalog():
    return deepcopy(PUBLICATIONS)
