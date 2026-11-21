# Outstanding clinical functionality

Optic classification and its report tables are temporarily disabled because the pinned full Classy image (`sha-6a3f603`) does not declare Optic model assets. Restore the commented output declaration in `modules/local/classifier_process/main.nf` and the four report-channel connections in `subworkflows/local/methylation_analysis/classy.nf` after pinning and validating an image with Optic support. All other classifier tables remain required. Tracked in Beads issue `oncoseq-33w`.

Ask the primary Oncoseq authors to update their Dorado image and benchmark the
current release with SUP v5.2.0 and its matching modification model. Tracked in
Beads issue `oncoseq-4f7`.
