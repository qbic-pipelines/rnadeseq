# qbic-pipelines/rnadeseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 2.6 - dev

### Added

- [#267](https://github.com/qbic-pipelines/rnadeseq/pull/267) Add EnhancedVolcano plots and export to output results
- [#274](https://github.com/qbic-pipelines/rnadeseq/pull/274) Added workflow section to report with mermaid graph

### Changed

- [#267](https://github.com/qbic-pipelines/rnadeseq/pull/267) Export static and interactive plotly volcano plots
- [#268](https://github.com/qbic-pipelines/rnadeseq/pull/268) Apply EnhancedVolcano color code to plotly volcano plots
- [#269](https://github.com/qbic-pipelines/rnadeseq/pull/269) Put gprofiler manhattan plots in tabs
- [#274](https://github.com/qbic-pipelines/rnadeseq/pull/274) Changed report headers, citations, links

### Fixed

- [#271](https://github.com/qbic-pipelines/rnadeseq/pull/271) Fix colors in plotly volcano plots
- [#273](https://github.com/qbic-pipelines/rnadeseq/pull/272) Fix duplicate row.names in pathway analysis and max dimensions error
- [#274](https://github.com/qbic-pipelines/rnadeseq/pull/274) Fix colors volcano plotly, sample dist. heatmap sampleNames, header trim whitespaces

## 2.5 - The Potato Eaters

### Added

- [#256](https://github.com/qbic-pipelines/rnadeseq/pull/256) Add trycatch to pathway enrichment plots so they are skipped when too large instead of throwing an error
- [#255](https://github.com/qbic-pipelines/rnadeseq/pull/255) Add usage docu for datasources, heatmaps_cluster_rows/cols and pathway_adj_pval_threshold params
- [#251](https://github.com/qbic-pipelines/rnadeseq/pull/251) Get raw gene count tables from either Salmon and RSEM analysis
- [#250](https://github.com/qbic-pipelines/rnadeseq/pull/250) Added clearer error message for incorrect contrast_pairs

### Changed

- [#260](https://github.com/qbic-pipelines/rnadeseq/pull/260) Release 2.5
- [#259](https://github.com/qbic-pipelines/rnadeseq/pull/259) Bump versions for release 2.5

### Fixed

- [#258](https://github.com/qbic-pipelines/rnadeseq/pull/258) Fixed some comments for release (removed excess checks for pathway_adj_pval_threshold, added default explanation of that param to Execute_report.R, fixed some whitespace)
- [#252](https://github.com/qbic-pipelines/rnadeseq/pull/252) Fixed github CI bug by updating actions/upload-artifact
- [#250](https://github.com/qbic-pipelines/rnadeseq/pull/250) Fixed incorrect reading and indexing of contrast_pairs

## 2.4 - A Pair of Shoes

### Added

- [#253](https://github.com/qbic-pipelines/rnadeseq/pull/253) Added separate param for adjusted p-value threshold for gprofiler
- [#245](https://github.com/qbic-pipelines/rnadeseq/pull/245) Added background gene list to pathway analysis output

### Changed

- [#249](https://github.com/qbic-pipelines/rnadeseq/pull/249) Release 2.4
- [#248](https://github.com/qbic-pipelines/rnadeseq/pull/248) Version bumps for release 2.4
- [#247](https://github.com/qbic-pipelines/rnadeseq/pull/247) Changed the report text that shows when the gprofiler GMT file could not be downloaded

### Removed

- [#247](https://github.com/qbic-pipelines/rnadeseq/pull/247) Removed hard-coded comment about usage of gprofiler databases KEGG and REAC in the report

## 2.3 - Flowering Orchards

### Added

- [#237](https://github.com/qbic-pipelines/rnadeseq/pull/237) Added Euclidian distance/z-score info to report text/heatmaps; added tables next to heatmaps in case the heatmap plotting does not work
- [#229](https://github.com/qbic-pipelines/rnadeseq/pull/229) Added param for clustering (or not) the heatmaps
- [#226](https://github.com/qbic-pipelines/rnadeseq/pull/226) Added logic to read newer multiqc files from rnaseq 3.12
- [#225](https://github.com/qbic-pipelines/rnadeseq/pull/225) Added param for pathway analysis datasources
- [#221](https://github.com/qbic-pipelines/rnadeseq/pull/221) Added padj to volcano hovertext

### Changed

- [#244](https://github.com/qbic-pipelines/rnadeseq/pull/244) Release 2.3 to master, part II (found a bug before actual release)
- [#237](https://github.com/qbic-pipelines/rnadeseq/pull/237) Changed geom_jitter to geom_point for volcano plot so that points are accurate
- [#233](https://github.com/qbic-pipelines/rnadeseq/pull/233) Release 2.3 to master, part I
- [#225](https://github.com/qbic-pipelines/rnadeseq/pull/225) Shortened names of some output files

### Fixed

- [#243](https://github.com/qbic-pipelines/rnadeseq/pull/243) Fix missing QBiC logo and css styling in report by adding logo/css params
- [#242](https://github.com/qbic-pipelines/rnadeseq/pull/242) Fix path issue in report by adding report_file/references_file params
- [#241](https://github.com/qbic-pipelines/rnadeseq/pull/241) Undo fix branch protection bug in order to get precommit test to run
- [#240](https://github.com/qbic-pipelines/rnadeseq/pull/240) Fix branch protection bug --> third attempt
- [#239](https://github.com/qbic-pipelines/rnadeseq/pull/239) Fix branch protection bug --> second attempt
- [#238](https://github.com/qbic-pipelines/rnadeseq/pull/238) Fix branch protection bug
- [#237](https://github.com/qbic-pipelines/rnadeseq/pull/237) Fix version in usage docs and remove duplicate datasources in test_custom_gmt.config
- [#236](https://github.com/qbic-pipelines/rnadeseq/pull/236) Fixed new multiqc check (in case both the files of the old and new mqc version are present)
- [#234](https://github.com/qbic-pipelines/rnadeseq/pull/234) Fixed unnecessary file permission changes from #225
- [#228](https://github.com/qbic-pipelines/rnadeseq/pull/228) Fixed text in report
- [#229](https://github.com/qbic-pipelines/rnadeseq/pull/229) Fixed cutoff enrichment plot labels, fixed wrong plotMA function being called (also fixed this changelog)
- [#225](https://github.com/qbic-pipelines/rnadeseq/pull/225) Fixed too many devices error from tryCatch around normalized heatmaps
- [#221](https://github.com/qbic-pipelines/rnadeseq/pull/221) Fixed non-conformable arrays bug, fix wrong volcano colors when no DE genes

## 2.2 - Avenue of Poplars

### Added

- [#215](https://github.com/qbic-pipelines/rnadeseq/pull/215) Added gene_name to boxplot titles and filenames; increased threshold before overlapping PCA labels are hidden
- [#213](https://github.com/qbic-pipelines/rnadeseq/pull/213) Added smrnaseq input support
- [#212](https://github.com/qbic-pipelines/rnadeseq/pull/212) Added computational methods if no --software_versions
- [#206](https://github.com/qbic-pipelines/rnadeseq/pull/206) Added logic to decide between rlog and vst, added tryCatch for heatmap saving because this only works unreliably
- [#202](https://github.com/qbic-pipelines/rnadeseq/pull/202) Added background list to pathway analysis
- [#197](https://github.com/qbic-pipelines/rnadeseq/pull/197) Added gprofiler version string to report
- [#196](https://github.com/qbic-pipelines/rnadeseq/pull/196) Added optional `--quote` parameter

### Changed

- [#219](https://github.com/qbic-pipelines/rnadeseq/pull/219) Release 2.2 to master
- [#220](https://github.com/qbic-pipelines/rnadeseq/pull/220) Commit suggestions from PR review for release 2.2
- [#218](https://github.com/qbic-pipelines/rnadeseq/pull/218) Preparing release 2.2 with version bumps
- [#217](https://github.com/qbic-pipelines/rnadeseq/pull/217) For rsem/salmon-imports, when a GTF is sometimes missing gene_names, these are now replaced by gene_ids instead
- [#215](https://github.com/qbic-pipelines/rnadeseq/pull/215) Boxplots of genes are now generated from rlog/vst counts instead of raw counts. Also, if batch effect correction is enabled, boxplots will be generated before and after correction
- [#214](https://github.com/qbic-pipelines/rnadeseq/pull/214) Added smrnaseq input support -->follow-up to check that the CI test works. Modified smrnaseq testdata to be like QBiC datasets. Updated container env and fixed a resulting bug with include_graphics. Heatmaps are now equally ordered in zip and report
- [#211](https://github.com/qbic-pipelines/rnadeseq/pull/21) Replaced heatmaply with pheatmap for static plots and removed kaleido and reticulate from container
- [#209](https://github.com/qbic-pipelines/rnadeseq/pull/209) Template update to 2.9, Chromium Falcon; exchanged file.exists for nf-core validation checks; changed round_DE param to int
- [#206](https://github.com/qbic-pipelines/rnadeseq/pull/206) Changed < and > to <= and => for logF/pval comparisons, renamed gene_counts_tables/deseq2_library_scaled_gene_counts.tsv to deseq2_library_scaled_gene_counts.tsv and added entry to the folder explanation; renamed param pval_threshold to adj_pval_threshold
- [#205](https://github.com/qbic-pipelines/rnadeseq/pull/205) Template update
- [#204](https://github.com/qbic-pipelines/rnadeseq/pull/204) Changed relevel path in test_relevel.config to the qbic-pipelines repo
- [#203](https://github.com/qbic-pipelines/rnadeseq/pull/203) Switched from Dockerhub to GHCR
- [#200](https://github.com/qbic-pipelines/rnadeseq/pull/200) Made software_versions optional
- [#198](https://github.com/qbic-pipelines/rnadeseq/pull/198) Changed heatmaps to scale in size automatically

### Fixed

- [#212](https://github.com/qbic-pipelines/rnadeseq/pull/212) Fixed movability of interactive gostplots
- [#208](https://github.com/qbic-pipelines/rnadeseq/pull/208) Fixed relevel bug, the function should now finally work!
- [#207](https://github.com/qbic-pipelines/rnadeseq/pull/207) Fixed check of samples in counts vs metadata
- [#206](https://github.com/qbic-pipelines/rnadeseq/pull/206) Added correct plot titles to meanSdPlot (depending on normalization)
- [#195](https://github.com/qbic-pipelines/rnadeseq/pull/195) Fixed section error in report

## 2.1 - Wheat Fields

### Added

- [#192](https://github.com/qbic-pipelines/rnadeseq/pull/192) Added software_version yml functionality
- [#191](https://github.com/qbic-pipelines/rnadeseq/pull/191) Added test software_versions.yml files for rsem and salmon
- [#188](https://github.com/qbic-pipelines/rnadeseq/pull/188) Added titles to static heatmaps, added labels to static PCA plots
- [#176](https://github.com/qbic-pipelines/rnadeseq/pull/176) Added output description to report
- [#175](https://github.com/qbic-pipelines/rnadeseq/pull/175) Added nf-core citation to report
- [#173](https://github.com/qbic-pipelines/rnadeseq/pull/173) Added GMT file to testdata dir to test #172
- [#172](https://github.com/qbic-pipelines/rnadeseq/pull/172) Added option to provide custom gost GMT, for online gost, GMT is downloaded
- [#151](https://github.com/qbic-pipelines/rnadeseq/pull/151) Added session info to report
- [#149](https://github.com/qbic-pipelines/rnadeseq/pull/149) Added gene names to PA tables
- [#148](https://github.com/qbic-pipelines/rnadeseq/pull/148) Added KEGG/REAC versions to report
- [#147](https://github.com/qbic-pipelines/rnadeseq/pull/147) Added check for contrast list/metadata comparison
- [#145](https://github.com/qbic-pipelines/rnadeseq/pull/145) Added pval threshold param
- [#136](https://github.com/qbic-pipelines/rnadeseq/pull/136) Added pytest checks and md5sums to make sure that output stays consistent

### Changed

- [#180](https://github.com/qbic-pipelines/rnadeseq/pull/180) Bump version to 2.1 in some more files
- [#179](https://github.com/qbic-pipelines/rnadeseq/pull/179) Release 2.1
- [#178](https://github.com/qbic-pipelines/rnadeseq/pull/178) Bump version to 2.1
- [#177](https://github.com/qbic-pipelines/rnadeseq/pull/177) Updated usage documentation
- [#174](https://github.com/qbic-pipelines/rnadeseq/pull/174) Template update, changed param --metadata to --input
- [#169](https://github.com/qbic-pipelines/rnadeseq/pull/169) Changed skip_pathway_analysis to run_pathway_analysis, default false
- [#165](https://github.com/qbic-pipelines/rnadeseq/pull/165) -fw entries in multiqc stats are now merged
- [#164](https://github.com/qbic-pipelines/rnadeseq/pull/164) Boxplots are now only generated for contrasts in list/matrix file if provided
- [#163](https://github.com/qbic-pipelines/rnadeseq/pull/163) Template update, re-added limma, annotationdbi, colorbrewr to env (were previously incorrectly deleted), switched container to mamba
- [#159](https://github.com/qbic-pipelines/rnadeseq/pull/159) Changed error messages for non-existing rsem/salmon files
- [#151](https://github.com/qbic-pipelines/rnadeseq/pull/151) PCA plots and heatmap are now interactive, volcano and enrichment plots are cleaned up in their layout
- [#145](https://github.com/qbic-pipelines/rnadeseq/pull/145) Renamed versions to software_versions

### Fixed

- [#188](https://github.com/qbic-pipelines/rnadeseq/pull/188) Fixed cut-off enrichment legends and cut-off volcano ylabs
- [#181](https://github.com/qbic-pipelines/rnadeseq/pull/181) Fixed nf-core version test in ci.yml, updated schema.yml
- [#167](https://github.com/qbic-pipelines/rnadeseq/pull/167) Corrected 5sum for batcheffect after container update
- [#152](https://github.com/qbic-pipelines/rnadeseq/pull/152) Fixed empty qlist for gost query

### Removed

- [#150](https://github.com/qbic-pipelines/rnadeseq/pull/150) Removed human and mouse db from env.yml, both are now installed during pipeline execution if needed

## 2.0.1 - Olive Trees hotfix I

### Added

- [#134](https://github.com/qbic-pipelines/rnadeseq/pull/134) Corrected some versions
- [#132](https://github.com/qbic-pipelines/rnadeseq/pull/132) Bump version to 2.0.1
- [#131](https://github.com/qbic-pipelines/rnadeseq/pull/131) Added design_batcheffect.txt
- [#130](https://github.com/qbic-pipelines/rnadeseq/pull/130) Added test_batcheffect to github tests

### Changed

- [#131](https://github.com/qbic-pipelines/rnadeseq/pull/131) Changed Sample_preparations.tsv by adding batch column

### Fixed

- [#130](https://github.com/qbic-pipelines/rnadeseq/pull/130) Fixed batch effect bug

## 2.0 - Olive Trees

### Added

- [#128](https://github.com/qbic-pipelines/rnadeseq/pull/128) Bump versions to 2.0
- [#125](https://github.com/qbic-pipelines/rnadeseq/pull/125) Added test_relevel to github tests, added explanation to usage.md that --species is not necessary if skipping pathway analysis
- [#123](https://github.com/qbic-pipelines/rnadeseq/pull/123) Export report volcano plots as SVG; save all plots additionally as PDF
- [#122](https://github.com/qbic-pipelines/rnadeseq/pull/122) Add searchable/sortable tables to report
- [#118](https://github.com/qbic-pipelines/rnadeseq/pull/118) Add parameter "--input_type" (and change --rawcounts to --gene_counts) to process featurecounts, rsem and salmon output from the new rnaseq; add igenomes.config to process different species
- [#104](https://github.com/qbic-pipelines/rnadeseq/pull/104) Add parameter "--skip_pathway_analysis"
- Bump versions to 1.4.0dev
- Add parameter "--min_DE_genes"
- [#97](https://github.com/qbic-pipelines/rnadeseq/pull/97) Update pipeline to DSL2
- [#107](https://github.com/qbic-pipelines/rnadeseq/pull/107) Add parameter "--skip_rlog"
- [#111](https://github.com/qbic-pipelines/rnadeseq/pull/111) Added enhanced volcano plots
- [#93](https://github.com/qbic-pipelines/rnadeseq/pull/93/) Add parameter "--nsubgenes"

### Changed

- [#115](https://github.com/qbic-pipelines/rnadeseq/pull/115) Template update
- [#117](https://github.com/qbic-pipelines/rnadeseq/pull/117) Turned LabID optional for report output in rnadeseq_report.Rmd
- Removed assets/report_options.yml
- [#110](https://github.com/qbic-pipelines/rnadeseq/pull/110) Changed report to use rlog normalization by default, vst is used if --skip_rlog is enabled

### Fixed

- [#127](https://github.com/qbic-pipelines/rnadeseq/pull/127) Allgenes files are not introduced in the PA report section anymore except for volcano plots
- [#126](https://github.com/qbic-pipelines/rnadeseq/pull/126) Allgenes files are not published in results anymore. Intermediate results are not zipped and published anymore
- [#125](https://github.com/qbic-pipelines/rnadeseq/pull/125) Fixed relevel bug
- [#124](https://github.com/qbic-pipelines/rnadeseq/pull/124) Combined the different scripts into the report. Manually added the changes from #126 and 127 to the 1script branch as at least the code part is combined from multiple previous scrips into a single one
- [#118](https://github.com/qbic-pipelines/rnadeseq/pull/118) Removed blacklist parameter and config and instead added trycatch to ignore pathways with errors
- [#105](https://github.com/qbic-pipelines/rnadeseq/pull/105) Fixed relevel and added test_relevel.config
- [#106](https://github.com/qbic-pipelines/rnadeseq/pull/106) Fixed `--logFCthreshold` bug
- [#108](https://github.com/qbic-pipelines/rnadeseq/pull/108) Fixed blacklist file not working
- [#88](https://github.com/qbic-pipelines/rnadeseq/issues/88) Fixed volcano plot axis

## 1.3.2 - Almond Blossoms hotfix II

### Added

- Bump versions to 1.3.2
- write deseq2 table to file

### Fixed

- Contrast names in report plots
- Convert species name to lower case also in report
- LogFC is also reported in the report and set in volcano plots

## 1.3.1 - Almond Blossoms hotfix

### Added

- Bump versions to 1.3.1

### Fixed

- Fix bug plots requested boxplots

### Changed

## 1.3.0 - Almond Blossoms

### Added

- Bump versions to 1.3.0dev
- [#74](https://github.com/qbic-pipelines/rnadeseq/pull/74) Add option to provide KEGG pathway blacklist
- [#74](https://github.com/qbic-pipelines/rnadeseq/pull/74) Make quote param optional
- [#74](https://github.com/qbic-pipelines/rnadeseq/pull/74) Make report options optional (default in assets)
- [#74](https://github.com/qbic-pipelines/rnadeseq/pull/74) DE gene and pathway summary table

### Fixed

- [#75](https://github.com/qbic-pipelines/rnadeseq/pull/75) Boxplot of normalized counts can also be done from non-DE genes.
- [#75](https://github.com/qbic-pipelines/rnadeseq/pull/75) More comprehensive variable names and comments

### Changed

- [#74](https://github.com/qbic-pipelines/rnadeseq/pull/74) Pathway analysis only perfomed if at least 2 DE genes

## 1.2.0 - starry night [11-09-2020]

### Added

- Added report_options.yml in assets/.

### Fixed

- Skipping pathway analysis for contrasts with no found DE genes.
- Fixed report pvalue typo.

## 1.1.0 - sunrise

### Added

- Major changes in handling contrasts.
- Major improvements to report.

## 1.0.0 - candlelight

Initial pre-release of qbic-pipelines/rnadeseq, created with the [nf-core](http://nf-co.re/) template.
