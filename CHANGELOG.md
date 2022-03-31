# qbic-pipelines/rnadeseq: Changelog

## 2.0.0 - dev

### Added

- [#104](https://github.com/qbic-pipelines/rnadeseq/pull/104) Add parameter "--skip_pathway_analysis"
- Bump versions to 1.4.0dev
- Add parameter "--min_DE_genes"
- [#97](https://github.com/qbic-pipelines/rnadeseq/pull/97) Update pipeline to DSL2
- [#107](https://github.com/qbic-pipelines/rnadeseq/pull/107) Add parameter "--skip_rlog"
- [#111]((https://github.com/qbic-pipelines/rnadeseq/pull/111) Added enhanced volcano plots
- [#93](https://github.com/qbic-pipelines/rnadeseq/pull/93/) Add parameter "--nsubgenes"


### Changed

- Removed assets/report_options.yml
- [#110](https://github.com/qbic-pipelines/rnadeseq/pull/110) Changed report to use rlog normalization by default, vst is used if --skip_rlog is enabled

### Fixed

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
