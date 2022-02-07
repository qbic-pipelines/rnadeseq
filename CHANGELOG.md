# qbic-pipelines/rnadeseq: Changelog

## 1.4.0 - dev

### Added

- Bump versions to 1.4.0dev
- Add parameter "--min_DE_genes"

### Fixed

- Fixed `--logFCthreshold` issue [#100](https://github.com/qbic-pipelines/rnadeseq/issues/100)

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
