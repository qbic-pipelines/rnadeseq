# qbic-pipelines/rnadeseq: Changelog

## 1.3.0 - dev

### Added

- Bump versions to 1.3.0dev
- Add option to provide KEGG pathway blacklist
- Make quote param optional
- Make report options optional (default in assets)

### Fixed

### Changed

- Pathway analysis only perfomed if at least 2 DE genes

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
