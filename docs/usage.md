# qbic-pipelines/rnadeseq: Usage

## Table of contents

<!-- Install Atom plugin markdown-toc-auto for this ToC to auto-update on save -->
<!-- TOC START min:2 max:3 link:true asterisk:true update:true -->
* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
* [Pre-requisites](#pre-requisites)
* [Running the pipeline](#running-the-pipeline)
  * [Updating the pipeline](#updating-the-pipeline)
  * [Reproducibility](#reproducibility)
* [Mandatory arguments](#Mandatory-arguments)
  * [`-profile`](#-profile)
  * [`--rawcounts`](#--rawcounts)
  * [`--metadata`](#--metadata)
  * [`--model`](#--model)
  * [`--species`](#--species)
  * [`--project_summary`](#--project-summary)
  * [`--multiqc`](#--multiqc)
  * [`--versions`](#--versions)
  * [`--report_options`](#--report_options)
* [Contrasts](#contrasts)
  * [`default`](#default)
  * [`--relevel`](#--relevel)
  * [`--contrast_matrix`](#--contrast_matrix)
  * [`--contrast_list`](#--contrast_list)
  * [`--contrast_pairs`](#--contrast_pairs)
* [Optional arguments](#Optional-arguments)
  * [`--logFCthreshold`](#--logFCthreshold)
  * [`--genelist`](#--genelist)
  * [`--batch_effect`](#--batch_effect)
  * [`--quote`](#--quote)
* [Special cases](#Special-cases)
  * [Controlling for batch effects](#Controlling-for-batch-effects)
* [AWS Batch specific parameters](#aws-batch-specific-parameters)
  * [`--awsqueue`](#--awsqueue)
  * [`--awsregion`](#--awsregion)
* [Job resources](#job-resources)
  * [Automatic resubmission](#automatic-resubmission)
  * [Custom resource requests](#custom-resource-requests)
* [Other command line parameters](#other-command-line-parameters)
  * [`--outdir`](#--outdir)
  * [`--email`](#--email)
  * [`-name`](#-name)
  * [`-resume`](#-resume)
  * [`-c`](#-c)
  * [`--custom_config_version`](#--custom_config_version)
  * [`--custom_config_base`](#--custom_config_base)
  * [`--max_memory`](#--max_memory)
  * [`--max_time`](#--max_time)
  * [`--max_cpus`](#--max_cpus)
  * [`--plaintext_email`](#--plaintext_email)
  * [`--monochrome_logs`](#--monochrome_logs)
  * [`--multiqc_config`](#--multiqc_config)
<!-- TOC END -->

## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

<!-- TODO qbic-pipelines: Document required command line parameters to run the pipeline-->

## Pre-requisites

The `qbic-pipelines/rnadeseq` pipeline relies on the output from the `nf-core/rnaseq` pipeline. To be able to match the results of the `nf-core/rnaseq` pipeline with the metadata sheet containing the experimental design for the differential expression analysis, **the filenames of the fastq files used as input to the `qbic-pipelines/rnadeseq` pipeline, need to start with the corresponding QBiC codes!**. *E.g. QBICKXXXXX_original_file_name.fastq*. Once the filenames are corrected if necessary, you can run the `qbic-pipelines/rnadeseq` pipeline as usual.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run qbic-pipelines/rnadeseq -r 1.1.0 -profile docker \
--rawcounts 'merged_gene_counts.txt' \
--metadata 'QXXXX_sample_preparations.tsv' \
--model 'linear_model.txt' \
--contrast_matrix 'contrasts.tsv' \
--project_summary 'QXXXX_summary.tsv' \
--multiqc 'MultiQC.zip' \
--quote 'QXXXX_signed_offer.pdf' \
--versions 'software_versions.csv' \
--report_options 'report_options.yml' \
--species Hsapiens
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull qbic-pipelines/rnadeseq
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [qbic-pipelines/rnadeseq releases page](https://github.com/qbic-pipelines/rnadeseq/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Mandatory arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker, test` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

* `awsbatch`
  * A generic configuration profile to be used with AWS Batch.
* `conda`
  * A generic configuration profile to be used with [conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from dockerhub: [`nfcore/rnadeseq`](http://hub.docker.com/r/nfcore/rnadeseq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub: [`nfcore/rnadeseq`](http://hub.docker.com/r/nfcore/rnadeseq/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

<!-- TODO qbic-pipelines: Document required command line parameters -->

### `--rawcounts`

Raw count table (TSV). Column names must start with the QBiC code. Columns are samples and rows are genes. For example:

```bash
--rawcounts 'path/to/raw_count_table.tsv'
```

### `--metadata`

Metadata table is the "Sample_preparations_sheet.tsv" that can be directly downloaded from the qPortal --> Browser. Rows are samples and columns contain sample grouping. Important columns are:

* **QBiC Code**: is needed to match metadata with the raw counts.
* **Secondary Name**, samples will be named with the pattern: QBiC code + Secondary name.
* **Condition: tag**: a separated column for each of the conditions. The headers of this columns start with "Condition: ". The values of these columns should not contain spaces.

### `--model`

Linear model function to calculate the contrasts (TXT). Variable names should be "condition_tag", where the tag matches the "Condition: tag" headers in the metadata file. E.g.

```txt
~ condition_genotype + condition_treatment
```

### `--species`

Species name. For example: Hsapiens, Mmusculus. To include new species, please open an issue with the species full scientific name.

### `--project_summary`

Project summary as downloaded from the portal: User database portlet, Projects tab, select your project and "Download Project Information". Please check first that this information is correct in the project Browser.

### `--multiqc`

Path to the MultiQC zipped folder containing the MultiQC plots generated by the RNAseq pipeline, the MultiQC html report, and the MultiQC raw data.

### `--versions`

Path to the `Software_versions.csv` file generated by the RNAseq pipeline.

### `--report_options`

Configuration file describing the sections that need to be present in the report, text can also be specified that is desired to be added to the outlook section.
Check here an [example file](https://raw.githubusercontent.com/qbic-pipelines/rnadeseq/dev/testdata/report_options.yml).

## Contrasts

There are three different parameters that can be used to define contrasts, which are explained in the following sections. One or multiple contrast input files can be provided, if multiple are provided, the contrasts in the multiple files will be added to the report.

### Default

By default, DESeq2 will calculate some pairwise contrasts given the linear model file. If you do not provide any contrast files, the differential gene expression analysis will be performed with the default contrasts as calculated by DESeq2. Try this option first, if you are unsure about your contrasts.

### `--relevel`

Sometimes condition factors have obvious levels. By default, the base level of a factor will be the first one when sorted alphabetically. In order to manually sort factors, to set them as the denominator for your contrasts, you can use this parameter to provide a tsv file for releveling the factor(s) for your contrast(s). An example input tsv file is shown here:

```tsv
factor  level
condition_genotype  wild_type
condition_treatment control

```

### `--contrast_matrix`

Table in tsv format indicating which contrasts to consider. Each contrast is specified in one column, each row corresponds to the each of the expanded terms of the linear model. If you are unsure about how the linear model expanded terms look like, run the pipeline once without specifying contrasts, then the coefficient terms for the provided model will be stored under "differential_gene_expression/metadata/DESeq2_coefficients.tsv". An example input tsv file is shown here:

```tsv
coefficient treatment_treated_vs_control  genotype_KO_vs_WT
intercept 0 0
condition_treatment_treated 1 0
condition_genotype_KO 0  1

```

### `--contrast_list`

Table in tsv format indicating which contrasts to consider. Each contrast is specified in one row. The columns correspond to the factor, the level to be considered in the numerator and the level to be considered in the denominator of the contrast. For example:

```tsv
factor  numerator denominator
condition_treatment treated control
condition_genotype  KO  WT

```

### `--contrast_pairs``

Table in tsv format indicating pairs of contrasts to consider. This is used to calculate interaction effects between contrasts. Each row corresponds to an interaction effect. The first column indicates the desired contrast name, the second column the first contrast in the numerator and the third column the contrast in the denominator, of the interaction.

```tsv
contrast_name contrast_numerator contrast_denominator
interaction_effect condition_treatment_treated_vs_control  condition_genotype_KO_vs_WT

```

## Optional arguments

### `--logFCthreshold`

Threshold (int) to apply to Log 2 Fold Change to consider a gene as differentially expressed. There is no threshold applied by default to Log2 Fold Change.

### `--genelist`

List of genes (one per line) of which to plot heatmaps for normalized counts across all samples.

### `--batch_effect`

Option needed to account for batch effects in the data. Please check the section `Controlling for batch effects` to do so.

### `--quote`

Path to the signed copy of the QBiC offer as pdf, to be included in the report.

## Special cases

### Controlling for batch effects

To control for batch effects follow ALL these steps:

* Include the batch effect in the metadata file in a column with the header `batch`.
* Your design file needs to additionally include the batch effect in the linear model. E.g.:

  ```R
  ~ batch + condition_genotype
  ```

* Use the `--batch_effect` option when running the pipeline to generate an extra PCA plot with the corrected batch effects.

Then the DESeq2 script calculates the contrasts as usual, the batch effect just needs to be considered during the design definition.
For more information, please check the [DESeq2 vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

## AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use the `-awsbatch` profile and then specify all of the following parameters.

### `--awsqueue`

The JobQueue that you intend to use on AWS Batch.

### `--awsregion`

The AWS region to run your job in. Default is set to `eu-west-1` but can be adjusted to your needs.

The syntax for this reference configuration is as follows:

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `qbic-pipelines` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://qbic-pipelines-invite.herokuapp.com/).

## Other command line parameters

<!-- TODO qbic-pipelines: Describe any other command line flags here -->

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default is set to `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.

### `--multiqc_config`

Specify a path to a custom MultiQC configuration file.
