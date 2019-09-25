# qbicsoftware/rnadeseq: Usage

## Table of contents

<!-- Install Atom plugin markdown-toc-auto for this ToC to auto-update on save -->
<!-- TOC START min:2 max:3 link:true asterisk:true update:true -->
* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
* [Running the pipeline](#running-the-pipeline)
  * [Updating the pipeline](#updating-the-pipeline)
  * [Reproducibility](#reproducibility)
* [Mandatory arguments](#Mandatory-arguments)
  * [`-profile`](#-profile)
  * [`--rawcounts`](#--rawcounts)
  * [`--metadata`](#--metadata)
  * [`--model`](#--model)
  * [`--contrasts`](#--contrasts)
  * [`--species`](#--species)
  * [`--project_summary`](#--project-summary)
  * [`--multiqc`](#--multiqc)
  * [`--versions`](#--versions)
  * [`--quote`](#--quote)
  * [`--config`](#--config)
* [Optional arguments](#Optional-arguments)
  * [`--logFCthreshold`](#--logFCthreshold)
  * [`--genelist`](#--genelist)
  * [`--fastqc`](#--fastqc)
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

<!-- TODO qbicsoftware: Document required command line parameters to run the pipeline-->

## Running the pipeline
The typical command for running the pipeline is as follows:

```bash
nextflow run qbicsoftware/rnadeseq --reads '*_R{1,2}.fastq.gz' -profile docker
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
nextflow pull qbicsoftware/rnadeseq
```

### Reproducibility
It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [qbicsoftware/rnadeseq releases page](https://github.com/qbicsoftware/rnadeseq/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

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

<!-- TODO qbicsoftware: Document required command line parameters -->

### `--rawcounts`
Raw count table (TSV). Column names must start with the QBiC code. Columns are samples and rows are genes. For example:

```bash
--rawcounts 'path/to/raw_count_table.tsv'
```

Please note the following requirements:
1. The path must be enclosed in quotes

### `--metadata`
Metadata table is the "Sample_preparations_sheet.tsv" that can be directly downloaded from the qPortal --> Browser. Rows are samples and columns contain sample grouping. Important columns are:
* **QBiC Code**: is needed to match metadata with the raw counts.
* **Secondary Name**, samples will be named with the pattern: QBiC code + Secondary name.
* **Condition: xxxxx**: a separated column for each of the conditions.


```bash
--metadata 'path/to/sample_preparations.tsv'
```

### `--model`
Linear model function to calculate the contrasts (TXT). Variable names should be columns in metadata file.

### `--contrasts`
Table indicating which contrasts to consider. 1 or 0 for every variable specified in the design. Alternatively you can set the parameter "defaultcontrasts".

### `--species`
Species name. For example: Hsapiens, Mmusculus.

### `--project_summary`
Project summary as downloaded from the portal: User database portlet, Projects tab, select your project and "Download Project Information". Please check first that this information is correct in the project Browser.

### `--multiqc`
MultiQC zipped folder containing the MultiQC plots generated by the RNAseq pipeline, the MultiQC html report, and the MultiQC raw data.

### `--versions`
Software_versions.csv file generated by the RNAseq pipeline.

### `--quote`
Signed copy of the QBiC offer as pdf.

### `--config`
Configuration file describing the sections that need to be present in the report; it also contains text that one might want to add to the outlook section.

## Optional arguments

### `--logFCthreshold`
Threshold (int) to apply to Log 2 Fold Change to consider a gene as differentially expressed. There is no threshold applied by default to Log2 Fold Change.

### `--genelist`
List of genes (one per line) of which to plot heatmaps for normalized counts across all samples.

### `--fastqc`
Zipped folder containing the fastqc reports.

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

If you are likely to be running `qbicsoftware` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://qbicsoftware-invite.herokuapp.com/).

## Other command line parameters

<!-- TODO qbicsoftware: Describe any other command line flags here -->
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
