# qbic-pipelines/rnadeseq: Usage

## Table of contents

<!-- Install Atom plugin markdown-toc-auto for this ToC to auto-update on save -->
<!-- TOC START min:2 max:3 link:true asterisk:true update:true -->

- [qbic-pipelines/rnadeseq: Usage](#qbic-pipelinesrnadeseq-usage)
- [Table of contents](#table-of-contents)
- [Introduction](#introduction)
- [Pre-requisites](#pre-requisites)
- [Running the pipeline](#running-the-pipeline)
  - [Updating the pipeline](#updating-the-pipeline)
  - [Testing the pipeline](#testing-the-pipeline)
  - [Reproducibility](#reproducibility)
- [Mandatory arguments](#mandatory-arguments)
  - [`--gene_counts`](#--gene_counts)
  - [`--input_type`](#--input_type)
  - [`--input`](#--input)
  - [`--model`](#--model)
- [Contrasts](#contrasts)
  - [Default](#default)
  - [`--relevel`](#--relevel)
  - [`--contrast_matrix`](#--contrast_matrix)
  - [`--contrast_list`](#--contrast_list)
  - [`--contrast_pairs`](#--contrast_pairs)
- [Optional arguments](#optional-arguments)
  - [`--logFC_threshold`](#--logFC_threshold)
  - [`--adj_pval_threshold`](#--adj_pval_threshold)
  - [`--genelist`](#--genelist)
  - [`--batch_effect`](#--batch_effect)
  - [`--min_DEG_pathway`](#--min_deg_pathway)
  - [`--norm_method`](#--norm_method)
  - [`--vst_genes_number`](#--vst_genes_number)
  - [`--round_DE`](#--round_DE)
  - [`--run_pathway_analysis`](#--run_pathway_analysis)
  - [`--pathway_adj_pval_threshold`](#--pathway_adj_pval_threshold)
  - [`--heatmaps_cluster_rows`](#--heatmaps_cluster_rows)
  - [`--heatmaps_cluster_cols`](#--heatmaps_cluster_cols)
  - [`--input_type`](#--input_type)
  - [`--multiqc`](#--multiqc)
  - [`--project_summary`](#--project_summary)
  - [`--software_versions`](#--software_versions)
  - [`--citest`](#--citest)
- [Reference genome options](#reference-genome-options)
  - [`--genome`](#--genome)
  - [`--gtf`](#--gtf)
  - [`--organism`](#--organism)
  - [`--species_library`](#--species_library)
  - [`--keytype`](#--keytype)
  - [`--custom_gmt`](#--custom_gmt)
  - [`--set_background`](#--set_background)
  - [`--custom_background`](#--custom_background)
  - [`--datasources`](#--datasources)
  - [`--igenomes_base`](#--igenomes_base)
  - [`--igenomes_ignore`](#--igenomes_ignore)
- [Job resources](#job-resources)
  - [Automatic resubmission](#automatic-resubmission)
  - [Custom resource requests](#custom-resource-requests)
- [Other command line parameters](#other-command-line-parameters)
  - [`--outdir`](#--outdir)
  - [`--email`](#--email)
  - [`-name`](#-name)
  - [`-resume`](#-resume)
  - [`-c`](#-c)
  - [`--custom_config_version`](#--custom_config_version)
  - [`--custom_config_base`](#--custom_config_base)
  - [`--max_memory`](#--max_memory)
  - [`--max_time`](#--max_time)
  - [`--max_cpus`](#--max_cpus)
  - [`--plaintext_email`](#--plaintext_email)
  - [`--monochrome_logs`](#--monochrome_logs)
  <!-- TOC END -->

## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

<!-- TODO qbic-pipelines: Document required command line parameters to run the pipeline-->

## Pre-requisites

The `qbic-pipelines/rnadeseq` pipeline relies on the output from the `nf-core/rnaseq` pipeline or the `nf-core/smrnaseq` pipeline. To be able to match the results of the initial pipeline with the metadata samplesheet containing the experimental design for the differential expression analysis, **the filenames of the fastq files used as input to the `qbic-pipelines/rnadeseq` pipeline need to start with the corresponding QBiC codes!**. _E.g. QBICKXXXXX_original_file_name.fastq_. Once the filenames are corrected if necessary, you can run the `qbic-pipelines/rnadeseq` pipeline as usual.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run qbic-pipelines/rnadeseq -r 2.4 -profile docker \
--gene_counts 'merged_gene_counts.txt' \
--input 'QXXXX_sample_preparations.tsv' \
--model 'linear_model.txt' \
--contrast_matrix 'contrasts.tsv' \
--project_summary 'QXXXX_summary.tsv' \
--multiqc 'MultiQC.zip' \
--software_versions 'software_versions.csv' \
--genome GRCh37
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
results             # Finished results (configurable, see below)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> ⚠️ Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run qbic-pipelines/rnadeseq -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
gene_counts: 'merged_gene_counts.txt'
model: 'linear_model.txt'
<...>
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull qbic-pipelines/rnadeseq
```

### Testing the pipeline

A number of test profiles are prepared to allow for easy execution of the pipeline with different parameters to check if these work. Two of these, test_star_rsem and test_star_salmon, do not work by simply running from qbic-pipelines as both require a param pointing to a folder, not a file, and these folders cannot be accessed via raw.githubusercontent.com. Instead you have to either clone the repo, change into the clone and run the pipeline locally, e.g. for the rsem profile, like this:

```bash
nextflow run . -profile docker,test_star_rsem
```

Or you have to download the testdata/QDESQ folder manually, like so:

```bash
curl https://codeload.github.com/qbic-pipelines/rnadeseq/tar.gz/master | tar -xz --strip=2 rnadeseq-master/testdata/QDESQ && mkdir testdata && mv QDESQ testdata
```

Afterwards, you should be able to also run test_star_rsem and test_star_salmon from qbic-pipelines/rnadeseq without manually cloning, e.g.:

```bash
nextflow run qbic-pipelines/rnadeseq -r 2.4 -profile docker,test_star_salmon
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [qbic-pipelines/rnadeseq releases page](https://github.com/qbic-pipelines/rnadeseq/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Mandatory arguments

<!-- TODO qbic-pipelines: Document required command line parameters -->

### `--gene_counts`

Gene counts. Can be:

- a raw count table (TSV), column names must start with the QBiC code, columns are samples and rows are genes, e.g.:

```bash
--gene_counts 'path/to/raw_count_table.tsv'
```

```tsv
Geneid  gene_name   QBICK00001_Sample1  QBICK00002_Sample2
ENSG00000000003  TSPAN6  150   3000
ENSG00000000005   TNMD    80  6
```

- OR a folder containing rsem output files (folder/sampleXXX.genes.results)
- OR a folder containing subfolders with salmon output (folder/sampleXXX/quant.sf).
  For rsem and salmon, the --metadata file "QBiC Code" column must provide the name of each sample (i.e. the respective folder/file name), and the --input_type parameter must be set to 'rsem' or 'salmon'.
  For example:

```bash
--gene_counts 'path/to/salmon_folder' --metadata 'path/to/salmon_metadata.tsv'

tree 'path/to/salmon_folder'
path/to/salmon_folder
├── QDESQ081AU
│   └── quant.sf
├── QDESQ082A4
│   └── quant.sf
├── QDESQ083AC
│   └── quant.sf
├── QDESQ084AK
│   └── quant.sf

head 'path/to/salmon_metadata.tsv'
```

```tsv
QBiC Code       Secondary Name  Condition: cellline
QDESQ081AU      Sample1 GM12878
QDESQ082A4      Sample2 GM12878
QDESQ083AC      Sample3 K562
QDESQ084AK      Sample4 K562
```

- OR a folder containing files with smrnaseq output (folder/sampleXXX_mature_hairpin.sorted.idxstats, folder/sampleXXX_mature_sorted.idxstats).
  For smrnaseq, the --input_type parameter must be set to 'smrnaseq'.
  For example:

```bash
--gene_counts 'path/to/smrnaseq_folder' --metadata 'path/to/smrnaseq_metadata.tsv'

tree testdata/smrnaseq/counts/
testdata/smrnaseq/counts/
├── Clone1_NN1_mature.sorted.idxstats
├── Clone1_NN1_mature_hairpin.sorted.idxstats
├── Clone1_NN3_mature.sorted.idxstats
├── Clone1_NN3_mature_hairpin.sorted.idxstats
├── Clone9_NN1_mature.sorted.idxstats
├── Clone9_NN1_mature_hairpin.sorted.idxstats
.
.
.

head 'path/to/smrnaseq_metadata.tsv'
```

```tsv
QBiC Code       Secondary Name  Condition: cellline
QDESQ081AU      Sample1 GM12878
QDESQ082A4      Sample2 GM12878
QDESQ083AC      Sample3 K562
QDESQ084AK      Sample4 K562
```

### `--input`

Metadata table/samplesheet (TSV) is the "Sample_preparations_sheet.tsv" that can be directly downloaded from the qPortal --> Browser. Rows are samples and columns contain sample grouping. Important columns are:

- **QBiC Code**: is needed to match metadata/samplesheet with the raw counts.
- **Secondary Name**, samples will be named with the pattern: QBiC code + Secondary name.
- **Condition: tag**: a separated column for each of the conditions. The headers of this columns start with "Condition: ". The values of these columns should not contain spaces.

```tsv
QBiC Code   Secondary Name  Condition: treatment
QBICK00001  Sample1 treated
QBICK00002  Sample2 untreated
```

### `--model`

Linear model function to calculate the contrasts (TXT). Variable names should be "condition_tag", where the tag matches the "Condition: tag" headers in the metadata/samplesheet file. E.g.

```txt
~ condition_genotype + condition_treatment
```

## Contrasts

Contrasts represent 2 or more conditions to be compared. One should compare the **Experiment** versus **Control** for correct differential expression analysis (and not the other way round). There are three different parameters that can be used to define contrasts, which are explained in the following sections. One or multiple contrast input files can be provided, if multiple are provided, the contrasts in the multiple files will be added to the report.

### Default

By default, DESeq2 will calculate pairwise contrasts given a linear model file. If you do not provide any contrast files, the differential gene expression analysis will be performed with the default contrasts as calculated by DESeq2. Try this option first, if you are unsure about your contrasts.

### `--relevel`

Sometimes condition factors have obvious levels. By default, the base level of a factor will be the first one when sorted alphabetically. In order to manually sort factors, to set them as the denominator for your contrasts, you can use this parameter to provide a tsv file for releveling the factor(s) for your contrast(s). An example input tsv file is shown here:

```tsv
factor  level
condition_genotype  wild_type
condition_treatment control

```

### `--contrast_matrix`

Table in tsv format indicating which contrasts to consider. Each contrast is specified in one column, the rows correspond to each of the expanded terms of the linear model. If you are unsure about how the linear model expanded terms look like, run the pipeline once without specifying contrasts, then the coefficient terms for the provided model will be stored under "differential_gene_expression/metadata/DESeq2_coefficients.tsv". An example input tsv file is shown here:

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

### `--contrast_pairs`

Table in tsv format indicating pairs of contrasts to consider. This is used to calculate interaction effects between contrasts. Each row corresponds to an interaction effect. The first column indicates the desired contrast name, the second column the first contrast in the numerator and the third column the contrast in the denominator, of the interaction.

```tsv
contrast_name contrast_numerator contrast_denominator
interaction_effect condition_treatment_treated_vs_control  condition_genotype_KO_vs_WT

```

## Optional arguments

### `--logFC_threshold`

Threshold (int) to apply to Log 2 Fold Change to consider a gene as differentially expressed. There is no threshold applied by default to Log2 Fold Change.

### `--adj_pval_threshold`

Adjusted p-value (float) to consider a gene as differentially expressed. The default value is 0.05.

### `--genelist`

List of genes of interest (one per line) for which additional plots are generated.
Heatmaps and boxplots use normalized counts across all samples.
If a genelist is provided, the volcano plots show only the labels of these requested genes.
The gene list should contain gene [HUGO](https://www.genenames.org/) symbols (no Ensemble IDs).

### `--batch_effect`

Option needed to account for batch effects in the data. To control for batch effects follow ALL these steps:

- Include the batch effect in the metadata file in a column with the header `batch`.
- Your design file needs to additionally include the batch effect in the linear model. E.g.:

  ```R
  ~ batch + condition_genotype
  ```

- Use the `--batch_effect` option when running the pipeline to generate an extra PCA plot with the corrected batch effects.

Then the DESeq2 script calculates the contrasts as usual, the batch effect just needs to be considered during the design definition.
For more information, please check the [DESeq2 vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).
Please note: Setting the `--batch_effect` option will NOT remove such effects from your data. The parameter is only used to visualize the effects by generating two sets of PCA and boxplots instead of only one; one set of plots will be generated from the input data, the other from the results of the batch correction in order to illustrate how the batches affect data.

### `--min_DEG_pathway`

Integer indicating how many genes in a pathway must be differentially expressed to be considered as enriched, and report these pathways in tables and the final report. The default value is 1.

### `--norm_method`

Set this parameter to either `vst`, `vst-force` or `rlog` (default) to control which transformation should be applied to the data. If the input data has too large size factor variances (>=0.05) the pipeline will override a user-specified `vst` transformation to use `rlog` instead. If this is enforced it will be stated in the report. If you still prefer to use the `vst` method regardless of the variances in the data, please use the `vst-force` parameter. While the `vst` transformation is much faster, the `rlog` is better suited to account for different sequencing depths between samples. Check here for more information on [count data transformations](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-data-transformations).

### `--vst_genes_number`

This is ignored if `--norm_method` is set to `rlog`. If using the `vst` transformation, consider using this parameter for small datasets and low numbers of genes, e.g. with small RNA-Seq data. The default number of genes for applying the `vst` function for varianceStabilizingTransformation in DESeq2 is 1000. For smaller datasets there will be an error. The solution is to reduce the number of genes to sample for the transformation ( < 1000 ). More information/solution here: [DESeq2 vst function error](https://www.biostars.org/p/456209/).

### `--round_DE`

Integer indicating to how many decimals to round the DE results (default: -1, indicating no rounding).

### `--run_pathway_analysis`

Set this flag to run pathway analysis, otherwise, this step will be skipped.

### `--pathway_adj_pval_threshold`

Use this param to specifically set the adjusted p value threshold for pathway enrichment analysis (will otherwise use the same adjusted p value threshold as for the DE analysis, as set with the param `--adj_pval_threshold`).

### `--heatmaps_cluster_rows`

Use this flag to set whether the heatmaps of gene expression in enriched pathways should by clustered row-wise (default true).

### `--heatmaps_cluster_cols`

Use this flag to set whether the heatmaps of gene expression in enriched pathways should by clustered column-wise (default false).

### `--input_type`

This tells the pipeline which type of input dataset is provided. Must be one of 'featurecounts', 'rsem', 'salmon', 'smrnaseq', default: featurecounts.

### `--multiqc`

Path to the MultiQC zipped folder containing the MultiQC plots generated by the RNAseq pipeline, the MultiQC html report, and the MultiQC raw data; optional.

### `--project_summary`

Project summary as downloaded from the portal: User database portlet, Projects tab, select your project and "Download Project Information". Please check first that this information is correct in the project Browser.

### `--software_versions`

Path to the `Software_versions.csv/.yml` file that is generated by the nf-core/rnaseq pipeline.

The CSV file should be tab-separated:

```bash
nf-core/rnaseq	v1.4dev
Nextflow	v19.01.0
FastQC	v0.11.8
Cutadapt	v2.1
Trim Galore!	v0.5.0
STAR	vSTAR_2.6.1d
```

Alternatively, the YML file should look like this:

```bash
BEDTOOLS_GENOMECOV:
  bedtools: 2.30.0
CUSTOM_DUMPSOFTWAREVERSIONS:
  python: 3.10.6
  yaml: '6.0'
CUSTOM_GETCHROMSIZES:
  custom: 1.15.1
```

### `--citest`

Developer param, tells the pipeline it is run on a github server which affects how output files are saved.

## Reference genome options

### `--genome`

Which genome to use for analysis, e.g. GRCh37; see `/conf/igenomes.config` for which genomes are available. When running the pipeline with rsem or salmon and/or with pathway analysis, this parameter is required unless you separately provide the parameters `--gtf` (if rsem/salmon), `--organism`, `--species_library` and `--keytype` (these three if pathway analysis is run). If your target genome has not been fully implemented (i.e. the entries for species_library, organism and keytype are missing), please open a new [issue](https://github.com/qbic-pipelines/rnadeseq/issues).

### `--gtf`

GTF file to be used for DESeq if input is rsem or salmon, not necessary for featurecounts or smrnaseq.

### `--organism`

Which organism name to use for pathway analysis, e.g. `hsapiens`, not necessary if `--run_pathway_analysis = false`.

### `--species_library`

Which bioconductor library to use for pathway analysis, e.g. org.Hs.eg.db, not necessary if `--run_pathway_analysis = false`.

### `--keytype`

Which keytype to use for pathway analysis, e.g. ENSEMBL, not necessary if `--run_pathway_analysis = false`.

### `--custom_gmt`

Path to custom GMT file to use instead of querying against the live gprofiler database, not necessary if `--run_pathway_analysis = false`.

### `--set_background`

Whether to restrict pathway analysis to a background gene list (default: true, will restrict to those genes with an expression > 0 in at least one sample), not necessary if `--run_pathway_analysis = false`.

### `--custom_background`

Path to custom background TXT file with one gene ID per line to use as background genes, not necessary if `--run_pathway_analysis = false` or `--set_background = false`.

### `--datasources`

Which datasources to use for pathway analysis, comma-separated string like 'KEGG,REAC'. See param 'sources' on https://rdrr.io/cran/gprofiler2/man/gost.html for a list of available sources. If not set, will use all sources. If set while a --custom_gmt is provided, will filter the GMT for these datasources (will not filter for the GO subtypes like GO:BP, just for GO).

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `qbic-pipelines` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://qbic-pipelines-invite.herokuapp.com/).

## Other command line parameters

<!-- TODO qbic-pipelines: Describe any other command line flags here -->

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default is set to `master`.

```bash
## Download and use config file with following git commit id
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

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data, therefore needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN). We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so based on the search results the file we want is `modules/nf-core/software/star/align/main.nf`. If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9). The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements. The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB. Providing you haven't set any other standard nf-core parameters to **cap** the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB. The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN' {
        memory = 100.GB
    }
}
```

> **NB:** We specify the full process name i.e. `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN` in the config file because this takes priority over the short name (`STAR_ALIGN`) and allows existing configuration using the full process name to be correctly overridden.
>
> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

### Updating containers

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon every time a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

   - For Docker:

   ```nextflow
   process {
       withName: PANGOLIN {
           container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
       }
   }
   ```

   - For Singularity:

   ```nextflow
   process {
       withName: PANGOLIN {
           container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
       }
   }
   ```

   - For Conda:

   ```nextflow
   process {
       withName: PANGOLIN {
           conda = 'bioconda::pangolin=3.0.5'
       }
   }
   ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
