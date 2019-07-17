# qbicsoftware/rnadeseq: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

<!-- TODO qbicsoftware: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [DE analysis](#DE analysis) - differential expression analysis with DESeq2
* [Pathway analysis](#Pathway analysis) - pathway analysis with gProfileR
* [Report](#Report) - QBiC Report describing results of the whole pipeline

## DE analysis

Differential expression analysis is perfomed with the [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) R package.

**Output directory: `results/DESeq2`**

* `metadata/`
  * `metadata.tsv`: metadata sheet used by the pipeline.
  * `contrasts.tsv`: contrasts used for DE analysis (if provided).
  * `linear_model.txt`: linear model used for DE analysis.
  * `gene_list.txt`: provided list of interesting genes (if provided).
* `gene_counts_tables/`
  * `raw_gene_counts.txt`: raw gene counts table from the nf-core/rnaseq pipeline and used for the differential gene expression analysis.
  * `rlog_transformed_gene_counts.tsv`: normalized gene counts with the "regularized logarithm" approach. This normalization is used prior to PCA analysis and heatmap plotting of the gene counts.
  * `vst_transformed_gene_counts.tsv`: normalized gene counts with the "variance stabilizing" transformation.
  * `sizeFactor_libraries.tsv`: size factors for each sample.
* `DE_genes_tables/`: folder containing one tab-separated table for each of the contrasts in the analysis. Each table contains a list of all differentially expressed genes in the contrast, specifying the mean gene expression across all samples (baseMean), and the log2 fold change value and p-adjusted values (padj) for this contrast.
* `final_gene_table/final_gene_list_DESeq2.tsv`: table containing a list of all genes considered in the analysis. Here a summary of the log2 Fold Change and p-adjusted values for all contrasts is displayed. Additionally, the column **filter** shows if this gene was differentially expressed (DE) in any of the contrasts, or not (not_DE). The column **contrast_vector** contains for each contrast considered in the analysis a 1 if the gene was differentially expressed for this contrast or a 0 if it was not.
* `plots/`:
  * `Heatmaps_of_distances.pdf/.svg`: heatmap of the pairwise euclidean distances among samples, when the rlog normalized gene counts are considered.
  * `PCA_plot.pdf`: PCA of the rlog normalized gene counts.
  * `boxplots_example_genes/`: boxplots of the normalized gene counts for each of the sample groups for example genes.
  * `boxplots_requested_genes/`: boxplots of the normalized gene counts for each of the sample groups for the list of requested genes.
  * `further_diagnostics_plots/`: plots for diagnostics of the differential gene analysis procedure.

## Pathway analysis

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `trim_galore` directory.

**Output directory: `results/fastqc`**

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files
* `zips/sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images


## MultiQC
[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/multiqc`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see [http://multiqc.info](http://multiqc.info)
