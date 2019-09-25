# qbicsoftware/rnadeseq: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

<!-- TODO qbicsoftware: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [DE analysis](#DE-analysis) - differential expression analysis with DESeq2
* [Pathway analysis](#Pathway-analysis) - pathway analysis with gProfileR
* [Report](#Report) - QBiC Report describing results of the whole pipeline

## DE analysis

Differential expression analysis is perfomed with the [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) R package.

**Output directory: `results/DESeq2`**

This directory contains the zipped results. When unzipping them, the following subfolders can be found.

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

Pathway analysis with [gProfileR](https://biit.cs.ut.ee/gprofiler/gost) R package.

**Output directory: `results/pathway_analysis/`**

This directory contains the zipped pathway analysis results (`gProfileR.zip`). When unzipping them, a subfolder for each contrast used for the differential gene expression analysis is found. Inside each contrast folder, there is the following output:

* `*_keg_pathway_enrichment_plot.pdf/png`
  * Barplots showing the proportion of differentially expressed genes in the pathway.
* `KEGG_pathways/`
  * Contains the KEGG pathways graphs with the log fold change of the differentially expressed genes.
* `pathway_heatmaps`
  * Contains heatmaps of the normalized gene counts for each of the differentially expressed pathways.

## Report

QBiC report for RNAseq analysis.

**Output directory: `results/report/`**

In this directory the zipped report is contained. This file is ready for upload to the project in the QBiC portal. Once unzipped, this directory contains:

* `Report.html`
  * QBiC report describing the RNAseq results.
* `YYYYMMDD_PiName_QXXXX_signed.pdf`
  * Signed offer for the QXXXX project
* `DESeq2/`
  * DESeq2 results to be attached to the report (see DESeq2 output description).
* `QC/`
  * Quality control results. Contains:
  * `fastqc.zip`: zipped fastqc results for each sample.
  * `multiqc_data`: data from the multiQC process of the RNAseq pipeline which is incorporated into the report.
  * `multiqc_plots`: plots from the multiQC process of the RNAseq pipeline which are incorporated into the report.
  * `qc_summary.tsv`: summary of the multiQC process in a table format.
