#!/usr/bin/env Rscript

library(rmarkdown)
library(optparse)

option_list = list(
    make_option(c("-r", "--report"), type="character", default=NULL, help="Report template file", metavar="character"),
    make_option(c("-o", "--output"), type="character", default="RNAseq_report.html", help="Output file name", metavar="character"),

    make_option(c("-y", "--input_type"), type="character", default="featurecounts", help="Which type of input data is provided; must be one of [featurecounts, rsem, salmon]", metavar="character"),
    make_option(c("-c", "--gene_counts"), type="character", default=NULL, help="Path to raw count table or directory", metavar="character"),
    make_option(c("-m", "--metadata"), type="character", default=NULL, help="Path to metadata table", metavar="character"),
    make_option(c("-d", "--model"), type="character", default=NULL, help="Path to linear model design file", metavar="character"),
    make_option(c("-t", "--gtf"), type="character", default=NULL, help="Path to gtf table if using salmon/rsem input", metavar="character"),

    make_option(c("-x", "--contrast_matrix"), type="character", default=NULL, help="Path to contrast matrix file", metavar="character"),
    make_option(c("-k", "--contrast_list"), type="character", default=NULL, help="Path to contrast list file", metavar="character"),
    make_option(c("-p", "--contrast_pairs"), type="character", default=NULL, help="Path to contrast pairs file", metavar="character"),
    make_option(c("-l", "--genelist"), type="character", default=NULL, help="Path to gene list file", metavar="character"),
    make_option(c("-e", "--relevel"), type="character", default=NULL, help="Path to factor relevel file", metavar="character"),
    make_option(c("-b", "--batch_effect"), action="store_true", default=FALSE, help="Batch effect correction."),
    make_option(c("-f", "--logFC_threshold"), type="double", default=NULL, help="Log Fold Change threshold to consider a gene DE."),
    make_option("--adj_pval_threshold", type="double", default=0.05, help="adjusted p value threshold to consider a gene DE."),
    make_option("--round_DE", type="integer", default=-1, help="How many decimals to keep after rounding the DE analysis values; if -1, will not round."),

    make_option(c("-z", "--norm_method"), type="character", default=NULL, help="Which transformation(s) to use."),
    make_option(c("-n", "--nsub_genes"), type="integer", default=NULL, help="Subset number of genes for vst."),

    make_option(c("-a", "--pathway_analysis"), action="store_true", default=FALSE, help="Whether to run pathway analysis."),
    make_option(c("-g", "--organism"), type="character", default=NULL, help="Organism, e.g. Hsapiens."),
    make_option(c("-i", "--species_library"), type="character", default=NULL, help="Library name. Example format: org.At.tair.db", metavar="character"),
    make_option(c("-u", "--keytype"), type="character", default=NULL, help="Keytype. Example format: TAIR (varies greatly depending on library!)", metavar="character"),
    make_option("--custom_gmt", type="character", default=NULL, help="Path to custom GMT file to query during pathway analysis.", metavar="character"),
    make_option("--set_background", action="store_true", default=TRUE, help="Whether to use a background list for pathway analysis; if true, will only consider expressed genes (i.e. mean counts > 0) for PA."),
    make_option("--custom_background", type="character", default=NULL, help="Path to a custom background list TXT for pathway analysis; if provided, will only consider genes in that list for PA."),
    make_option(c("-w", "--min_DEG_pathway"), type="integer", default=NULL, help="min. number of genes DE in a pathway for this pathway to be considered enriched.", metavar="integer"),
    make_option("--datasources", type="character", default=NULL, help="Which datasources to use for pathway analysis.", metavar="character"),
    make_option("--heatmaps_cluster_rows", action="store_true", default=FALSE, help="Whether to activate row clustering when generating heatmaps of gene expression in enriched pathways."),
    make_option("--heatmaps_cluster_cols", action="store_true", default=FALSE, help="Whether to activate column clustering when generating heatmaps of gene expression in enriched pathways."),
    make_option("--pathway_pval_threshold", type="double", default=-1, help="Which p value threshold to use for pathway analysis."),

    make_option(c("-s", "--proj_summary"), type="character", default=NULL, help="Project summary file", metavar="character"),
    make_option(c("--path_quote"), type="character", default=NULL, help="Path to the quote PDF", metavar="character"),
    make_option(c("-v", "--software_versions"), type="character", default=NULL, help="Versions file", metavar="character"),
    make_option(c("-j", "--revision"), type="character", default=NULL, help="Rnadeseq workflow revision", metavar="character"),
    make_option("--logo", type="character", default=NULL, help="Logo image file", metavar="character"),

    make_option("--citest", action="store_true", default=FALSE, help="Run github test and don't save pathway heatmaps.")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

wd=getwd()

if(!is.null(opt$genelist)){
    path_genelist = opt$genelist
} else {
    path_genelist = ''
}

rmarkdown::render(opt$report, output_file = opt$output, knit_root_dir = wd, output_dir = wd,
                    params = list(
                                input_type = opt$input_type,
                                path_gene_counts = opt$gene_counts,
                                path_metadata = opt$metadata,
                                path_design = opt$model,
                                path_gtf = opt$gtf,

                                path_contrast_matrix = opt$contrast_matrix,
                                path_contrast_list = opt$contrast_list,
                                path_contrast_pairs = opt$contrast_pairs,
                                path_genelist = path_genelist,
                                path_relevel = opt$relevel,
                                batch_effect = opt$batch_effect,
                                logFC_threshold = opt$logFC_threshold,
                                adj_pval_threshold = opt$adj_pval_threshold,
                                norm_method = opt$norm_method,
                                nsub_genes = opt$nsub_genes,
                                round_DE = opt$round_DE,

                                pathway_analysis = opt$pathway_analysis,
                                organism = opt$organism,
                                species_library = opt$species_library,
                                custom_gmt = opt$custom_gmt,
                                set_background = opt$set_background,
                                custom_background = opt$custom_background,
                                keytype = opt$keytype,
                                min_DEG_pathway = opt$min_DEG_pathway,
                                datasources = opt$datasources,
                                heatmaps_cluster_rows = opt$heatmaps_cluster_rows,
                                heatmaps_cluster_cols = opt$heatmaps_cluster_cols,
                                pathway_pval_threshold = opt$pathway_pval_threshold,

                                path_proj_summary = opt$proj_summary,
                                path_quote = opt$path_quote,
                                path_software_versions = opt$software_versions,
                                revision = opt$revision,
                                logo = opt$logo,

                                citest = opt$citest))
