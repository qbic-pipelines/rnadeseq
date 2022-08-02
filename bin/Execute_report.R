#!/usr/bin/env Rscript

library(rmarkdown)
library(optparse)

option_list = list(
    make_option(c("-r", "--report"), type="character", default=NULL, help="report template file", metavar="character"),
    make_option(c("-o", "--output"), type="character", default="RNAseq_report.html", help="output file name", metavar="character"),

    make_option(c("-y", "--input_type"), type="character", default="featurecounts", help="Which type of input data is provided; must be one of [featurecounts, rsem, salmon]", metavar="character"),
    make_option(c("-c", "--gene_counts"), type="character", default=NULL, help="path to raw count table or directory", metavar="character"),
    make_option(c("-m", "--metadata"), type="character", default=NULL, help="path to metadata table", metavar="character"),
    make_option(c("-d", "--model"), type="character", default=NULL, help="path to linear model design file", metavar="character"),
    make_option(c("-t", "--gtf"), type="character", default=NULL, help="path to gtf table if using salmon/rsem input", metavar="character"),

    make_option(c("-x", "--contrast_matrix"), type="character", default=NULL, help="path to contrast matrix file", metavar="character"),
    make_option(c("-k", "--contrast_list"), type="character", default=NULL, help="path to contrast list file", metavar="character"),
    make_option(c("-p", "--contrast_pairs"), type="character", default=NULL, help="path to contrast pairs file", metavar="character"),
    make_option(c("-l", "--genelist"), type="character", default=NULL, help="path to gene list file", metavar="character"),
    make_option(c("-e", "--relevel"), type="character", default=NULL, help="path to factor relevel file", metavar="character"),
    make_option(c("-b", "--batch_effect"), action="store_true", default=FALSE, help="Batch effect correction."),
    make_option(c("-f", "--log_FC_threshold"), type="double", default=NULL, help="Log Fold Change threshold to consider a gene DE."),

    make_option(c("-n", "--nsub_genes"), type="integer", default=NULL, help="subset number of genes for vst."),
    make_option(c("-z", "--rlog"), action="store_true", default=TRUE, help="Use rlog instead of vst normalization."),

    make_option(c("-a", "--pathway_analysis"), action="store_true", default=FALSE, help="Pathway analysis."),
    make_option(c("-g", "--organism"), type="character", default=NULL, help="Organism, e.g. Hsapiens."),
    make_option(c("-i", "--species_library"), type="character", default=NULL, help="Library name. Example format: org.At.tair.db", metavar="character"),
    make_option(c("-u", "--keytype"), type="character", default=NULL, help="Keytype. Example format: TAIR (varies greatly depending on library!)", metavar="character"),
    make_option(c("-w", "--min_DEG_pathway"), type="integer", default=NULL, help="min. number of genes DE in a pathway for this pathway to be considered enriched.", metavar="integer"),

    make_option(c("-s", "--proj_summary"), type="character", default=NULL, help="project summary file", metavar="character"),
    make_option(c("-v", "--versions"), type="character", default=NULL, help="versions file", metavar="character"),
    make_option(c("-j", "--revision"), type="character", default=NULL, help="rnadeseq workflow revision", metavar="character"),

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
                                log_FC_threshold = opt$log_FC_threshold,
                                nsub_genes = opt$nsub_genes,
                                rlog = opt$rlog,

                                pathway_analysis = opt$pathway_analysis,
                                organism = opt$organism,
                                species_library = opt$species_library,
                                keytype = opt$keytype,
                                min_DEG_pathway = opt$min_DEG_pathway,

                                path_proj_summary = opt$proj_summary,
                                path_versions = opt$versions,
                                revision = opt$revision,

                                citest = opt$citest))
