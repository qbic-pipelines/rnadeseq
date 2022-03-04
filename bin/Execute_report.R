#!/usr/bin/env Rscript

library(rmarkdown)
library(optparse)

option_list = list(
    make_option(c("-r", "--report"), type="character", default=NULL, help="report template file", metavar="character"),
    make_option(c("-o", "--output"), type="character", default="RNAseq_report.html", help="output file name", metavar="character"),
    make_option(c("-s", "--proj_summary"), type="character", default=NULL, help="project summary file", metavar="character"),
    make_option(c("-v", "--versions"), type="character", default=NULL, help="versions file", metavar="character"),
    make_option(c("-m", "--model"), type="character", default=NULL, help="linear model file", metavar="character"),
    make_option(c("-k", "--contrasts"), type="character", default=NULL, help="contrasts file", metavar="character"),
    make_option(c("-l", "--genelist"), type="character", default=NULL, help="path to gene list file", metavar="character"),
    make_option(c("-q", "--quote"), type="character", default=NULL, help="path to the signed quote PDF file", metavar="character"),
    make_option(c("-g", "--organism"), type="character", default=NULL, help="Organism, e.g. Hsapiens."),
    make_option(c("-b", "--batch_effect"), action="store_true", default=FALSE, help="Batch effect correction."),
    make_option(c("-f", "--log_FC"), type="double", default=NULL, help="Log Fold Change threshold to consider a gene DE."),
    make_option(c("-x", "--revision"), type="character", default=NULL, help="rnadeseq workflow revision", metavar="character"),
    make_option(c("-p", "--min_DEG_pathway"), type="integer", default=NULL, help="min. number of genes DE in a pathway for this pathway to be considered enriched.", metavar="integer"),
    make_option(c("-a", "--pathway_analysis"), action="store_true", default=FALSE, help="Pathway analysis."),
    make_option(c("-y", "--rlog"), action="store_true", default=FALSE, help="Use rlog instead of vst normalization.")
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
                                params = list(path_proj_summary = opt$proj_summary,
                                path_versions = opt$versions,
                                path_design = opt$model,
                                path_wd = wd,
                                path_contrasts = opt$contrasts,
                                path_genelist = path_genelist,
                                path_min_DEG = opt$min_DEG_pathway,
                                path_quote = opt$quote,
                                organism = opt$organism,
                                batch_effect = opt$batch_effect,
                                log_FC = opt$log_FC,
                                revision = opt$revision,
                                pathway_analysis = opt$pathway_analysis,
                                rlog = opt$rlog))
