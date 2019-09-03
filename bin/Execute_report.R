#!/usr/bin/env Rscript

library(rmarkdown)
library(optparse)

option_list = list(
  make_option(c("-r", "--report"), type="character", default=NULL, help="report template file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="RNAseq_report.html", help="output file name", metavar="character"),
  make_option(c("-s", "--proj_summary"), type="character", default=NULL, help="project summary file", metavar="character"),
  make_option(c("-v", "--versions"), type="character", default=NULL, help="versions file", metavar="character"),
  make_option(c("-m", "--model"), type="character", default=NULL, help="linear model file", metavar="character"),
  make_option(c("-c", "--config"), type="character", default=NULL, help="report config file", metavar="character"),
  make_option(c("-k", "--contrasts"), type="character", default=NULL, help="contrasts file", metavar="character"),
  make_option(c("-l", "--genelist"), type="character", default=NULL, help="path to gene list file", metavar="character"),
  make_option(c("-q", "--quote"), type="character", default=NULL, help="path to the signed quote PDF file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

wd=getwd()
if(!is.null(opt$contrasts)){
  path_contrasts = opt$contrasts
} else {
  path_contrasts = ''
}
if(!is.null(opt$genelist)){
  path_genelist = opt$genelist
} else {
  path_genelist = ''
}

rmarkdown::render(opt$report, output_file = opt$output, knit_root_dir = wd, output_dir = wd,
                  params = list(path_proj_summary = opt$proj_summary,
                                path_versions = opt$versions,
                                path_design = opt$model,
                                path_config = opt$config,
                                path_wd = wd,
                                path_contrasts = path_contrasts,
                                path_genelist = path_genelist
                                path_quote = opt$quote))