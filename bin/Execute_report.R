#!/usr/bin/env Rscript

library(rmarkdown)
library(optparse)

option_list = list(
  make_option(c("-r", "--report"), type="character", default=NULL, help="report template file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="RNAseq_report.html", help="output file name", metavar="character"),
  make_option(c("-s", "--qc_summary"), type="character", default=NULL, help="qc summary file", metavar="character"),
  make_option(c("-v", "--versions"), type="character", default=NULL, help="versions file", metavar="character"),
  make_option(c("-m", "--model"), type="character", default=NULL, help="linear model file", metavar="character"),
  make_option(c("-c", "--config"), type="character", default=NULL, help="report config file", metavar="character"),
  make_option(c("-k", "--contrast"), type="character", default=NULL, help="contrasts file", metavar="character"),
  make_option(c("-f", "--fastqc"), type="character", default=NULL, help="fastqc zipped folder", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
#qc_summary_ = paste('"', opt$qc_summary, '"', sep='')


rmarkdown::render(opt$report, output_file = opt$output, knit_root_dir = '.', 
                  params = list(path_qc_summary = paste('\"', opt$qc_summary, '\"', sep=''),
                                path_versions = opt$versions,
                                path_design = opt$model,
                                path_config = opt$config,
                                path_contrast = opt$contrast,
                                path_fastqc = opt$fastqc))