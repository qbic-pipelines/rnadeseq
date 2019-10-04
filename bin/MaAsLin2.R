#!/usr/bin/env Rscript

# Multivariate Association with Linear Models (MaAsLin2) for normalized data (e.g. relative abundance)
# Author: Daniel Straub
# Template by: Stefan Czemmel
# QBiC 2019; MIT License

library("Maaslin2")
library("optparse")

# #1)check input data path

# provide these files as arguments:
option_list = list(
  make_option(c("-d", "--data"), type="character", default=NULL, help="path to normalized table", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL, help="path to metadata table", metavar="character"),
  make_option(c("-c", "--cores"), type="numeric", default=1, help="The number of R processes to run in parallel", metavar="character"),
  make_option(c("-o", "--output_folder"), type="character", default=1, help="the output folder, will be created", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$data)){
  print_help(opt_parser)
  stop("Normalized data table needs to be provided!")
} else {
  path_normalized_table = opt$data
}
if (is.null(opt$metadata)){
  print_help(opt_parser)
  stop("Metadata table needs to be provided!")
} else {
  metadata_path = opt$metadata
}
number_of_cores = opt$cores
output_folder = opt$output_folder

##2) load normalized table

norm.table <- read.table(path_normalized_table, header = T ,sep = "\t", na.strings =c("","NA"), quote=NULL, stringsAsFactors=F, dec=".", fill=TRUE, comment.char = "")#, row.names=1)
#put first column as row names
row.names(norm.table) = norm.table[, 1]
norm.table[, 1] = NULL
#rename columns to match metadata
names(norm.table) = gsub("_Abundance","",names(norm.table))
names(norm.table) = gsub(".RPKs","",names(norm.table))

print(head(norm.table, 5))

###3) load metadata: sample preparations tsv file from qPortal

m <- read.table(metadata_path, sep="\t", header=TRUE,na.strings =c("","NaN"),quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE,row.names = 1)
#system(paste("mv ",metadata_path," metadata.tsv",sep=""))

#Are all sample names from normalized table in metadata?
print(names(norm.table))
print(names(m))
stopifnot(names(norm.table) %in% row.names(m))

#make sure m is factor where needed
names(m) = gsub("Condition..","condition_",names(m))
conditions = names(m)[grepl("condition_",names(m))]
#sapply(m,class)
for (i in conditions) {
  m[,i] = as.factor(m[,i])
}
#sapply(m,class)

#extract only conditions
m_conditions = m[, conditions]

#extract samples names in normalized data table from metadata
m_conditions = m_conditions[row.names(m_conditions) %in% names(norm.table),]

#extract only conditions that have more than 1 level
m_conditions = Filter(function(x)(length(unique(x))>1), m_conditions)

print(head(m_conditions))

###4) Perform tests
fit_data <- Maaslin2(
    norm.table, 
    m_conditions,
    output_folder, 
    standardize = TRUE,
    plot_scatter = FALSE,
    cores = number_of_cores)

#write to file
#write.table(result, "DESeq2/final_gene_table/final_gene_list_DESeq2.tsv", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,  col.names = T, qmethod = c("escape", "double"))

#end of script
####-------------save Sessioninfo
fn <- paste(output_folder,"/sessionInfo_",format(Sys.Date(), "%d_%m_%Y"),".txt",sep="")
sink(fn)
sessionInfo()
sink()
####---------END----save Sessioninfo
