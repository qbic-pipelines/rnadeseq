#!/usr/bin/env Rscript

# Differential expression analysis from raw read count table using DESeq2
# Author: Stefan Czemmel
# Contributors: Gisela Gabernet
# QBiC 2019; MIT License

library(RColorBrewer)
library(reshape2)
library(genefilter)
library(DESeq2)
library(ggplot2)
library(plyr)
library(vsn)
library(gplots)
library(pheatmap)
library(optparse)
library(svglite)
library(extrafont)
library(limma)

# clean up graphs
graphics.off()

# Load fonts and set plot themes
theme_set(theme_bw(base_family = "ArialMT") + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(family="ArialMT")))
extrafont::font_import()
extrafont::loadfonts()

# create directories needed
ifelse(!dir.exists("differential_gene_expression"), dir.create("differential_gene_expression"), FALSE)
dir.create("differential_gene_expression/metadata")
dir.create("differential_gene_expression/plots")
dir.create("differential_gene_expression/plots/boxplots_example_genes")
dir.create("differential_gene_expression/plots/boxplots_requested_genes")
dir.create("differential_gene_expression/plots/further_diagnostics_plots")
dir.create("differential_gene_expression/gene_counts_tables")
dir.create("differential_gene_expression/DE_genes_tables")
dir.create("differential_gene_expression/final_gene_table")

# check input data path
# provide these files as arguments:
option_list = list(
  make_option(c("-c", "--counts"), type="character", default=NULL, help="path to raw count table", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL, help="path to metadata table", metavar="character"),
  make_option(c("-d", "--design"), type="character", default=NULL, help="path to linear model design file", metavar="character"),
  make_option(c("-x", "--contrasts_matrix"), type="character", default=NULL, help="path to contrasts matrix file", metavar="character"),
  make_option(c("-r", "--relevel"), type="character", default=NULL, help="path to factor relevel file", metavar="character"),
  make_option(c("-k", "--contrasts_list"), type="character", default=NULL, help="path to contrasts list file", metavar="character"),
  make_option(c("-p", "--contrasts_pairs"), type="character", default=NULL, help="path to contrasts pairs file", metavar="character"),
  make_option(c("-l", "--genelist"), type="character", default=NULL, help="path to gene list file", metavar="character"),
  make_option(c("-t", "--logFCthreshold"), type="integer", default=0, help="Log 2 Fold Change threshold for DE genes", metavar="character"),
  make_option(c("-b", "--batchEffect"), default=FALSE, action="store_true", help="Whether to consider batch effects in the DESeq2 analysis", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Validate and read input
if (is.null(opt$counts)){
  print_help(opt_parser)
  stop("Counts table needs to be provided!")
} else {
  path_count_table = opt$counts
}
if (is.null(opt$metadata)){
  print_help(opt_parser)
  stop("Metadata table needs to be provided!")
} else {
  metadata_path = opt$metadata
}
if (is.null(opt$design)){
  print_help(opt_parser)
  stop("Linear model design file needs to be provided!")
} else {
  path_design = opt$design
}
if(!is.null(opt$relevel)){
  path_relevel = opt$relevel
}
if(!is.null(opt$contrasts_matrix)){
  if(!is.null(opt$contrasts_list) & !is.null(opt$contrasts_pairs)) {
    stop("Provide only one of contrasts_matrix / contrasts_list / contrasts pairs!")
  }
  path_contrasts_matrix = opt$contrasts_matrix
}
if(!is.null(opt$contrasts_list)){
  if(!is.null(opt$contrasts_matrix) & !is.null(opt$contrasts_pairs)) {
    stop("Provide only one of contrasts_matrix / contrasts_list / contrasts pairs!")
  }
  path_contrasts_list = opt$contrasts_list
}
if(!is.null(opt$contrasts_pairs)){
  if(!is.null(opt$contrasts_matrix) & !is.null(opt$contrasts_list)) {
    stop("Provide only one of contrasts_matrix / contrasts_list / contrasts pairs!")
  }
  path_contrasts_pairs = opt$contrasts_pairs
}
if(!is.null(opt$genelist)){
  requested_genes_path = opt$genelist
}

####### LOADING AND PROCESSING COUNT TABLE AND METADATA TABLE #####################################

# Load count table
count.table <- read.table(path_count_table,  header = T,sep = "\t",na.strings =c("","NA"),quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE,row.names=1)
count.table$Ensembl_ID <- row.names(count.table)
drop <- c("Ensembl_ID","gene_name")
gene_names <- count.table[,drop]

# Reduce sample names to QBiC codes in count table
names(count.table) <- substr(names(count.table), 1, 10)
count.table <- count.table[ , !(names(count.table) %in% drop)]

# Remove lines with "__" from HTSeq, not needed for featureCounts (will not harm here)
count.table <- count.table[!grepl("^__",row.names(count.table)),]
# Do some hard filtering for genes with 0 expression
count.table = count.table[rowSums(count.table)>0,]

# Load metadata: sample preparations tsv file from qPortal
metadata <- read.table(metadata_path, sep="\t", header=TRUE,na.strings =c("","NaN"), quote=NULL, stringsAsFactors=F, dec=".", fill=TRUE, row.names=1)
system(paste("mv ",metadata_path," differential_gene_expression/metadata/metadata.tsv",sep=""))

# Make sure metadata is factor where needed
names(metadata) = gsub("Condition..","condition_",names(metadata))
conditions = names(metadata)[grepl("condition_",names(metadata))]
for (i in conditions) {
  metadata[,i] = as.factor(metadata[,i])
}

# Need to order columns in count.table
count.table <- count.table[, order(names(count.table))]

print(names(count.table))
print(row.names(metadata))

# check metadata and count table sorting, (correspond to QBiC codes): if not in the same order stop
stopifnot(identical(names(count.table),row.names(metadata)))

# process secondary names and change row names in metadata
metadata$Secondary.Name <- gsub(" ; ", "_", metadata$Secondary.Name)
metadata$Secondary.Name <- gsub(" ", "_", metadata$Secondary.Name)
metadata$sampleName = paste(row.names(metadata),metadata$Secondary.Name,sep="_")
names(count.table) = metadata$sampleName
row.names(metadata) = metadata$sampleName

stopifnot(identical(names(count.table),row.names(metadata)))

#to get all possible pairwise comparisons, make a combined factor

conditions <- grepl(colnames(metadata),pattern = "condition_")
metadata$combfactor <- apply(as.data.frame(metadata[ ,conditions]),1,paste, collapse = "_")

# Read design file
design <- read.csv(path_design, sep="\t", header = F)
write.table(design, file="differential_gene_expression/metadata/linear_model.txt", sep="\t", quote=F, col.names = F, row.names = F)

################## RUN DESEQ2 ######################################

# Apply relevel if provided to metadata
if (!is.null(opt$relevel)) {
  relevel <- read.table(path_relevel, sep="\t", header = T, colClasses = "character")
  write.table(relevel, file="differential_gene_expression/metadata/relevel.tsv")

  for (i in c(1:nrow(relevel))) {
    relev <- relevel[i,]
    metadata[,relev[1]] <- relevel(metadata[,relev[1]], relev[2])
  }
}

# Run DESeq function
cds <- DESeqDataSetFromMatrix( countData =count.table, colData =metadata, design = eval(parse(text=as.character(design[[1]]))))
cds <- DESeq(cds,  parallel = FALSE)

# SizeFactors(cds) as indicator of library sequencing depth
sizeFactors(cds)
write.table(sizeFactors(cds),paste("differential_gene_expression/gene_counts_tables/sizeFactor_libraries.tsv",sep=""), append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,  col.names = F, qmethod = c("escape", "double"))

# Write raw counts to file
count_table_names <- merge(x=gene_names, y=count.table, by.x = "Ensembl_ID", by.y="row.names")
write.table(count_table_names, paste("differential_gene_expression/gene_counts_tables/raw_gene_counts.tsv",sep=""), append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F, qmethod = c("escape", "double"))

# Contrasts
coefficients <- resultsNames(cds)
coef_tab <- data.frame(coef=coefficients)
write.table(coef_tab,file="differential_gene_expression/metadata/DESeq2_coefficients.tsv", sep="\t", quote=F, col.names = T, row.names = F)

bg = data.frame(bg = character(nrow(cds)))

if (!is.null(opt$contrasts_matrix)){
  contrasts <- read.table(path_contrasts_matrix, sep="\t", header = T)
  write.table(contrasts, file="differential_gene_expression/metadata/contrasts.tsv", sep="\t", quote=F, col.names = T, row.names = F)
  
  if(length(coefficients) != nrow(contrasts)){
    print(coefficients)
    stop("Error: Your contrast table has different number of rows than the number of coefficients in the DESeq2 model.")
  }

  ## Contrast calculation for contrast table
  for (i in c(1:ncol(contrasts))) {
    d1 <-results(cds, contrast=contrasts[[i]])

    contname <- names(contrasts[i])
    d1 <- as.data.frame(d1)
    print(contname)
    d1_name <- d1
    d1_name$Ensembl_ID = row.names(d1)
    d1_name <- merge(x=d1_name, y=gene_names, by.x = "Ensembl_ID", by.y="Ensembl_ID", all.x=T)
    d1_name = d1_name[,c(dim(d1_name)[2],1:dim(d1_name)[2]-1)]
    d1_name = d1_name[order(d1_name[,"Ensembl_ID"]),]
    d1DE <- subset(d1_name, padj < 0.05 & (log2FoldChange > opt$logFCthreshold | log2FoldChange < opt$logFCthreshold))
    d1DE <- d1DE[order(d1DE$padj),]
    write.table(d1DE, file=paste("differential_gene_expression/DE_genes_tables/DE_contrast_",contname,".tsv",sep=""), sep="\t", quote=F, col.names = T, row.names = F)
    names(d1) = paste(names(d1),contname,sep="_")
    bg = cbind(bg,d1)
  }
  write(colnames(contrasts), file="contrast_names.txt", sep="\t")

} else if (!is.null(opt$contrasts_list)) {
  contrasts <- read.table(path_contrasts_list, sep="\t", header=T, colClasses = "character")
  write.table(contrasts, file="differential_gene_expression/metadata/contrasts.tsv")

  # Contrast calculation for contrast list
  contnames <- c()
  for (i in c(1:nrow(contrasts))) {
    cont <- as.character(contrasts[i,])
    contname <- paste0(cont[1], "_", cont[2], "_vs_", cont[3])
# TODO: add checks if provided contnames and factors are in metadata
    d1 <- results(cds, contrast=cont)
    d1 <- as.data.frame(d1)
    print(contname)
    d1_name <- d1
    d1_name$Ensembl_ID = row.names(d1)
    d1_name <- merge(x=d1_name, y=gene_names, by.x = "Ensembl_ID", by.y="Ensembl_ID", all.x=T)
    d1_name = d1_name[,c(dim(d1_name)[2],1:dim(d1_name)[2]-1)]
    d1_name = d1_name[order(d1_name[,"Ensembl_ID"]),]
    d1DE <- subset(d1_name, padj < 0.05 & (log2FoldChange > opt$logFCthreshold | log2FoldChange < opt$logFCthreshold))
    d1DE <- d1DE[order(d1DE$padj),]
    write.table(d1DE, file=paste("differential_gene_expression/DE_genes_tables/DE_contrast_",contname,".tsv",sep=""), sep="\t", quote=F, col.names = T, row.names = F)
    names(d1) = paste(names(d1),contname,sep="_")
    bg = cbind(bg,d1)

    contnames <- append(contnames, contname)
  }
  write(contnames, file="contrast_names.txt", sep="\t")

} else if (!is.null(opt$contrasts_pairs)) {
    contrasts <- read.table(path_contrasts_pairs, sep="\t", header = T, colClasses = "character")
    write.table(contrasts, file="differential_gene_expression/metadata/contrasts.tsv", sep="\t", quote=F, col.names = T, row.names = F)

    # Contrast calculation for contrast pairs
    contnames <- c()
    for (i in c(1:nrow(contrasts))) {
      cont <- as.character(contrasts[i,])
      contname <- cont[0]
      if (!(cont[2] %in% coefficients & cont[3] %in% coefficients)){
        stop(paste0("Provided contrast name is invalid, it needs to be contained in ", coefficients))
      } 
      d1 <- results(cds, contrast=list(cont[1],cont[2]))
      d1 <- as.data.frame(d1)
      print(contname)
      d1_name <- d1
      d1_name$Ensembl_ID = row.names(d1)
      d1_name <- merge(x=d1_name, y=gene_names, by.x = "Ensembl_ID", by.y="Ensembl_ID", all.x=T)
      d1_name = d1_name[,c(dim(d1_name)[2],1:dim(d1_name)[2]-1)]
      d1_name = d1_name[order(d1_name[,"Ensembl_ID"]),]
      d1DE <- subset(d1_name, padj < 0.05 & (log2FoldChange > opt$logFCthreshold | log2FoldChange < opt$logFCthreshold))
      d1DE <- d1DE[order(d1DE$padj),]
      write.table(d1DE, file=paste("differential_gene_expression/DE_genes_tables/DE_contrast_",contname,".tsv",sep=""), sep="\t", quote=F, col.names = T, row.names = F)
      names(d1) = paste(names(d1),contname,sep="_")
      bg = cbind(bg,d1)

      contnames <- append(contnames, contname)
    }
    write(contnames, file="contrast_names.txt", sep="\t")

} else {
  for (contname in coefficients[2:length(coefficients)]) {
    d1 <- results(cds, name=contname)
    d1 <- as.data.frame(d1)
    print(contname)
    d1_name <- d1
    d1_name$Ensembl_ID = row.names(d1)
    d1_name <- merge(x=d1_name, y=gene_names, by.x ="Ensembl_ID", by.y="Ensembl_ID", all.x=T)
    d1_name = d1_name[,c(dim(d1_name)[2],1:dim(d1_name)[2]-1)]
    d1_name = d1_name[order(d1_name[,"Ensembl_ID"]),]
    d1DE <- subset(d1_name, padj < 0.05 & (log2FoldChange > opt$logFCthreshold | log2FoldChange < opt$logFCthreshold))
    write.table(d1DE, file=paste("differential_gene_expression/DE_genes_tables/DE_contrast_",contname,".tsv",sep=""), sep="\t", quote=F, col.names = T, row.names = F)
    names(d1) = paste(names(d1),contname,sep="_")
    bg = cbind(bg,d1)
  }
  write(coefficients[2:length(coefficients)], file="contrast_names.txt", sep="\t")
}

# Remove identical columns
bg$bg <- NULL
idx <- duplicated(t(bg))
bg <- bg[, !idx]
bg$Ensembl_ID <- row.names(bg)
bg <- bg[,c(dim(bg)[2],1:dim(bg)[2]-1)]
names(bg)[1:2] = c("Ensembl_ID","baseMean")

# Get DE genes from any contrast
padj=names(bg)[grepl("padj",names(bg))]
logFC = names(bg)[grepl("log2FoldChange", names(bg))]
logFC = bg[,logFC,drop=F]
padj = bg[,padj,drop=F]
padj[is.na(padj)] <- 1
padj_bin = data.matrix(ifelse(padj < 0.05, 1, 0))
logFC_bin = data.matrix(ifelse(abs(logFC) > opt$logFCthreshold, 1, 0))
DE_bin = padj_bin * logFC_bin
DE_bin = as.data.frame(DE_bin)
cols <- names(padj)

if (ncol(DE_bin)>1){
  DE_bin$contrast_vector <- apply(DE_bin[ ,cols],1,paste, collapse = "-")
  DE_bin$Ensembl_ID = row.names(padj)
} else {
  DE_bin$contrast_vector <- DE_bin[,1]
  DE_bin$Ensembl_ID = row.names(padj)
}

DE_bin = DE_bin[,c("Ensembl_ID","contrast_vector")]

# Make final data frame
bg1 = merge(bg,DE_bin,by.x="Ensembl_ID",by.y="Ensembl_ID")
stopifnot(identical(dim(bg1)[1],dim(assay(cds))[1]))
bg1$outcome = ifelse(grepl("1",bg1$contrast_vector),"DE","not_DE")
bg1 = merge(x=bg1, y=gene_names, by.x="Ensembl_ID", by.y="Ensembl_ID", all.x = T)
bg1 = bg1[,c(dim(bg1)[2],1:dim(bg1)[2]-1)]
bg1 = bg1[order(bg1[,"Ensembl_ID"]),]

#write to file
write.table(bg1, "differential_gene_expression/final_gene_table/final_DE_gene_list.tsv", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,  col.names = T, qmethod = c("escape", "double"))

############### BOXPLOTS GENE EXPRESSION PER CONDITION ##########################
# extract ID for genes to plot, make 20 plots:
kip <- subset(bg1, outcome == "DE")
kip = unique(kip$Ensembl_ID)

if (length(kip) > 20) {
  set.seed(10)
  kip1 = sample(kip,size = 2)
} else {
  kip1 = kip
  }

for (i in kip1){
  d <- plotCounts(cds, gene=i, intgroup=c("combfactor"), returnData=TRUE, normalized = T)
  d$variable = row.names(d)
  plot <- ggplot(data=d, aes(x=combfactor, y=count, fill=combfactor)) +
            geom_boxplot(position=position_dodge()) +
            geom_jitter(position=position_dodge(.8)) +
            ggtitle(paste("Gene ",i,sep="")) + xlab("") + ylab("Normalized gene counts") + theme_bw() +
            theme(text = element_text(size=12),
               axis.text.x = element_text(angle=45, vjust=1,hjust=1))
  ggsave(filename=paste("differential_gene_expression/plots/boxplots_example_genes/",i,".svg",sep=""), width=10, height=5, plot=plot)
  ggsave(filename=paste("differential_gene_expression/plots/boxplots_example_genes/",i,".png",sep=""), width=10, height=5, plot=plot)

  print(i)
}

# make plots of interesting genes in gene list
if (!is.null(opt$genelist)){
  gene_ids <- read.table(requested_genes_path, col.names = "requested_gene_name")
  write.table(gene_ids, file="differential_gene_expression/metadata/gene_list.txt", col.names=F, row.names=F, sep="\t")
  gene_ids$requested_gene_name <- sapply(gene_ids$requested_gene_name, toupper)
  bg1$gene_name <- sapply(bg1$gene_name, toupper)

  kip2 <- subset(bg1, gene_name %in% gene_ids$requested_gene_name)
  kip2_Ensembl <- kip2$Ensembl_ID
  kip2_gene_name <- kip2$gene_name
  for (i in c(1:length(kip2_Ensembl)))
  {
    d <- plotCounts(cds, gene=kip2_Ensembl[i], intgroup=c("combfactor"), returnData=TRUE,normalized = T)
    d$variable = row.names(d)
    plot <- ggplot(data=d, aes(x=combfactor, y=count, fill=combfactor)) +
      geom_boxplot(position=position_dodge()) +
      geom_jitter(position=position_dodge(.8)) +
      ggtitle(paste("Gene ",kip2_gene_name[i],sep="")) + xlab("") + ylab("Normalized gene counts") + theme_bw() +
      theme(text = element_text(size=12),
            axis.text.x = element_text(angle=45, vjust=1,hjust=1))
    ggsave(filename=paste("differential_gene_expression/plots/boxplots_requested_genes/",kip2_gene_name[i],"_",kip2_Ensembl[i],".svg",sep=""), plot=plot)
    ggsave(filename=paste("differential_gene_expression/plots/boxplots_requested_genes/",kip2_gene_name[i],"_",kip2_Ensembl[i],".png",sep=""), plot=plot)
    print(kip2_gene_name[i])
  }
}

############################## TRANSFORMED AND NORMALIZED COUNTS ###################
# rlog transformation
rld <- rlog(cds, blind=FALSE)
# vst transformation
vsd <- vst(cds, blind=FALSE)

# write normalized values to a file
rld_names <- merge(x=gene_names, y=assay(rld), by.x = "Ensembl_ID", by.y="row.names")
write.table(rld_names, "differential_gene_expression/gene_counts_tables/rlog_transformed_gene_counts.tsv", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,  qmethod = c("escape", "double"))
vsd_names <- merge(x=gene_names, y=assay(vsd), by.x = "Ensembl_ID", by.y="row.names")
write.table(vsd_names, "differential_gene_expression/gene_counts_tables/vst_transformed_gene_counts.tsv", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F, qmethod = c("escape", "double"))


##################  SAMPLE DISTANCES HEATMAP ##################
# Sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colours = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pdf("differential_gene_expression/plots/Heatmaps_of_distances.pdf")
par(oma=c(3,3,3,3))
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colours,fontsize=11)
dev.off()

svg("differential_gene_expression/plots/Heatmaps_of_distances.svg")
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colours,fontsize=11)
dev.off()


############################ PCA PLOTS ########################
pcaData <- plotPCA(vsd,intgroup=c("combfactor"),ntop = dim(vsd)[1], returnData=TRUE)
percentVar <- round(100*attr(pcaData, "percentVar"))
pca <- ggplot(pcaData, aes(PC1, PC2, color=combfactor)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2], "% variance")) +
  theme(legend.title = element_blank()) +
  coord_fixed()
ggsave(plot = pca, filename = "differential_gene_expression/plots/PCA_plot.pdf", device = "pdf", dpi = 300)
ggsave(plot = pca, filename = "differential_gene_expression/plots/PCA_plot.svg", device = "svg", dpi = 150)

########################### PCA PLOT with batch-corrected data ############
if(opt$batchEffect){
  assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)
  pcaData2 <- plotPCA(vsd, intgroup=c("combfactor"), ntop = dim(vsd)[1], returnData=TRUE)
  percentVar <- round(100*attr(pcaData, "percentVar"))
  pca2 <- ggplot(pcaData2, aes(PC1, PC2, color=combfactor)) +
            geom_point(size=3)+
            xlab(paste0("PC1: ", percentVar[1],"% variance")) +
            ylab(paste0("PC2: ", percentVar[2], "% variance")) +
            theme(legend.title = element_blank()) +
            coord_fixed()
  ggsave(plot = pca2, filename = "differential_gene_expression/plots/PCA_batch_corrected_plot.pdf", device = "pdf", dpi=300)
  ggsave(plot = pca2, filename = "differential_gene_expression/plots/PCA_batch_corrected_plot.svg", device = "svg", dpi = 150)
}


############################## DIAGNOSTICS AND QUALITY CONTROL PLOTS ###############################

# Cooks distances: get important for example when checking knock-out and overexpression studies
pdf("differential_gene_expression/plots/further_diagnostics_plots/Cooks-distances.pdf")
par(mar=c(10,3,3,3))
par( mfrow = c(1,2))
boxplot(log10(assays(cds)[["cooks"]]), range=0, las=2,ylim = c(-15, 15),main="log10-Cooks")
boxplot(log2(assays(cds)[["cooks"]]), range=0, las=2,ylim = c(-15, 15),main="log2-Cooks")
dev.off()

# The function plotDispEsts visualizes DESeqs dispersion estimates: 
pdf("differential_gene_expression/plots/further_diagnostics_plots/Dispersion_plot.pdf")
plotDispEsts(cds, ylim = c(1e-5, 1e8))
dev.off()

# Effects of transformations on the variance
notAllZero <- (rowSums(counts(cds))>0) 
pdf("differential_gene_expression/plots/further_diagnostics_plots/Effects_of_transformations_on_the_variance.pdf")
par(oma=c(3,3,3,3))
par(mfrow = c(1, 3))
meanSdPlot(log2(counts(cds,normalized=TRUE)[notAllZero,] + 1),ylab  = "sd raw count data")
meanSdPlot(assay(rld[notAllZero,]),ylab  = "sd rlog transformed count data")
meanSdPlot(assay(vsd[notAllZero,]),ylab  = "sd vst transformed count data")
dev.off()

# Further diagnostics plots
res=0
for (i in resultsNames(cds)[-1]) {
  res = results(cds,name = i)
  pdf(paste("differential_gene_expression/plots/further_diagnostics_plots/all_results_MA_plot_",i,".pdf",sep=""))
  plotMA(res,ylim = c(-4, 4))
  dev.off()
  # multiple hyptothesis testing
  qs <- c( 0, quantile(results(cds)$baseMean[res$baseMean > 0], 0:4/4 ))
  bins <- cut(res$baseMean, qs )
  # rename the levels of the bins using the middle point
  levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
  # calculate the ratio of p values less than .01 for each bin
  ratios <- tapply(res$pvalue, bins, function(p) mean(p < .01, na.rm=TRUE ))
  # plot these ratios
  pdf(paste("differential_gene_expression/plots/further_diagnostics_plots/dependency_small.pval_mean_normal.counts_",i,".pdf",sep=""))
  barplot(ratios, xlab="mean normalized count", ylab="ratio of small p values")
  dev.off()
  # plot number of rejections
  pdf(paste("differential_gene_expression/plots/further_diagnostics_plots/number.of.rejections_",i,".pdf",sep=""))
  plot(metadata(res)$filterNumRej,
       type="b", ylab="number of rejections",
       xlab="quantiles of filter")
  lines(metadata(res)$lo.fit, col="red")
  abline(v=metadata(res)$filterTheta)
  dev.off()
  # Histogram of passed and rejected hypothesis
  use <- res$baseMean > metadata(res)$filterThreshold
  table(use)
  h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
  h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
  colori <- c('do not pass'="khaki", 'pass'="powderblue")
  pdf(paste("differential_gene_expression/plots/further_diagnostics_plots/histogram_of_p.values",i,".pdf",sep=""))
  barplot(height = rbind(h1$density, h2$density), beside = FALSE,
          col = colori, space = 0, main = "", xlab="p value",ylab="frequency")
  text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
       adj = c(0.5,1.7), xpd=NA)
  legend("topleft", fill=rev(colori), legend=rev(names(colori)))
  dev.off()
  rm(res,qs,bins,ratios,use,h1,h2,colori)
}
