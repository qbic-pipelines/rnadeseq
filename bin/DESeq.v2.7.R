#!/usr/bin/env Rscript

# Analysis pipeline from raw read count table to differential expression lists using DESeq2
# Author: Stefan Czemmel
# Contributors: Gisela Gabernet
# MIT License

library("RColorBrewer")
library("reshape2")
library("genefilter")
library("DESeq2")
library("ggplot2")
library("plyr")
library("vsn")
library("gplots")
library("pheatmap")
library("optparse")
library("svglite")
library("extrafont")

#clean up graphs
graphics.off()

theme_set(theme_classic())
extrafont::font_import()
extrafont::loadfonts()

######################################
# create directories needed
ifelse(!dir.exists("DESeq2"), dir.create("DESeq2"), FALSE)
dir.create("DESeq2/raw_counts")
dir.create("DESeq2/results")
dir.create("DESeq2/metadata")
dir.create("DESeq2/results/plots")
dir.create("DESeq2/results/plots/plots_example_genes")
dir.create("DESeq2/results/plots/plots_requested_genes")
dir.create("DESeq2/results/count_tables")
dir.create("DESeq2/results/DE_genes_tables")
dir.create("DESeq2/results/final")
######################################

# #1)check input data path

# provide these files as arguments:
option_list = list(
  make_option(c("-c", "--counts"), type="character", default=NULL, help="raw count table path", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL, help="metadata table path", metavar="character"),
  make_option(c("-d", "--design"), type="character", default=NULL, help="design linear model path", metavar="character"),
  make_option(c("-k", "--contrasts"), type="character", default=NULL, help="contrast matrix file", metavar="character"),
  make_option(c("-l", "--genelist"), type="character", default=NULL, help="gene list file", metavar="character"),
  make_option(c("-t", "--logFCthreshold"), type="integer", default=0, help="Log 2 Fold Change threshold for DE genes", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

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
if(!is.null(opt$contrasts)){
  path_contrasts = opt$contrasts
}
if(!is.null(opt$genelist)){
  requested_genes_path = opt$genelist
}

##2) load count table

count.table <- read.table(path_count_table,  header = T,sep = "\t",na.strings =c("","NA"),quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE,row.names=1)
count.table$Ensembl_ID <- row.names(count.table)
drop <- c("Ensembl_ID","gene_name")
gene_names <- count.table[,drop]


##Need to reduce gene names to QBiC codes
names(count.table) <- gsub('([A-Z0-9]{10})\\Aligned\\.(.*)\\.(.*)','\\1', names(count.table))
write.table(count.table , file = "DESeq2/raw_counts/raw_counts.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,  col.names = T, qmethod = c("escape", "double"))

count.table <- count.table[ , !(names(count.table) %in% drop)]

###2.1) remove lines with "__" from HTSeq, not needed for featureCounts (will not harm here)
count.table <- count.table[!grepl("^__",row.names(count.table)),]
#do some hard filtering just based on number of 0's per row
count.table = count.table[rowSums(count.table)>0,]

###3) load metadata: sample preparations tsv file from qPortal

system(paste("cp ",metadata_path," DESeq2/metadata/",sep=""))
m <- read.table(metadata_path, sep="\t", header=TRUE,na.strings =c("","NaN"),quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE,row.names = 1)


#make sure m is factor where needed
names(m) = gsub("Condition..","condition_",names(m))
conditions = names(m)[grepl("condition_",names(m))]
sapply(m,class)
for (i in conditions) {
  m[,i] = as.factor(m[,i])
}
sapply(m,class)


## Need to order columns in count.table
count.table <- count.table[, order(names(count.table))]

#check m and count table sorting, (correspond to QBiC codes): if not in the same order stop
stopifnot(identical(names(count.table),row.names(m)))

#change row.names now:
m$Secondary.Name <- gsub(" ; ", "_", m$Secondary.Name)
m$nn = paste(row.names(m),m$Secondary.Name,sep="_")
names(count.table) = m$nn
row.names(m) = m$nn

stopifnot(identical(names(count.table),row.names(m)))


#to get all possible pairwise comparisons, make a combined factor

conditions <- grepl(colnames(m),pattern = "condition_")
m$x <- apply(as.data.frame(m[ ,conditions]),1,paste, collapse = "_")
print(m)


###4) run DESeq function
design <- read.csv(path_design, sep="\t", header = F)
cds <- DESeqDataSetFromMatrix( countData =count.table, colData =m, design = eval(parse(text=as.character(design[[1]]))))
cds <- DESeq(cds,  parallel = FALSE)


### 4.1) sizeFactors(cds) as indicator of library sequencing depth
sizeFactors(cds)
write.table(sizeFactors(cds),paste("DESeq2/results/count_tables/sizeFactor_libraries.tsv",sep=""), append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,  col.names = F, qmethod = c("escape", "double"))

#write raw counts to file
write.table(count.table, paste("DESeq2/results/count_tables/raw.read.counts.tsv",sep=""), append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,  col.names = NA, qmethod = c("escape", "double"))

###4.2) contrasts
coefficients <- resultsNames(cds)
bg = data.frame(bg = character(nrow(cds)))
if (!is.null(opt$contrasts)){
  contrasts <- read.table(path_contrasts, sep="\t", header = T)
  stopifnot(length(coefficients)==nrow(contrasts))

  ## Contrast calculation
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
    write.table(d1DE, file=paste("DESeq2/results/DE_genes_tables/DE_contrast_",contname,".tsv",sep=""), sep="\t", quote=F, col.names = T, row.names = F)
    names(d1) = paste(names(d1),contname,sep="_")
    bg = cbind(bg,d1)
  }
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
    write.table(d1DE, file=paste("DESeq2/results/DE_genes_tables/DE_contrast_",contname,".tsv",sep=""), sep="\t", quote=F, col.names = T, row.names = F)
    names(d1) = paste(names(d1),contname,sep="_")
    bg = cbind(bg,d1)
  }
}


#remove identical columns
bg$bg <- NULL
idx <- duplicated(t(bg))
bg <- bg[, !idx]
bg$Ensembl_ID <- row.names(bg)
bg <- bg[,c(dim(bg)[2],1:dim(bg)[2]-1)]
names(bg)[1:2] = c("Ensembl_ID","baseMean")
names(bg)

#4.3) get DE genes from any contrast
padj=names(bg)[grepl("padj",names(bg))]
logFC = names(bg)[grepl("log2FoldChange", names(bg))]
logFC = bg[,logFC]
padj = bg[,padj]
padj[is.na(padj)] <- 1
padj_bin = data.matrix(ifelse(padj < 0.05, 1, 0))
logFC_bin = data.matrix(ifelse(abs(logFC) > opt$logFCthreshold, 1, 0))
DE_bin = padj_bin * logFC_bin
DE_bin = as.data.frame(DE_bin)
cols <- names(padj)
DE_bin$filter <- apply(DE_bin[ ,cols],1,paste, collapse = "-")
DE_bin$Ensembl_ID = row.names(padj)
DE_bin = DE_bin[,c("Ensembl_ID","filter")]

#make final data frame
bg1 = merge(bg,DE_bin,by.x="Ensembl_ID",by.y="Ensembl_ID")
stopifnot(identical(dim(bg1)[1],dim(assay(cds))[1]))
bg1$filter1 = ifelse(grepl("1",bg1$filter),"DE","not_DE")
bg1 = merge(x=bg1, y=gene_names, by.x="Ensembl_ID", by.y="Ensembl_ID", all.x = T)
bg1 = bg1[,c(dim(bg1)[2],1:dim(bg1)[2]-1)]
bg1 = bg1[order(bg1[,"Ensembl_ID"]),]

#4.4) extract ID for genes to plot, make 20 plots:
kip <- subset(bg1, filter1 == "DE")
kip = unique(kip$Ensembl_ID)

if (length(kip) > 20) {
  set.seed(10)
  kip1 = sample(kip,size = 2)
} else {
  kip1 = kip
  }

for (i in kip1){
  d <- plotCounts(cds, gene=i, intgroup=c("x"), returnData=TRUE,normalized = T)
  d$variable = row.names(d)
  plot <- ggplot(data=d, aes(x=x, y=count, fill=x)) +
            geom_boxplot(position=position_dodge()) +
            geom_jitter(position=position_dodge(.8)) +
            ggtitle(paste("Gene ",i,sep="")) + xlab("") + ylab("Normalized gene counts") + theme_bw() +
            theme(text = element_text(size=12),
               axis.text.x = element_text(angle=45, vjust=1,hjust=1))
  ggsave(filename=paste("DESeq2/results/plots/plots_example_genes/",i,".svg",sep=""), width=10, height=5, plot=plot)
  print(i)
}

if (!is.null(opt$genelist)){
  #4.5) make plots of interesting genes
  gene_ids <- read.table(requested_genes_path, col.names = "requested_gene_name")
  gene_ids$requested_gene_name <- sapply(gene_ids$requested_gene_name, toupper)
  bg1$gene_name <- sapply(bg1$gene_name, toupper)

  kip2 <- subset(bg1, gene_name %in% gene_ids$requested_gene_name)
  kip2_Ensembl <- kip2$Ensembl_ID
  kip2_gene_name <- kip2$gene_name
  for (i in c(1:length(kip2_Ensembl)))
  {
    d <- plotCounts(cds, gene=kip2_Ensembl[i], intgroup=c("x"), returnData=TRUE,normalized = T)
    d$variable = row.names(d)
    plot <- ggplot(data=d, aes(x=x, y=count, fill=x)) +
      geom_boxplot(position=position_dodge()) +
      geom_jitter(position=position_dodge(.8)) +
      ggtitle(paste("Gene ",kip2_gene_name[i],sep="")) + xlab("") + ylab("Normalized gene counts") + theme_bw() +
      theme(text = element_text(size=12),
            axis.text.x = element_text(angle=45, vjust=1,hjust=1))
    ggsave(filename=paste("DESeq2/results/plots/plots_requested_genes/",kip2_gene_name[i],"_",kip2_Ensembl[i],".svg",sep=""), plot=plot)
    print(kip2_gene_name[i])
  }
}

#write to file
write.table(bg1, "DESeq2/results/final/final_list_DESeq2.tsv", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,  col.names = T, qmethod = c("escape", "double"))

##5) Data transformation
#rlog
rld <- rlog(cds)
##vst
vsd <- varianceStabilizingTransformation(cds)

#write normalized values to a file
write.table(assay(rld), "DESeq2/results/count_tables/rlog_transformed.read.counts.tsv", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = NA, qmethod = c("escape", "double"))
write.table(assay(vsd), "DESeq2/results/count_tables/vst_transformed.read.counts.tsv", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = NA, qmethod = c("escape", "double"))

##6) Diagnostic plots

#Cooks distances: get important for example when checking knock-out and overexpression studies
pdf("DESeq2/results/plots/Cooks-distances.pdf")
par(mar=c(10,3,3,3))
par( mfrow = c(1,2))
boxplot(log10(assays(cds)[["cooks"]]), range=0, las=2,ylim = c(-15, 15),main="log10-Cooks")
boxplot(log2(assays(cds)[["cooks"]]), range=0, las=2,ylim = c(-15, 15),main="log2-Cooks")
dev.off()

#The function plotDispEsts visualizes DESeqs dispersion estimates: 
pdf("DESeq2/results/plots/Dispersion_plot.pdf")
plotDispEsts(cds, ylim = c(1e-5, 1e8))
dev.off()

#Effects of transformations on the variance
notAllZero <- (rowSums(counts(cds))>0) 
pdf("DESeq2/results/plots/Effects_of_transformations_on_the_variance.pdf")
par(oma=c(3,3,3,3))
par(mfrow = c(1, 3))
meanSdPlot(log2(counts(cds,normalized=TRUE)[notAllZero,] + 1),ylab  = "sd raw count data")
meanSdPlot(assay(rld[notAllZero,]),ylab  = "sd rlog transformed count data")
meanSdPlot(assay(vsd[notAllZero,]),ylab  = "sd vst transformed count data")
dev.off()

###Sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
colours = colorRampPalette(rev(brewer.pal(9, "Blues")))(255) 

### Visualization of distance using Heatmaps
pdf("DESeq2/results/plots/Heatmaps_of_distances.pdf")
par(oma=c(3,3,3,3))
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colours,fontsize=6)
dev.off()

svg("DESeq2/results/plots/Heatmaps_of_distances.svg")
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colours,fontsize=6)
dev.off()

#### Visualization of distance using PCA plots
pcaData <- plotPCA(rld,intgroup=c("x"),ntop = dim(rld)[1], returnData=TRUE)
percentVar <- round(100*attr(pcaData, "percentVar"))
pca <- ggplot(pcaData, aes(PC1, PC2, color=x)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2], "% variance")) +
  coord_fixed()
ggsave(plot = pca, filename = "DESeq2/results/plots/PCA_plot.pdf", device = "pdf", dpi = 300)
ggsave(plot = pca, filename = "DESeq2/results/plots/PCA_plot.svg", device = "svg", dpi = 150)

###Gene clustering
topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE), 50)
pdf("DESeq2/results/plots/heatmap_of_top50_genes_with_most_variance_across_samples.pdf")
par(oma=c(3,3,3,3))
heatmap.2(assay(rld)[topVarGenes, ],scale="row",trace="none",dendrogram="col",col=colorRampPalette( rev(brewer.pal(9, "RdBu")))(255),cexRow=0.5,cexCol=0.5)
dev.off()

svg("DESeq2/results/plots/heatmap_of_top50_genes_with_most_variance_across_samples.svg")
par(oma=c(3,3,3,3))
heatmap.2(assay(rld)[topVarGenes, ],scale="row",trace="none",dendrogram="col",col=colorRampPalette( rev(brewer.pal(9, "RdBu")))(255),cexRow=0.5,cexCol=0.5)
dev.off()

#further diagnostics plots
dir.create("DESeq2/results/plots/further_diagnostics")
res=0
for (i in resultsNames(cds)[-1]) {
  res = results(cds,name = i)
  pdf(paste("DESeq2/results/plots/further_diagnostics/all_results_MA_plot_",i,".pdf",sep=""))
  plotMA(res,ylim = c(-4, 4))
  dev.off()
  #multiple hyptothesis testing
  qs <- c( 0, quantile(results(cds)$baseMean[res$baseMean > 0], 0:4/4 ))
  bins <- cut(res$baseMean, qs )
  # rename the levels of the bins using the middle point
  levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
  # 3) calculate the ratio of p values less than .01 for each bin
  ratios <- tapply(res$pvalue, bins, function(p) mean(p < .01, na.rm=TRUE ))
  # 4) plot these ratios
  pdf(paste("DESeq2/results/plots/further_diagnostics/dependency_small.pval_mean_normal.counts_",i,".pdf",sep=""))
  barplot(ratios, xlab="mean normalized count", ylab="ratio of small p values")
  dev.off()
  #plot number of rejections
  pdf(paste("DESeq2/results/plots/further_diagnostics/number.of.rejections_",i,".pdf",sep=""))
  plot(metadata(res)$filterNumRej,
       type="b", ylab="number of rejections",
       xlab="quantiles of filter")
  lines(metadata(res)$lo.fit, col="red")
  abline(v=metadata(res)$filterTheta)
  dev.off()
  #Histogram of passed and rejected hypothesis
  use <- res$baseMean > metadata(res)$filterThreshold
  table(use)
  h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
  h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
  colori <- c('do not pass'="khaki", 'pass'="powderblue")
  pdf(paste("DESeq2/results/plots/further_diagnostics/histogram_of_p.values",i,".pdf",sep=""))
  barplot(height = rbind(h1$density, h2$density), beside = FALSE,
          col = colori, space = 0, main = "", xlab="p value",ylab="frequency")
  text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
       adj = c(0.5,1.7), xpd=NA)
  legend("topleft", fill=rev(colori), legend=rev(names(colori)))
  dev.off()
  rm(res,qs,bins,ratios,use,h1,h2,colori)
}

#end of script
####-------------save Sessioninfo
fn <- paste("DESeq2/sessionInfo_",format(Sys.Date(), "%d_%m_%Y"),".txt",sep="")
sink(fn)
sessionInfo()
sink()
####---------END----save Sessioninfo
