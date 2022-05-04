#!/usr/bin/env Rscript

# Differential expression analysis from raw read count table using DESeq2
# Author: Gisela Gabernet, Stefan Czemmel
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
library(dplyr)

library(tximeta)
library(tximport)
library(SummarizedExperiment)
library(impute)

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
    make_option(c("-y", "--input_type"), type="character", default="rawcounts", help="Which type of input data is provided; must be one of [rawcounts, rsem, salmon]", metavar="character"),
    make_option(c("-c", "--gene_counts"), type="character", default=NULL, help="path to raw count table", metavar="character"),
    make_option(c("-m", "--metadata"), type="character", default=NULL, help="path to metadata table", metavar="character"),
    make_option(c("-d", "--model"), type="character", default=NULL, help="path to linear model design file", metavar="character"),
    make_option(c("-x", "--contrasts_matrix"), type="character", default=NULL, help="path to contrasts matrix file", metavar="character"),
    make_option(c("-f", "--gtf"), type="character", default=NULL, help="path to gtf table if using salmon/rsem input", metavar="character"),
    make_option(c("-r", "--relevel"), type="character", default=NULL, help="path to factor relevel file", metavar="character"),
    make_option(c("-k", "--contrasts_list"), type="character", default=NULL, help="path to contrasts list file", metavar="character"),
    make_option(c("-p", "--contrasts_pairs"), type="character", default=NULL, help="path to contrasts pairs file", metavar="character"),
    make_option(c("-l", "--genelist"), type="character", default=NULL, help="path to gene list file", metavar="character"),
    make_option(c("-t", "--logFCthreshold"), type="integer", default=0, help="Log 2 Fold Change threshold for DE genes", metavar="character"),
    make_option(c("-g", "--rlog"), type="logical", default=TRUE, help="if TRUE, perform rlog transformation", metavar="character"),
    make_option(c("-b", "--batchEffect"), default=FALSE, action="store_true", help="Whether to consider batch effects in the DESeq2 analysis", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (!(opt$input_type %in% c("rawcounts", "rsem", "salmon"))){
    stop(paste0("Wrong input type ", opt$input_type, ", must be one of [rawcounts, rsem, salmon]!"))
}
if (opt$input_type %in% c("rsem", "salmon") && is.null(opt$gtf)){
    stop(paste0("For input type salmon, gtf file needs to be provided!"))
}

# Validate and read input
if (is.null(opt$gene_counts)){
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
    path_design = opt$model
}
if (!is.null(opt$relevel)){
    path_relevel = opt$relevel
}
if (!is.null(opt$contrasts_matrix)){
    if (!is.null(opt$contrasts_list) & !is.null(opt$contrasts_pairs)) {
        stop("Provide only one of contrasts_matrix / contrasts_list / contrasts pairs!")
    }
    path_contrasts_matrix = opt$contrasts_matrix
}
if (!is.null(opt$contrasts_list)){
    if (!is.null(opt$contrasts_matrix) & !is.null(opt$contrasts_pairs)) {
        stop("Provide only one of contrasts_matrix / contrasts_list / contrasts pairs!")
    }
    path_contrasts_list = opt$contrasts_list
}
if (!is.null(opt$contrasts_pairs)){
    if (!is.null(opt$contrasts_matrix) & !is.null(opt$contrasts_list)) {
        stop("Provide only one of contrasts_matrix / contrasts_list / contrasts pairs!")
    }
    path_contrasts_pairs = opt$contrasts_pairs
}
if (!is.null(opt$genelist)){
    requested_genes_path = opt$genelist
}



####### LOADING AND PROCESSING COUNT TABLE AND METADATA TABLE #####################################
# Load metadata: sample preparations tsv file from qPortal
metadata <- read.table(metadata_path, sep="\t", header=TRUE,na.strings =c("","NaN"), quote=NULL, stringsAsFactors=F, dec=".", fill=TRUE, row.names=1)
system(paste("mv ",metadata_path," differential_gene_expression/metadata/metadata.tsv",sep=""))
# TODO: Is this fine, or should I use another column? Sample names maybe?
qbicCodes <- rownames(metadata)
# Make sure metadata is factor where needed
names(metadata) = gsub("Condition..","condition_",names(metadata))
conditions = names(metadata)[grepl("condition_",names(metadata))]
for (i in conditions) {
    metadata[,i] = as.factor(metadata[,i])
}

# Load count table
if (opt$input_type == "rawcounts"){
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
    # Need to order columns in count.table
    count.table <- count.table[, order(names(count.table))]

    print("Count table column headers are:")
    print(names(count.table))
    print("Metadata table row names are:")
    print(row.names(metadata))

    print("If count table headers do not exactly match the metadata table row names the pipeline will stop.")

    # check metadata and count table sorting, (correspond to QBiC codes): if not in the same order stop
    stopifnot(identical(names(count.table),row.names(metadata)))
}

# process secondary names and change row names in metadata
metadata$Secondary.Name <- gsub(" ; ", "_", metadata$Secondary.Name)
metadata$Secondary.Name <- gsub(" ", "_", metadata$Secondary.Name)
metadata$sampleName = paste(row.names(metadata),metadata$Secondary.Name,sep="_")
row.names(metadata) = metadata$sampleName
if (opt$input_type == "rawcounts"){
    names(count.table) = metadata$sampleName
    stopifnot(identical(names(count.table),row.names(metadata)))
}

# to get all possible pairwise comparisons, make a combined factor
conditions <- grepl(colnames(metadata),pattern = "condition_")
metadata$combfactor <- apply(as.data.frame(metadata[ ,conditions]),1,paste, collapse = "_")

# Read design file
design <- read.csv(path_design, sep="\t", header = F)
write.table(design, file="differential_gene_expression/metadata/linear_model.txt", sep="\t", quote=F, col.names = F, row.names = F)

################## RUN DESEQ2 ######################################
#Apply relevel if provided to metadata
if (!is.null(opt$relevel)) {
    relevel_table <- read.table(path_relevel, sep="\t", header = T, colClasses = "character")
    write.table(relevel_table, file="differential_gene_expression/metadata/relevel.tsv")

    for (i in c(1:nrow(relevel_table))) {
        relev <- relevel_table[i,]
        cds[[paste(relev[1])]] <- relevel(cds[[paste(relev[1])]], paste(relev[2]))
    }
}

# Run DESeq function
if (opt$input_type == "rawcounts") {
    cds <- DESeqDataSetFromMatrix( countData =count.table, colData =metadata, design = eval(parse(text=as.character(design[[1]]))))
    cds <- DESeq(cds,  parallel = FALSE)
} else if (opt$input_type %in% c("rsem", "salmon")) {
    ## Create a dataframe which consists of both the gene id and the transcript name
    gtf <- rtracklayer::import(opt$gtf)
    gtf <- as.data.frame(gtf, header=T)
    # TODO: For some gtf files, transcript_id does not work!!
    tx2gene_gtf <- gtf[c("transcript_id", "gene_id")]
    tx2gene_gtf <- distinct(tx2gene_gtf)
    tx2gene_gtf[] <- lapply(tx2gene_gtf, function(x) gsub("\\.\\d+", "", x))
    colnames(tx2gene_gtf) <- c("transcript_id", "gene_id") #, "TXID"
    gene_names <- gtf[c("gene_id", "gene_name")]
    colnames(gene_names) <- c("Ensembl_ID", "gene_name")
    gene_names <- distinct(gene_names)
    rownames(gene_names) <- gene_names[,1]

    if (opt$input_type == "rsem") {
        files <- file.path(gsub("/$", "", path_count_table), paste0(qbicCodes, ".genes.results"))
        if (!(all(file.exists(files)))) {
            stop("DESeq2.R could not find all of the specified .genes.results files!")
        }
        #Extract condition columns and other info for tximeta
        condition_names <- unlist(strsplit(design[,1], split = " "))
        condition_names <- grep("condition", condition_names, value=T)
        sampleconditions <- data.frame(metadata[,condition_names])
        colnames(sampleconditions) <- condition_names
        coldata <- data.frame(files = files, names= qbicCodes, sampleconditions)
        coldata$combfactor <- metadata$combfactor
        rownames(coldata) <- NULL

        #Do tximeta, this is necessary to run DESeq on rsem
        se <- tximeta(coldata, type="rsem", txIn=FALSE, txOut=FALSE, skipMeta=TRUE)
        assays(se)$length[ assays(se)$length == 0] <- NA # set these as missing
        #Impute lengths for the 0-length values:
        length_imp <- impute.knn(assays(se)$length)
        assays(se)$length <- length_imp$data
        #dds from SummarizedExperiment <se>, then run DESeq
        cds <- DESeqDataSet(se, design = as.formula(eval(parse(text=as.character(design[[1]])))))
        cds <- DESeq(cds)
    } else if (opt$input_type == "salmon") {
        files <- file.path(gsub("/$", "", path_count_table), qbicCodes, "quant.sf")
        if (!(all(file.exists(files)))) {
            stop("DESeq2.R could not find all of the specified quant.sf files!")
        }
        ## Import all of the samples information and transform the identifiers
        txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene_gtf, ignoreTxVersion = T)
        # Run cds with tximport depending on whether rsem or salmon was used
                #Extract condition columns and other info for tximeta
        condition_names <- unlist(strsplit(design[,1], split = " "))
        condition_names <- grep("condition", condition_names, value=T)
        sampleconditions <- data.frame(metadata[,condition_names])
        colnames(sampleconditions) <- condition_names
        coldata <- data.frame(files = files, names= qbicCodes, sampleconditions)
        coldata$combfactor <- metadata$combfactor
        rownames(coldata) <- qbicCodes
        cds <- DESeqDataSetFromTximport(txi=txi.salmon, colData =coldata, design = eval(parse(text=as.character(design[[1]]))))
        cds <- DESeq(cds)
    }
} else {
    stop("Input type must be one of [rawcounts, rsem, salmon]!")
}

# SizeFactors(cds) as indicator of library sequencing depth
write.table(sizeFactors(cds),paste("differential_gene_expression/gene_counts_tables/sizeFactor_libraries.tsv",sep=""), append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,  col.names = F, qmethod = c("escape", "double"))

# Write cds assay table to file
write.table(counts(cds, normalized=T), paste("differential_gene_expression/gene_counts_tables/deseq2_table.tsv", sep=""), append=F, quote = F, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T, qmethod = c("escape", "double"))
if (opt$input_type == "rawcounts"){
    # Write raw counts to file
    count_table_names <- merge(x=gene_names, y=count.table, by.x = "Ensembl_ID", by.y="row.names")
    write.table(count_table_names, paste("differential_gene_expression/gene_counts_tables/raw_gene_counts.tsv",sep=""), append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F, qmethod = c("escape", "double"))
}
# Contrasts coefficient table write in metadata
coefficients <- resultsNames(cds)
coef_tab <- data.frame(coef=coefficients)
write.table(coef_tab,file="differential_gene_expression/metadata/DESeq2_coefficients.tsv", sep="\t", quote=F, col.names = T, row.names = F)

# Start variables to store DE genes for all contrasts
DE_genes_df = data.frame(DE_genes_df = character(nrow(cds)))
contrast_names <- c()

if (!is.null(opt$contrasts_matrix)){
    contrasts <- read.table(path_contrasts_matrix, sep="\t", header = T, row.names = 1)
    write.table(contrasts, file="differential_gene_expression/metadata/contrast_matrix.tsv", sep="\t", quote=F, col.names = T, row.names = F)

    # Check that contrast matrix is valid
    if (length(coefficients) != nrow(contrasts)){
        stop("Error: Your contrast table has different number of rows than the number of coefficients in the DESeq2 model.")
    }

    ## Contrast calculation for contrast matrix
    for (i in c(1:ncol(contrasts))) {
        results_DEseq_contrast <-results(cds, contrast=contrasts[[i]])

        contname <- names(contrasts[i])
        results_DEseq_contrast <- as.data.frame(results_DEseq_contrast)
        print("Analyzing contrast:")
        print(contname)
        # Add gene name in table
        DE_genes_contrast_genename <- results_DEseq_contrast
        DE_genes_contrast_genename$Ensembl_ID = row.names(results_DEseq_contrast)
        DE_genes_contrast_genename <- merge(x=DE_genes_contrast_genename, y=gene_names, by.x = "Ensembl_ID", by.y="Ensembl_ID", all.x=T)
        DE_genes_contrast_genename = DE_genes_contrast_genename[,c(dim(DE_genes_contrast_genename)[2],1:dim(DE_genes_contrast_genename)[2]-1)]
        DE_genes_contrast_genename = DE_genes_contrast_genename[order(DE_genes_contrast_genename[,"Ensembl_ID"]),]
        # Save all DE genes (even if not significant) to separate files for the volcano plot
        write.table(DE_genes_contrast_genename, file=paste("differential_gene_expression/DE_genes_tables/DE_contrast_allgenes_",contname,".tsv",sep=""), sep="\t", quote=F, col.names = T, row.names = F)
        # Select only significantly DE
        DE_genes_contrast <- subset(DE_genes_contrast_genename, padj < 0.05 & abs(log2FoldChange) > opt$logFCthreshold)
        DE_genes_contrast <- DE_genes_contrast[order(DE_genes_contrast$padj),]
        # Save table
        write.table(DE_genes_contrast, file=paste("differential_gene_expression/DE_genes_tables/DE_contrast_",contname,".tsv",sep=""), sep="\t", quote=F, col.names = T, row.names = F)
        names(results_DEseq_contrast) = paste(names(results_DEseq_contrast),contname,sep="_")
        # Append to DE genes table for all contrasts
        DE_genes_df = cbind(DE_genes_df,results_DEseq_contrast)
    }
    contrast_names <- append(contrast_names, colnames(contrasts))
}

if (!is.null(opt$contrasts_list)) {
    contrasts <- read.table(path_contrasts_list, sep="\t", header=T, colClasses = "character")
    write.table(contrasts, file="differential_gene_expression/metadata/contrast_list.tsv", quote=F)

    ## Contrast calculation for contrast list
    for (i in c(1:nrow(contrasts))) {
        cont <- as.character(contrasts[i,])
        contname <- paste0(cont[1], "_", cont[2], "_vs_", cont[3])
    # TODO: add checks if provided contrast_names and factors are in metadata
        results_DEseq_contrast <- results(cds, contrast=cont)
        results_DEseq_contrast <- as.data.frame(results_DEseq_contrast)
        print(contname)
        # Add gene name in table
        DE_genes_contrast_genename <- results_DEseq_contrast
        DE_genes_contrast_genename$Ensembl_ID = row.names(results_DEseq_contrast)
        DE_genes_contrast_genename <- merge(x=DE_genes_contrast_genename, y=gene_names, by.x = "Ensembl_ID", by.y="Ensembl_ID", all.x=T)
        DE_genes_contrast_genename = DE_genes_contrast_genename[,c(dim(DE_genes_contrast_genename)[2],1:dim(DE_genes_contrast_genename)[2]-1)]
        DE_genes_contrast_genename = DE_genes_contrast_genename[order(DE_genes_contrast_genename[,"Ensembl_ID"]),]
        # Save all DE genes (even if not significant) to separate files for the volcano plot
        write.table(DE_genes_contrast_genename, file=paste("differential_gene_expression/DE_genes_tables/DE_contrast_allgenes_",contname,".tsv",sep=""), sep="\t", quote=F, col.names = T, row.names = F)
        # Select only significantly DE
        DE_genes_contrast <- subset(DE_genes_contrast_genename, padj < 0.05 & abs(log2FoldChange) > opt$logFCthreshold)
        DE_genes_contrast <- DE_genes_contrast[order(DE_genes_contrast$padj),]
        # Save table
        write.table(DE_genes_contrast, file=paste("differential_gene_expression/DE_genes_tables/DE_contrast_",contname,".tsv",sep=""), sep="\t", quote=F, col.names = T, row.names = F)
        names(results_DEseq_contrast) = paste(names(results_DEseq_contrast),contname,sep="_")
        # Append to DE genes table for all contrasts
        DE_genes_df = cbind(DE_genes_df,results_DEseq_contrast)

        contrast_names <- append(contrast_names, contname)
    }
}

if (!is.null(opt$contrasts_pairs)) {
    contrasts <- read.table(path_contrasts_pairs, sep="\t", header = T, colClasses = "character")
    write.table(contrasts, file="differential_gene_expression/metadata/contrast_pairs.tsv", sep="\t", quote=F, col.names = T, row.names = F)

    # Contrast calculation for contrast pairs
    for (i in c(1:nrow(contrasts))) {
        cont <- as.character(contrasts[i,])
        contname <- cont[0]
        if (!(cont[2] %in% coefficients & cont[3] %in% coefficients)){
            stop(paste0("Provided contrast name is invalid, it needs to be contained in ", coefficients))
        }
        results_DEseq_contrast <- results(cds, contrast=list(cont[1],cont[2]))
        results_DEseq_contrast <- as.data.frame(results_DEseq_contrast)
        print("Analyzing contrast:")
        print(contname)
        # Add gene name in table
        DE_genes_contrast_genename <- results_DEseq_contrast
        DE_genes_contrast_genename$Ensembl_ID = row.names(results_DEseq_contrast)
        DE_genes_contrast_genename <- merge(x=DE_genes_contrast_genename, y=gene_names, by.x = "Ensembl_ID", by.y="Ensembl_ID", all.x=T)
        DE_genes_contrast_genename = DE_genes_contrast_genename[,c(dim(DE_genes_contrast_genename)[2],1:dim(DE_genes_contrast_genename)[2]-1)]
        DE_genes_contrast_genename = DE_genes_contrast_genename[order(DE_genes_contrast_genename[,"Ensembl_ID"]),]
        # Save all DE genes (even if not significant) to separate files for the volcano plot
        write.table(DE_genes_contrast_genename, file=paste("differential_gene_expression/DE_genes_tables/DE_contrast_allgenes_",contname,".tsv",sep=""), sep="\t", quote=F, col.names = T, row.names = F)
        # Select only significantly DE
        DE_genes_contrast <- subset(DE_genes_contrast_genename, padj < 0.05 & abs(log2FoldChange) > opt$logFCthreshold)
        DE_genes_contrast <- DE_genes_contrast[order(DE_genes_contrast$padj),]
        # Save table
        write.table(DE_genes_contrast, file=paste("differential_gene_expression/DE_genes_tables/DE_contrast_",contname,".tsv",sep=""), sep="\t", quote=F, col.names = T, row.names = F)
        names(results_DEseq_contrast) = paste(names(results_DEseq_contrast),contname,sep="_")
        # Append to DE genes table for all contrasts
        DE_genes_df = cbind(DE_genes_df,results_DEseq_contrast)

        contrast_names <- append(contrast_names, contname)
    }
}

# Calculating DE genes for default contrasts (no contrast matrix or list or pairs provided)
if (is.null(opt$contrasts_matrix) & is.null(opt$contrasts_list) & is.null(opt$contrasts_pairs)) {
    contrast_names <- coefficients[2:length(coefficients)]
    for (contname in contrast_names) {
        results_DEseq_contrast <- results(cds, name=contname)
        results_DEseq_contrast <- as.data.frame(results_DEseq_contrast)
        print("Analyzing contrast:")
        print(contname)

        # Adding gene name to table
        DE_genes_contrast_genename <- results_DEseq_contrast
        DE_genes_contrast_genename$Ensembl_ID = row.names(results_DEseq_contrast)
        DE_genes_contrast_genename <- merge(x=DE_genes_contrast_genename, y=gene_names, by.x ="Ensembl_ID", by.y="Ensembl_ID", all.x=T)
        DE_genes_contrast_genename = DE_genes_contrast_genename[,c(dim(DE_genes_contrast_genename)[2],1:dim(DE_genes_contrast_genename)[2]-1)]
        DE_genes_contrast_genename = DE_genes_contrast_genename[order(DE_genes_contrast_genename[,"Ensembl_ID"]),]
        # Save all DE genes (even if not significant) to separate files for the volcano plot
        write.table(DE_genes_contrast_genename, file=paste("differential_gene_expression/DE_genes_tables/DE_contrast_allgenes_",contname,".tsv",sep=""), sep="\t", quote=F, col.names = T, row.names = F)
        # Select only significantly DE
        DE_genes_contrast <- subset(DE_genes_contrast_genename, padj < 0.05 & abs(log2FoldChange) > opt$logFCthreshold)
        DE_genes_contrast <- DE_genes_contrast[order(DE_genes_contrast$padj),]
        # Save table
        write.table(DE_genes_contrast, file=paste("differential_gene_expression/DE_genes_tables/DE_contrast_",contname,".tsv",sep=""), sep="\t", quote=F, col.names = T, row.names = F)
        names(results_DEseq_contrast) = paste(names(results_DEseq_contrast),contname,sep="_")
        # Append to DE genes table for all contrasts
        DE_genes_df = cbind(DE_genes_df,results_DEseq_contrast)
    }
}

# Write contrast names to file
write(contrast_names, file="contrast_names.txt", sep="\t")

# Remove identical columns of DE_genes_df
DE_genes_df$DE_genes_df <- NULL
idx <- duplicated(t(DE_genes_df))
DE_genes_df <- DE_genes_df[, !idx]
DE_genes_df$Ensembl_ID <- row.names(DE_genes_df)
DE_genes_df <- DE_genes_df[,c(dim(DE_genes_df)[2],1:dim(DE_genes_df)[2]-1)]
names(DE_genes_df)[1:2] = c("Ensembl_ID","baseMean")

# Get DE genes from any contrast
padj_cols=names(DE_genes_df)[grepl("padj",names(DE_genes_df))]
logFC_cols = names(DE_genes_df)[grepl("log2FoldChange", names(DE_genes_df))]
logFC = DE_genes_df[,logFC_cols,drop=F]
padj = DE_genes_df[,padj_cols,drop=F]
padj[is.na(padj)] <- 1
# Convert to binary (1/0) matrix if padj < 0.05 or not, respectively
padj_bin = data.matrix(ifelse(padj < 0.05, 1, 0))
# Convert to binary (1/0) matrix if logFC is bigger or smaller than threshold or not, respectively
logFC_bin = data.matrix(ifelse(abs(logFC) > opt$logFCthreshold, 1, 0))
# Multiply the two bin matrices --> if padj matrix value or LogFC matrix value is 0, will be 0
DE_bin = padj_bin * logFC_bin

# Save as data frame
DE_bin = as.data.frame(DE_bin)
cols <- names(padj)

# Contrast vector column -> contains 1 or 0 if gene was DE for each contrast
if (ncol(DE_bin)>1){
    DE_bin$contrast_vector <- apply(DE_bin[ ,cols],1,paste, collapse = "-")
    DE_bin$Ensembl_ID = row.names(padj)
} else {
    DE_bin$contrast_vector <- DE_bin[,1]
    DE_bin$Ensembl_ID = row.names(padj)
}

DE_bin = DE_bin[,c("Ensembl_ID","contrast_vector")]

# Add contrast vector to final DE genes data frame
DE_genes_final_table = merge(DE_genes_df,DE_bin,by.x="Ensembl_ID",by.y="Ensembl_ID")
stopifnot(identical(dim(DE_genes_final_table)[1],dim(assay(cds))[1]))

# Calculate outcome --> if gene is DE in any contrast, annotate as DE
DE_genes_final_table$outcome = ifelse(grepl("1",DE_genes_final_table$contrast_vector),"DE","not_DE")
DE_genes_final_table = merge(x=DE_genes_final_table, y=gene_names, by.x="Ensembl_ID", by.y="Ensembl_ID", all.x = T)

DE_genes_final_table = DE_genes_final_table[,c(dim(DE_genes_final_table)[2],1:dim(DE_genes_final_table)[2]-1)]
DE_genes_final_table = DE_genes_final_table[order(DE_genes_final_table[,"Ensembl_ID"]),]

#write to file
write.table(DE_genes_final_table, "differential_gene_expression/final_gene_table/final_DE_gene_list.tsv", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,  col.names = T, qmethod = c("escape", "double"))

############################## TRANSFORMED AND NORMALIZED COUNTS ###################
if (opt$rlog){
    # rlog transformation
    rld <- rlog(cds, blind=FALSE)
    rld_names <- merge(x=gene_names, y=assay(rld), by.x = "Ensembl_ID", by.y="row.names")
    write.table(rld_names, "differential_gene_expression/gene_counts_tables/rlog_transformed_gene_counts.tsv", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,  qmethod = c("escape", "double"))
} else {
    # vst transformation
    vsd <- vst(cds, blind=FALSE)
    vsd_names <- merge(x=gene_names, y=assay(vsd), by.x = "Ensembl_ID", by.y="row.names")
    write.table(vsd_names, "differential_gene_expression/gene_counts_tables/vst_transformed_gene_counts.tsv", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F, qmethod = c("escape", "double"))
}



############### BOXPLOTS GENE EXPRESSION PER CONDITION ##########################
# extract ID for genes to plot, make 20 plots:
DE_genes_plot <- subset(DE_genes_final_table, outcome == "DE")
DE_genes_plot = unique(DE_genes_plot$Ensembl_ID)

if (length(DE_genes_plot) > 20) {
    set.seed(10)
    random_DE_genes_plot = sample(DE_genes_plot,size = 2)
} else {
    random_DE_genes_plot = DE_genes_plot
}

for (i in random_DE_genes_plot){
    boxplot_counts <- plotCounts(cds, gene=i, intgroup=c("combfactor"), returnData=TRUE, normalized = T)
    boxplot_counts$variable = row.names(boxplot_counts)
    plot <- ggplot(data=boxplot_counts, aes(x=combfactor, y=count, fill=combfactor)) +
                geom_boxplot(position=position_dodge()) +
                geom_jitter(position=position_dodge(.8)) +
                ggtitle(paste("Gene ",i,sep="")) + xlab("") + ylab("Normalized gene counts") + theme_bw() +
                theme(text = element_text(size=12),
                axis.text.x = element_text(angle=45, vjust=1,hjust=1))
    ggsave(filename=paste("differential_gene_expression/plots/boxplots_example_genes/",i,".svg",sep=""), width=10, height=5, plot=plot)
    ggsave(filename=paste("differential_gene_expression/plots/boxplots_example_genes/",i,".png",sep=""), width=10, height=5, plot=plot)
    ggsave(filename=paste("differential_gene_expression/plots/boxplots_example_genes/",i,".pdf",sep=""), width=10, height=5, plot=plot)
}

# make boxplots of interesting genes in gene list
if (!is.null(opt$genelist)){
    gene_ids <- read.table(requested_genes_path, col.names = "requested_gene_name")
    write.table(gene_ids, file="differential_gene_expression/metadata/requested_gene_list.txt", col.names=F, row.names=F, sep="\t")
    gene_ids$requested_gene_name <- sapply(gene_ids$requested_gene_name, toupper)
    gene_names$gene_name <- sapply(gene_names$gene_name, toupper)
    # get Ensemble IDs from requested genes
    requested_genes_plot <- subset(gene_names, gene_name %in% gene_ids$requested_gene_name)

    # Check that genes are in the cds table
    requested_genes_plot <- subset(requested_genes_plot, requested_genes_plot$Ensembl_ID %in% row.names(cds))

    requested_genes_plot_Ensembl <- requested_genes_plot$Ensembl_ID
    requested_genes_plot_gene_name <- requested_genes_plot$gene_name


    for (i in seq_along(requested_genes_plot_Ensembl)) {
        boxplot_counts <- plotCounts(cds, gene=requested_genes_plot_Ensembl[i], intgroup=c("combfactor"), returnData=TRUE, normalized = T)
        boxplot_counts$variable = row.names(boxplot_counts)
        plot <- ggplot(data=boxplot_counts, aes(x=combfactor, y=count, fill=combfactor)) +
        geom_boxplot(position=position_dodge()) +
        geom_jitter(position=position_dodge(.8)) +
        ggtitle(paste("Gene ",requested_genes_plot_gene_name[i],sep="")) + xlab("") + ylab("Normalized gene counts") + theme_bw() +
        theme(text = element_text(size=12),
                axis.text.x = element_text(angle=45, vjust=1,hjust=1))
        ggsave(filename=paste("differential_gene_expression/plots/boxplots_requested_genes/",requested_genes_plot_gene_name[i],"_",requested_genes_plot_Ensembl[i],".svg",sep=""), width=10, height=5, plot=plot)
        ggsave(filename=paste("differential_gene_expression/plots/boxplots_requested_genes/",requested_genes_plot_gene_name[i],"_",requested_genes_plot_Ensembl[i],".png",sep=""), width=10, height=5, plot=plot)
        ggsave(filename=paste("differential_gene_expression/plots/boxplots_requested_genes/",requested_genes_plot_gene_name[i],"_",requested_genes_plot_Ensembl[i],".pdf",sep=""), width=10, height=5, plot=plot)
    }
}

##################  SAMPLE DISTANCES HEATMAP ##################
# Sample distances
sampleDists <- dist(t(assay(if (opt$rlog) rld else vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colours = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pdf("differential_gene_expression/plots/Heatmaps_of_distances.pdf")
par(oma=c(3,3,3,3))
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colours,fontsize=10)
dev.off()

svg("differential_gene_expression/plots/Heatmaps_of_distances.svg")
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colours,fontsize=10)
dev.off()

png("differential_gene_expression/plots/Heatmaps_of_distances.png")
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colours,fontsize=10)
dev.off()

############################ PCA PLOTS ########################
pcaData <- plotPCA(if (opt$rlog) rld else vsd,intgroup=c("combfactor"),ntop = dim(if (opt$rlog) rld else vsd)[1], returnData=TRUE)
percentVar <- round(100*attr(pcaData, "percentVar"))
pca <- ggplot(pcaData, aes(PC1, PC2, color=combfactor)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2], "% variance")) +
    theme(legend.title = element_blank()) +
    coord_fixed()
ggsave(plot = pca, filename = "differential_gene_expression/plots/PCA_plot.pdf", device = "pdf", dpi = 300)
ggsave(plot = pca, filename = "differential_gene_expression/plots/PCA_plot.svg", device = "svg", dpi = 150)
ggsave(plot = pca, filename = "differential_gene_expression/plots/PCA_plot.png", device = "png", dpi = 150)

########################### PCA PLOT with batch-corrected data ############
if(opt$batchEffect){
    assay(if (opt$rlog) rld else vsd) <- limma::removeBatchEffect(assay(if (opt$rlog) rld else vsd), (if (opt$rlog) rld else vsd)$batch)
    pcaData2 <- plotPCA(if (opt$rlog) rld else vsd, intgroup=c("combfactor"), ntop = dim(if (opt$rlog) rld else vsd)[1], returnData=TRUE)
    percentVar <- round(100*attr(pcaData, "percentVar"))
    pca2 <- ggplot(pcaData2, aes(PC1, PC2, color=combfactor)) +
                geom_point(size=3)+
                xlab(paste0("PC1: ", percentVar[1],"% variance")) +
                ylab(paste0("PC2: ", percentVar[2], "% variance")) +
                theme(legend.title = element_blank()) +
                coord_fixed()
    ggsave(plot = pca2, filename = "differential_gene_expression/plots/PCA_batch_corrected_plot.pdf", device = "pdf", dpi=300)
    ggsave(plot = pca2, filename = "differential_gene_expression/plots/PCA_batch_corrected_plot.svg", device = "svg", dpi = 150)
    ggsave(plot = pca2, filename = "differential_gene_expression/plots/PCA_batch_corrected_plot.png", device = "png", dpi = 150)
}


############################## DIAGNOSTICS AND QUALITY CONTROL PLOTS ###############################

# Cooks distances: get important for example when checking knock-out and overexpression studies
pdf("differential_gene_expression/plots/further_diagnostics_plots/Cooks-distances.pdf")
par(mar=c(10,3,3,3))
par( mfrow = c(1,2))
boxplot(log10(assays(cds)[["cooks"]]), range=0, las=2,ylim = c(-15, 15),main="log10-Cooks")
boxplot(log2(assays(cds)[["cooks"]]), range=0, las=2,ylim = c(-15, 15),main="log2-Cooks")
dev.off()

png("differential_gene_expression/plots/further_diagnostics_plots/Cooks-distances.png")
par(mar=c(10,3,3,3))
par( mfrow = c(1,2))
boxplot(log10(assays(cds)[["cooks"]]), range=0, las=2,ylim = c(-15, 15),main="log10-Cooks")
boxplot(log2(assays(cds)[["cooks"]]), range=0, las=2,ylim = c(-15, 15),main="log2-Cooks")
dev.off()

svg("differential_gene_expression/plots/further_diagnostics_plots/Cooks-distances.svg")
par(mar=c(10,3,3,3))
par( mfrow = c(1,2))
boxplot(log10(assays(cds)[["cooks"]]), range=0, las=2,ylim = c(-15, 15),main="log10-Cooks")
boxplot(log2(assays(cds)[["cooks"]]), range=0, las=2,ylim = c(-15, 15),main="log2-Cooks")
dev.off()

# The function plotDispEsts visualizes DESeqs dispersion estimates:
pdf("differential_gene_expression/plots/further_diagnostics_plots/Dispersion_plot.pdf")
plotDispEsts(cds, ylim = c(1e-5, 1e8))
dev.off()

png("differential_gene_expression/plots/further_diagnostics_plots/Dispersion_plot.png")
plotDispEsts(cds, ylim = c(1e-5, 1e8))
dev.off()

svg("differential_gene_expression/plots/further_diagnostics_plots/Dispersion_plot.svg")
plotDispEsts(cds, ylim = c(1e-5, 1e8))
dev.off()

# Effects of transformations on the variance
notAllZero <- (rowSums(counts(cds))>0)
pdf("differential_gene_expression/plots/further_diagnostics_plots/Effects_of_transformations_on_the_variance.pdf")
par(oma=c(3,3,3,3))
par(mfrow = c(1, 3))
meanSdPlot(log2(counts(cds,normalized=TRUE)[notAllZero,] + 1),ylab  = "sd raw count data")
meanSdPlot(assay((if (opt$rlog) rld else vsd)[notAllZero,]),ylab  = "sd rlog transformed count data")
meanSdPlot(assay((if (opt$rlog) rld else vsd)[notAllZero,]),ylab  = paste("sd ", if (opt$rlog) "rld" else "vsd" ," transformed count data"))
dev.off()

notAllZero <- (rowSums(counts(cds))>0)
png("differential_gene_expression/plots/further_diagnostics_plots/Effects_of_transformations_on_the_variance.png")
par(oma=c(3,3,3,3))
par(mfrow = c(1, 3))
#Should this be done for salmon and rsem as well?
meanSdPlot(log2(counts(cds,normalized=TRUE)[notAllZero,] + 1),ylab  = "sd raw count data")
meanSdPlot(assay((if (opt$rlog) rld else vsd)[notAllZero,]),ylab  = "sd rlog transformed count data")
meanSdPlot(assay((if (opt$rlog) rld else vsd)[notAllZero,]),ylab  = paste("sd ", if (opt$rlog) "rld" else "vsd" ," transformed count data"))
dev.off()

notAllZero <- (rowSums(counts(cds))>0)
svg("differential_gene_expression/plots/further_diagnostics_plots/Effects_of_transformations_on_the_variance.svg")
par(oma=c(3,3,3,3))
par(mfrow = c(1, 3))
meanSdPlot(log2(counts(cds,normalized=TRUE)[notAllZero,] + 1),ylab  = "sd raw count data")
meanSdPlot(assay((if (opt$rlog) rld else vsd)[notAllZero,]),ylab  = "sd rlog transformed count data")
meanSdPlot(assay((if (opt$rlog) rld else vsd)[notAllZero,]),ylab  = paste("sd ", if (opt$rlog) "rld" else "vsd" ," transformed count data"))
dev.off()

# Further diagnostics plots
res=0
for (i in resultsNames(cds)[-1]) {
    res = results(cds,name = i)
    pdf(paste("differential_gene_expression/plots/further_diagnostics_plots/all_results_MA_plot_",i,".pdf",sep=""))
    plotMA(res,ylim = c(-4, 4))
    dev.off()
    png(paste("differential_gene_expression/plots/further_diagnostics_plots/all_results_MA_plot_",i,".png",sep=""))
    plotMA(res,ylim = c(-4, 4))
    dev.off()
    svg(paste("differential_gene_expression/plots/further_diagnostics_plots/all_results_MA_plot_",i,".svg",sep=""))
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
    png(paste("differential_gene_expression/plots/further_diagnostics_plots/dependency_small.pval_mean_normal.counts_",i,".png",sep=""))
    barplot(ratios, xlab="mean normalized count", ylab="ratio of small p values")
    dev.off()
    svg(paste("differential_gene_expression/plots/further_diagnostics_plots/dependency_small.pval_mean_normal.counts_",i,".svg",sep=""))
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
    png(paste("differential_gene_expression/plots/further_diagnostics_plots/number.of.rejections_",i,".png",sep=""))
    plot(metadata(res)$filterNumRej,
        type="b", ylab="number of rejections",
        xlab="quantiles of filter")
    lines(metadata(res)$lo.fit, col="red")
    abline(v=metadata(res)$filterTheta)
    dev.off()
    svg(paste("differential_gene_expression/plots/further_diagnostics_plots/number.of.rejections_",i,".svg",sep=""))
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
    png(paste("differential_gene_expression/plots/further_diagnostics_plots/histogram_of_p.values",i,".png",sep=""))
    barplot(height = rbind(h1$density, h2$density), beside = FALSE,
            col = colori, space = 0, main = "", xlab="p value",ylab="frequency")
    text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
        adj = c(0.5,1.7), xpd=NA)
    legend("topleft", fill=rev(colori), legend=rev(names(colori)))
    dev.off()
    svg(paste("differential_gene_expression/plots/further_diagnostics_plots/histogram_of_p.values",i,".svg",sep=""))
    barplot(height = rbind(h1$density, h2$density), beside = FALSE,
            col = colori, space = 0, main = "", xlab="p value",ylab="frequency")
    text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
        adj = c(0.5,1.7), xpd=NA)
    legend("topleft", fill=rev(colori), legend=rev(names(colori)))
    dev.off()
    rm(res,qs,bins,ratios,use,h1,h2,colori)
}

