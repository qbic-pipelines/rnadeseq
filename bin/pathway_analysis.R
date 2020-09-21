#!/usr/bin/env Rscript

# Pathway analysis from differentially expressed gene lists using DESeq2
# Author: Gisela Gabernet
# QBiC 2019; MIT License

library(gprofiler2)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(pathview)
library(AnnotationDbi)
library(optparse)

# Need to load library for your species
library(org.Mm.eg.db) #Mmusculus
library(org.Hs.eg.db) #Hsapiens

# Blacklist pathways: some pathways are corrupted in KEGG and produce errors. Add the pathway here if you have this kind of error:
blacklist_pathways <- c("mmu05206", "mmu04215", "hsa05206", "mmu04723")

# ------------------
# Reading parameters
# ------------------

option_list = list(
  make_option(c("-c", "--dirContrasts"), type="character", default=".", help="directory with DE gene list for each contrast", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL, help="path to metadata table", metavar="character"),
  make_option(c("-d", "--model"), type="character", default=NULL, help="path to linear model file", metavar="character"),
  make_option(c("-n", "--normCounts"), type="character", default=NULL, help="path to normalized counts", metavar="character"),
  make_option(c("-s", "--species"), type="character", default=NULL, help="Species name. Example format: Hsapiens", metavar="character"),
  make_option(c("-g", "--genelist"), type="character", default=NULL, help="Path to gene list for heatmap plot.", metavar="character"),
  make_option(c("-b", "--kegg_blacklist"), type="character", default=NULL, help="Path to KEGG pathway blacklist.", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Metadata
if (is.null(opt$metadata)){
  print_help(opt_parser)
  stop("Metadata table needs to be provided!")
} else {
  metadata_path = opt$metadata
}

# Linear model
if (is.null(opt$model)){
  print_help(opt_parser)
  stop("Linear model file needs to be provided!")
} else {
  path_design = opt$model
}

# Count table
if(is.null(opt$normCounts)){
  print_help(opt_parser)
  stop("Normalized counts file needs to be provided")
} else {
  path_norm_counts = opt$normCounts
}

# Pathway blacklist append if provided
if(!is.null(opt$KEGG_blacklist)){
  KEGG_blacklist <- read.table(file=KEGG_blacklist, sep = "\t", header = F, quote="")
  blacklist_pathways <- append(blacklist_pathways, KEGG_blacklist$V1)
}

# Need to provide short and long names and libraries for your species
if(is.null(opt$species)){
  print_help(opt_parser)
  stop("Species needs to be provided")
} else if (tolower(opt$species) == "hsapiens") {
  organism <- "hsapiens"
  short_organism_name <- "hsa"
  library <- org.Hs.eg.db
} else if (tolower(opt$species) == "mmusculus") {
  organism <- "mmusculus"
  short_organism_name <- "mmu"
  library <- org.Mm.eg.db
} else {
  stop("Species name is unknown, check for typos or contact responsible person to add your species.")
}

# Contrast files
contrast_files <- list.files(path=opt$dirContrasts)
path_contrasts <- opt$dirContrasts

# Reading count table
norm_counts <- read.table(file = path_norm_counts, header = T, row.names = 1, sep = "\t", quote = "")

# Reading metadata table
metadata <- read.table(file=metadata_path, sep = "\t", header = T, quote="")
metadata$Secondary.Name <- gsub(" ; ", "_", metadata$Secondary.Name) # Remove blank spaces and ; from secondary name
metadata$Secondary.Name <- gsub(" ", "_", metadata$Secondary.Name) # Remove blank spaces if there from secondary name


# Create output directory
outdir <- "pathway_analysis"
dir.create(outdir)
pathway_heatmaps_dir <- "pathway_heatmaps"
kegg_pathways_dir <- "KEGG_pathways"

# ------------------
# Set default params
# ------------------

# gprofiler pathway / term sources parameters
datasources <- c("KEGG", "REAC")

# Set theme for graphs
theme_set(theme_classic())

# ----------------------
# Start pathway analysis
# ----------------------

# For each contrast do pathway analysis

# Defining summary variables
contrast <- c()
number_DE_genes <- c()
number_enriched_pathways <- c()
list_DE_genes <- list()
list_enriched_pathways <- list()

for (file in contrast_files){
  #Reading DE genes list
  fname <- tools::file_path_sans_ext(basename(file))
  
  dir.create(paste(outdir, fname, sep="/"))
  dir.create(paste(outdir, fname, pathway_heatmaps_dir, sep="/"))
  dir.create(paste(outdir, fname, kegg_pathways_dir, sep="/"))
  
  DE_genes <- read.csv(file = paste0(path_contrasts, file), sep="\t", header = T)
  DE_genes <- as.data.frame(DE_genes)

  # Skip pathway analysis for the contrast if not 2 or more DE genes were found
  if (nrow(DE_genes) < 2){ 
    print(paste0("Not enough DE genes to allow for a pathway analysis for contrast: ", fname))
    next
  }

  # Define list of Ensemble IDs (q) to run the pathway analysis
  q = as.character(DE_genes$Ensembl_ID)

  # gost query
  gostres <- gost(query=q,
                  organism=organism,
                  significant=TRUE,
                  correction_method="fdr",
                  sources=datasources,
                  evcodes=TRUE,
                  user_threshold=0.05)

  # Make data frame of gost result
  pathway_gostres <- gostres$result
  # Select only significantly enriched pathways (according to adjusted p-value)
  pathway_gostres <- as.data.frame(pathway_gostres[which(pathway_gostres$significant==TRUE),])

  # Plot pathways if there were any
  if (nrow(pathway_gostres) > 0){

    # annotate query size (number of DE genes in contrast)
    pathway_gostres$original_query_size <- rep(length(q), nrow(pathway_gostres))

    # Generate non-interactive pathway dotplots in the folder
    pg <- gostplot(gostres, capped=T, interactive=F)
    ggsave(pg, filename = paste0(outdir, "/", fname, "_gost_pathway_enrichment_plot.pdf"), 
          device="pdf", 
          height=10, width=15, units="cm", limitsize=F)
    ggsave(pg, filename = paste0(outdir, "/", fname, "_gost_pathway_enrichment_plot.png"), 
          device="png", 
          height=10, width=15, units="cm", dpi=300, limitsize=F)
  }

  # Remove parents column to be able to save the table in tsv format
  pathway_gostres_table = pathway_gostres
  pathway_gostres_table$parents <- NULL

  # Save pathway enrichment table in tsv format
  write.table(pathway_gostres_table, 
              file = paste0(outdir, "/", fname, "/", fname, "_pathway_enrichment_results.tsv"), 
              sep="\t", quote = F, col.names = T, row.names = F)

  # Collecting summary variables
  contrast <- append(contrast, fname)
  number_DE_genes <- append(number_DE_genes, length(DE_genes$Ensembl_ID))
  number_enriched_pathways <- append(number_enriched_pathways, summary(as.factor(pathway_gostres_table$source)))
  list_DE_genes <- append(list_DE_genes, DE_genes$Ensembl_ID)
  list_enriched_pathways <- append(list_enriched_pathways, pathway_gostres_table$term_id)

  # Printing summary variables
  print("------------------------------------")
  print(fname)
  print("Number of genes in query:")
  print(length(DE_genes$Ensembl_ID))
  print("Number of pathways found:")
  print(summary(as.factor(pathway_gostres_table$source)))
  print("------------------------------------")
  
  if (nrow(pathway_gostres) > 0){ #if there are enriched pathways
    # Splitting results according to pathway resources (KEGG / REACTOME / GO)
    res <- split(pathway_gostres, pathway_gostres$source)
    for (df in res){
      db_source <- df$source[1]
      df$short_name <- sapply(df$term_name, substr, start=1, stop=50)

      # Plotting results for df
      df_subset <- data.frame(Pathway_name = df$short_name, Pathway_code = df$term_id, DE_genes = df$intersection_size, Pathway_size = df$term_size, Fraction_DE = (df$intersection_size / df$term_size), Padj = df$p_value)
      write.table(df_subset, 
              file = paste0(outdir, "/", fname, "/", fname, "_", db_source, "_pathway_enrichment_results.tsv"), 
              sep="\t", quote = F, col.names = T, row.names = F)

      # Enriched pathways horizontal barplots of padj values
      p <- ggplot(df_subset, aes(x=reorder(Pathway_name, Fraction_DE), y=Fraction_DE)) +
        geom_bar(aes(fill=Padj), stat="identity", width = 0.7) +
        geom_text(aes(label=paste0(df_subset$DE_genes, "/", df_subset$Pathway_size)), vjust=0.4, hjust=-0.5, size=3) +
        coord_flip() +
        scale_y_continuous(limits = c(0.00, 1.00)) +
        scale_fill_continuous(high = "#132B43", low = "#56B1F7") +
        ggtitle("Enriched pathways") +
        xlab("") + ylab("Gene fraction (DE genes / Pathway size)")
      ggsave(p, filename = paste0(outdir, "/", fname, "/", fname, "_", db_source, "_pathway_enrichment_plot.pdf"), device = "pdf", height = 2+0.5*nrow(df_subset), units = "cm", limitsize=F)
      ggsave(p, filename = paste0(outdir, "/", fname, "/", fname,"_", db_source, "_pathway_enrichment_plot.png"), device = "png", height = 2+0.5*nrow(df_subset), units = "cm", dpi = 300, limitsize=F)

      # Plotting heatmaps and KEGG pathways for all pathways
      print("Plotting heatmaps...")
      if (nrow(df) <= 100 & nrow(df) > 0) {
        conditions <- grepl("Condition", colnames(metadata))
        metadata_cond <- as.data.frame(metadata[,conditions])
        metadata_name <- metadata[,c("QBiC.Code", "Secondary.Name")]
        row.names(metadata_cond) <- apply(metadata_name,1,paste, collapse = "_")

        for (i in c(1:nrow(df))){
          pathway <- df[i,]
          gene_list <- unlist(strsplit(pathway$intersection, ","))
          mat <- norm_counts[gene_list, ]
          rownames(mat) <- mat$gene_name
          mat$gene_name <- NULL
          mat <- data.matrix(mat)

          if (nrow(mat)>1){
            png(filename = paste(outdir, "/",fname, "/", pathway_heatmaps_dir, "/", "Heatmap_normalized_counts_", pathway$source, "_", pathway$term_id, "_",fname, ".png", sep=""), width = 2500, height = 3000, res = 300)
            pheatmap(mat = mat, annotation_col = metadata_cond, main = paste(pathway$short_name, "(",pathway$source,")",sep=" "), scale = "row", cluster_cols = F, cluster_rows = T )
            dev.off()

            pdf(paste(outdir, "/", fname, "/", pathway_heatmaps_dir, "/", "Heatmap_normalized_counts_", pathway$source, "_", pathway$term_id, "_", fname, ".pdf", sep=""))
            pheatmap(mat = mat, annotation_col = metadata_cond, main = paste(pathway$short_name, "(",pathway$source,")",sep=" "), scale = "row", cluster_cols = F, cluster_rows = T )
            dev.off()
          }

          # Plotting pathway view only for kegg pathways
          if (pathway$source == "KEGG"){
            pathway_kegg <- sapply(pathway$term_id, function(x) paste0(short_organism_name, unlist(strsplit(as.character(x), ":"))[2]))
            
            # Avoid KEGG pathways in blacklist. This pathway graphs contain errors and pathview crashes if plotting them.
            if (pathway_kegg %in% blacklist_pathways) {
              print(paste0("Skipping pathway: ",pathway_kegg,". This pathway file has errors in KEGG database."))
            } else {
              print(paste0("Plotting pathway: ", pathway_kegg))
              gene.data = DE_genes
              gene.data.subset = gene.data[gene.data$Ensembl_ID %in% gene_list, c("Ensembl_ID","log2FoldChange")]
              
              entrez_ids = mapIds(library, keys=as.character(gene.data.subset$Ensembl_ID), column = "ENTREZID", keytype="ENSEMBL", multiVals="first")
              
              gene.data.subset <- gene.data.subset[!(is.na(entrez_ids)),]

              if (length(entrez_ids)!=length(unique(entrez_ids))) {
                print(paste0("Skipping pathway: ", pathway_kegg,". This pathway has multiple IDs with same name."))
              } else {
                row.names(gene.data.subset) <- entrez_ids[!is.na(entrez_ids)]
                
                gene.data.subset$Ensembl_ID <- NULL
                pathview(gene.data  = gene.data.subset,
                        pathway.id = pathway_kegg,
                        species    = short_organism_name,
                        out.suffix=paste(fname,sep="_"))
                mv_command <- paste0("mv *.png *.xml ","./",outdir, "/",fname, "/", kegg_pathways_dir, "/")
                rm_command <- paste0("rm ","./",outdir, "/",fname, "/", kegg_pathways_dir, "/", "*.xml")
                system(mv_command)
              }
            }
          }
        }
      }
    }
  }
}

# Writing DE gene summary table
df_summary <- data.frame(Contrast = contrast, Number_DE_genes = number_DE_genes, Number_enriched_pathways = number_enriched_pathways, DE_genes_list = list_DE_genes, Enriched_pathways_list = list_enriched_pathways)
write.table(df_summary, 
        file = paste0(outdir, "/", fname, "/", fname, "_pathway_enrichment_summary.tsv"), 
        sep="\t", quote = F, col.names = T, row.names = F)


# Plotting heatmap for provided gene list

if (!is.null(opt$genelist)){

  genelist_path = opt$genelist

  genelist_heatmaps_dir <- "heatmap_gene_list"
  dir.create(paste(outdir, genelist_heatmaps_dir, sep="/"))


  print("Plotting heatmaps...")
  conditions <- grepl("Condition", colnames(metadata))
  condition <- metadata[,conditions]
  metadata_cond <- as.data.frame(condition)
  metadata_name <- metadata[,c("QBiC.Code", "Secondary.Name")]
  row.names(metadata_cond) <- apply(metadata_name,1,paste, collapse = "_")

  gene_list_tab <- read.table(file=genelist_path, sep = "\t", header = F, quote="")
  gene_list_unique_tab <- data.frame(gene_list=unique(gene_list_tab$V1))

  norm_counts <- read.table(file = path_norm_counts, header = T, sep = "\t", quote = "")
  norm_counts$gene_name <- toupper(norm_counts$gene_name)

  IDs <- norm_counts[,c("Ensembl_ID","gene_name")]

  genestoEnsmbl <- merge(x=gene_list_unique_tab, y=IDs, by.x="gene_list", by.y="gene_name", all.x=T)
  gene_list <- genestoEnsmbl$Ensembl_ID
  
  # Omit genes not present in the count table (NAs)
  gene_list <- na.omit(gene_list)

  rownames(norm_counts) <- norm_counts$Ensembl_ID

  mat <- norm_counts[gene_list, ]
  rownames(mat) <- mat$gene_name
  mat$gene_name <- NULL
  mat$Ensembl_ID <- NULL
  mat <- data.matrix(mat)

  if (nrow(mat)>1){
    png(filename = paste(outdir, "/", genelist_heatmaps_dir, "/", "Heatmap_normalized_counts_gene_list.png", sep=""), width = 2500, height = 3000, res = 300)
    pheatmap(mat = mat, annotation_col = metadata_cond, main = "", scale = "row", cluster_cols = F, cluster_rows = T )
    dev.off()

    pdf(paste(outdir, "/", genelist_heatmaps_dir, "/", "Heatmap_normalized_counts_gene_list.pdf", sep=""))
    pheatmap(mat = mat, annotation_col = metadata_cond, main = "", scale = "row", cluster_cols = F, cluster_rows = T )
    dev.off()
  }
}
