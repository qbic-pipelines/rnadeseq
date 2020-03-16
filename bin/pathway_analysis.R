#!/usr/bin/env Rscript

# Pathway analysis from differentially expressed gene lists using DESeq2
# Author: Gisela Gabernet
# QBiC 2019; MIT License

library(gprofiler2)
library(gProfileR)
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

# Reading parameters

option_list = list(
  make_option(c("-c", "--dirContrasts"), type="character", default=".", help="directory with DE gene list for each contrast", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL, help="path to metadata table", metavar="character"),
  make_option(c("-d", "--model"), type="character", default=NULL, help="path to linear model file", metavar="character"),
  make_option(c("-n", "--normCounts"), type="character", default=NULL, help="path to normalized counts", metavar="character"),
  make_option(c("-s", "--species"), type="character", default=NULL, help="Species name. Example format: Hsapiens", metavar="character"),
  make_option(c("-g", "--genelist"), type="character", default=NULL, help="Gene list for heatmap plot.", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$metadata)){
  print_help(opt_parser)
  stop("Metadata table needs to be provided!")
} else {
  metadata_path = opt$metadata
}

if (is.null(opt$model)){
  print_help(opt_parser)
  stop("Linear model file needs to be provided!")
} else {
  path_design = opt$model
}

if(is.null(opt$normCounts)){
  print_help(opt_parser)
  stop("Normalized counts file needs to be provided")
} else {
  path_norm_counts = opt$normCounts
}

# Need to provide short and long names and libraries for your species
if(is.null(opt$species)){
  print_help(opt_parser)
  stop("Species needs to be provided")
} else if (opt$species == "hsapiens") {
  organism <- "hsapiens"
  short_organism_name <- "hsa"
  library <- org.Hs.eg.db
} else if (opt$species == "mmusculus") {
  organism <- "mmusculus"
  short_organism_name <- "mmu"
  library <- org.Mm.eg.db
} else {
  stop("Species name is unknown, check for typos or contact responsible person to add your species.")
}

contrast_files <- list.files(path=opt$dirContrasts)
path_contrasts <- opt$dirContrasts

outdir <- "pathway_analysis"

# Reading count table
norm_counts <- read.table(file = path_norm_counts, header = T, row.names = 1, sep = "\t", quote = "")

# Reading metadata table
metadata <- read.table(file=metadata_path, sep = "\t", header = T, quote="")
metadata$Secondary.Name <- gsub(" ; ", "_", metadata$Secondary.Name) # Remove blank spaces and ; from secondary name
metadata$Secondary.Name <- gsub(" ", "_", metadata$Secondary.Name) # Remove blank spaces if there from secondary name

# Search params
datasources <- c("KEGG", "REAC")
min_set_size <- 1
max_set_size <- 500
min_isect_size <- 1

# Create output directory
dir.create(outdir)
pathway_heatmaps_dir <- "pathway_heatmaps"
kegg_pathways_dir <- "KEGG_pathways"

# Set theme for graphs
theme_set(theme_classic())

for (file in contrast_files){
  #Reading DE genes list
  fname <- tools::file_path_sans_ext(basename(file))
  
  dir.create(paste(outdir, fname, sep="/"))
  dir.create(paste(outdir, fname, pathway_heatmaps_dir, sep="/"))
  dir.create(paste(outdir, fname, kegg_pathways_dir, sep="/"))
  
  DE_genes <- read.csv(file = paste0(path_contrasts, file), sep="\t", header = T)
  q = as.character(DE_genes$Ensembl_ID)
  
 
  #gprofiler query
  #path_enrich <- gprofiler(query = q, organism=organism, 
  #                        significant = T, 
  #                        correction_method = "fdr",
  #                        min_set_size = min_set_size, 
  #                        max_set_size = max_set_size, 
  #                        min_isect_size = min_isect_size,
  #                        src_filter = datasources)

  #gost query
  gostres <- gost(query=q,
                  organism=organism,
                  significant=TRUE,
                  correction_method="fdr",
                  sources=datasources,
                  evcodes=TRUE,
                  user_threshold=0.05)

  path_gostres<- gostres$result
  path_gostres <- path_gostres[which(path_gostres$significant==TRUE),]

  if (nrow(path_gostres) > 0){
    path_gostres$original.query.size <- rep(length(q), nrow(path_gostres))

    pg <- gostplot(gostres, capped=T, interactive=F)
    ggsave(pg, filename = paste0(outdir, "/", fname, "_gost_pathway_enrichment_plot.pdf"), 
          device="pdf", 
          height=10, width=15, units="cm", limitsize=F)
    ggsave(pg, filename = paste0(outdir, "/", fname, "_gost_pathway_enrichment_plot.png"), 
          device="png", 
          height=10, width=15, units="cm", dpi=300, limitsize=F)
  }
  
  if (nrow(path_gostres) > 0){
    path_gostres$original_query_size <- rep(length(q), nrow(path_gostres))
  }
  #write.table(path_enrich, 
  #            file = paste0(outdir, "/", fname, "/",fname, "_pathway_enrichment_results.tsv"), 
  #            sep = "\t", quote = F, col.names = T, row.names = F )
  path_gostres_table = path_gostres
  path_gostres_table$parents <- NULL
  write.table(path_gostres_table, 
              file = paste0(outdir, "/", fname, "/", fname, "_pathway_enrichment_results.tsv"), 
              sep="\t", quote = F, col.names = T, row.names = F)

  # Printing numbers
  print("------------------------------------")
  print(fname)
  print("Number of genes in query:")
  print(length(DE_genes$Ensembl_ID))
  print("Number of pathways found:")
  print(summary(as.factor(path_gostres_table$source)))
  print("------------------------------------")
  
  if (nrow(path_gostres) > 0){ #if there are enriched pathways
    # Splitting results according to tools
    res <- split(path_gostres, path_gostres$source)
    for (df in res){
      db_source <- df$source[1]
      df$short_name <- sapply(df$term_name, substr, start=1, stop=50)

      # Plotting results for df
      df_subset <- data.frame(Pathway_name = df$short_name, Pathway_code = df$term_id, DE_genes = df$intersection_size, Pathway_size = df$term_size, Fraction_DE = (df$intersection_size / df$term_size), Padj = df$p_value)
      write.table(df_subset, 
              file = paste0(outdir, "/", fname, "/", fname, "_", db_source, "_pathway_enrichment_results.tsv"), 
              sep="\t", quote = F, col.names = T, row.names = F)

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


      # Plotting heatmaps and pathways for all pathways
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
            # KEGG pathway blacklist. This pathway graphs contain errors and pathview crashes if plotting them.
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
                system(mv_command)
              }
            }
          }
        }
      }
    }
  }
}

# Plotting heatmap for provided gene list

if (!is.null(opt$metadata)){

  genelist_heatmaps_dir <- "heatmap_gene_list"
  dir.create(paste(outdir, genelist_heatmaps_dir, sep="/"))


  print("Plotting heatmaps...")
  conditions <- grepl("Condition", colnames(metadata))
  metadata_cond <- as.data.frame(metadata[,conditions])
  metadata_name <- metadata[,c("QBiC.Code", "Secondary.Name")]
  row.names(metadata_cond) <- apply(metadata_name,1,paste, collapse = "_")

  gene_list <- read.table(file=metadata_path, sep = "\t", header = F, quote="")
  gene_list <- unlist(gene_list)
  print(gene_list)

  mat <- norm_counts[gene_list, ]
  rownames(mat) <- mat$gene_name
  mat$gene_name <- NULL
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




