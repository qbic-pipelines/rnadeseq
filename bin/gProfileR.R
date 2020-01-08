#!/usr/bin/env Rscript

# Pathway analysis from differentially expressed gene lists using DESeq2
# Author: Gisela Gabernet
# QBiC 2019; MIT License

library("gProfileR")
library("ggplot2")
library("reshape2")
library("pheatmap")
library("pathview")
library("AnnotationDbi")
library("optparse")

# Need to load library for your species
library(org.Mm.eg.db) #Mmusculus
library(org.Hs.eg.db) #Hsapiens

# Reading parameters

option_list = list(
  make_option(c("-c", "--dirContrasts"), type="character", default=".", help="directory with DE gene list for each contrast", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL, help="path to metadata table", metavar="character"),
  make_option(c("-d", "--model"), type="character", default=NULL, help="path to linear model file", metavar="character"),
  make_option(c("-n", "--normCounts"), type="character", default=NULL, help="path to normalized counts", metavar="character"),
  make_option(c("-s", "--species"), type="character", default=NULL, help="Species name. Example format: Hsapiens", metavar="character")
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
} else if (opt$species == "Hsapiens") {
  organism <- "hsapiens"
  short_organism_name <- "hsa"
  library <- org.Hs.eg.db
} else if (opt$species == "Mmusculus") {
  organism <- "mmusculus"
  short_organism_name <- "mmu"
  library <- org.Mm.eg.db
} else {
  stop("Species name is unknown, check for typos or contact responsible person to add your species.")
}

contrast_files <- list.files(path=opt$dirContrasts)
path_contrasts <- opt$dirContrasts

outdir <- "gProfileR"

norm_counts <- read.table(file = path_norm_counts, header = T, row.names = 1, sep = "\t", quote = "")
metadata <- read.table(file=metadata_path, sep = "\t", header = T, quote="")


#Search params

datasources <- c("KEGG","REAC")
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
  path_enrich <- gprofiler(query = q, organism=organism, 
                           significant = T, correction_method = "fdr",
                           min_set_size = min_set_size, max_set_size = max_set_size, min_isect_size = min_isect_size,
                           src_filter = datasources)
  
  if (nrow(path_enrich) > 0){
    path_enrich$original.query.size <- rep(length(q), nrow(path_enrich))
  }
  write.table(path_enrich, file = paste0(outdir, "/", fname, "/",fname, "_pathway_enrichment_results.tsv"), sep = "\t", quote = F, col.names = T, row.names = F )
  
  # Printing numbers
  print("------------------------------------")
  print(fname)
  print("Number of genes in query:")
  print(length(DE_genes$Ensembl_ID))
  print("Number of pathways found:")
  print(summary(as.factor(path_enrich$domain)))
  print("------------------------------------")
  
  if (nrow(path_enrich) > 0){ #if there are enriched pathways
    # Splitting results according to tools
    res <- split(path_enrich, path_enrich$domain)
    for (df in res){
      db_source <- df$domain[1]
      print(db_source)
      df$short_name <- sapply(df$term.name, substr, start=1, stop=50)

      # Plotting results for df
      df_subset <- data.frame(Pathway_name = df$short_name, Query = df$overlap.size, Pathway = df$term.size, Fraction = (df$overlap.size / df$term.size), Pval = df$p.value)

      p <- ggplot(df_subset, aes(x=reorder(Pathway_name, Fraction), y=Fraction)) +
        geom_bar(aes(fill=Pval), stat="identity", width = 0.7) +
        geom_text(aes(label=paste0(df_subset$Query, "/", df_subset$Pathway)), vjust=0.4, hjust=-0.5, size=3) +
        coord_flip() +
        scale_y_continuous(limits = c(0.00, 1.00)) +
        scale_fill_continuous(high = "#132B43", low = "#56B1F7") +
        ggtitle("Enriched pathways") +
        xlab("") + ylab("Gene fraction (Query / Pathway)")
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

          if (nrow(mat)>1){
            png(filename = paste(outdir, "/",fname, "/", pathway_heatmaps_dir, "/", "Heatmap_normalized_counts_", pathway$domain, "_", pathway$term.id, "_",fname, ".png", sep=""), width = 2500, height = 3000, res = 300)
            tryCatch(withCallingHandlers(pheatmap(mat = mat, annotation_col = metadata_cond, main = paste(pathway$short_name, "(",pathway$domain,")",sep=" "), scale = "row", cluster_cols = F, cluster_rows = T ), 
                    error=function(e) {print(paste0("Skipping heatmap plot due to problem:"))},
                    warning=function(w) {print(paste0("Warning for heatmap plot."))
                                invokeRestart("muffleWarning")}), 
            error = function(e) { print(paste0("Heatmap plot finished")) })
            dev.off()
          }

          # Plotting pathway view only for kegg pathways
          if (pathway$domain == "keg"){
            pathway_kegg <- sapply(pathway$term.id, function(x) paste0(short_organism_name, unlist(strsplit(as.character(x), ":"))[2]))
            print(paste0("Plotting pathway: ", pathway_kegg))
            # KEGG pathway blacklist. This pathway graphs contain errors and pathview crashes if plotting them.
            if (pathway_kegg %in% c("mmu05206", "mmu04215", "hsa05206") ) {
              print(paste0("Skipping pathway: ",pathway_kegg,". This pathway file has errors in KEGG database."))
            } else {
              gene.data = DE_genes
              gene.data.subset = gene.data[gene.data$Ensembl_ID %in% gene_list, c("Ensembl_ID","log2FoldChange")]
              
              entrez_ids = mapIds(library, keys=as.character(gene.data.subset$Ensembl_ID), column = "ENTREZID", keytype="ENSEMBL", multiVals="first")
              
              gene.data.subset <- gene.data.subset[!(is.na(entrez_ids)),]
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



