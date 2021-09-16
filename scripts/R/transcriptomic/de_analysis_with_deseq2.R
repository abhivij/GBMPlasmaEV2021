library(tidyverse)
library(readxl)
library(edgeR)
library(caret)
library("DESeq2")


base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

source("scripts/R/utils.R")
source("scripts/R/plot_data.R")

umi_counts <- read.csv("Data/RNA/umi_counts.csv", row.names = 1)
metadata <- read.csv("Data/metadata.csv") %>%
  column_to_rownames("SUBJECT_ORIGINAL")
umi_counts <- umi_counts[, rownames(metadata)]

all(rownames(metadata) == colnames(umi_counts))

dds <- DESeqDataSetFromMatrix(countData = umi_counts,
                              colData = metadata,
                              design = ~ GROUP_Q1to6)
dds

print(dim(counts(dds)))

# keep <- rowSums(counts(dds)) >= 30
# dds <- dds[keep,]

print(dim(counts(dds)))

dds <- DESeq(dds)

contrast_vec <- c("MET", "HC")
comparison_num <- 2
get_comparison_results <- function(contrast_vec, comparison_num, umi_counts, metadata){
  res <- results(dds, contrast = c("GROUP_Q1to6", contrast_vec[1], contrast_vec[2]), alpha = 0.05)
  res
  print(summary(res))
  
  result <- data.frame(res) %>%
    filter(!is.na(padj)) %>%
    rownames_to_column("rna") %>%
    select(rna, log2FoldChange, pvalue, padj) %>%
    dplyr::rename(Molecule = rna, adjPVal = padj, pVal = pvalue, logFC = log2FoldChange)
  
  title <- paste(contrast_vec, collapse = " Vs ")
  
  file_name <- paste("de", "comp", comparison_num, 
                     paste(contrast_vec, collapse = "_"), sep = "_")
  file_name <- paste(file_name, "png", sep = ".")
  
  
  plot_transcriptomic_data(umi_counts[(result %>% filter(adjPVal < 0.05))$Molecule, ], 
                           metadata, contrast_vec, "Adj pval < 0.05 transcripts",
                           file_name = "adjpval")
  plot_transcriptomic_data(umi_counts[(result %>% filter(adjPVal < 0.05))$Molecule, ], 
                           metadata, contrast_vec, "Adj pval < 0.05 transcripts", dim_red = "umap",
                           file_name = "adjpval")
  
  plot_transcriptomic_data(umi_counts[(result %>% filter(pVal < 0.05))$Molecule, ], 
                           metadata, contrast_vec, "pval < 0.05 transcripts",
                           file_name = "pval")
  plot_transcriptomic_data(umi_counts[(result %>% filter(pVal < 0.05))$Molecule, ], 
                           metadata, contrast_vec, "pval < 0.05 transcripts", dim_red = "umap",
                           file_name = "pval")
  
  
  create_volcano_plot(result, 
                      title = title, 
                      file_name = paste("volcano", file_name, sep = "_"), 
                      dir_path = "plots/de/transcriptomic/deseq2", logFC_cutoff = 1)
  create_volcano_plot(result, 
                      title = title, 
                      file_name = paste("volcano_pval", file_name, sep = "_"), 
                      dir_path = "plots/de/transcriptomic/deseq2", logFC_cutoff = 1,
                      use_p_val = TRUE)  
  
}

  

get_comparison_results(c("PREOPE", "MET"), 2, umi_counts, metadata)
get_comparison_results(c("PREOPE", "HC"), 2, umi_counts, metadata)
get_comparison_results(c("MET", "HC"), 2, umi_counts, metadata)




plot_transcriptomic_data <- function(data, metadata, categories, title, dim_red = "pca",
                                     file_name = ""){
  metadata <- metadata %>%
    filter(GROUP_Q1to6 %in% categories)
  print("metadata size")
  print(dim(metadata))
  
  data <- data[, rownames(metadata)]
  print(all(rownames(metadata) == colnames(data)))
  groups <- metadata$GROUP_Q1to6
  
  file_name <- paste("umi_counts", file_name, paste(categories, collapse = "_"), sep = "_")
  file_name <- paste(file_name, "png", sep = ".")
  
  plot_data(t(data), file_name, title, groups = groups, dim_red = dim_red,
            colour_label = "Labels")
}

plot_transcriptomic_data(umi_counts, metadata, c("PREOPE", "MET"), "All transcripts")
plot_transcriptomic_data(umi_counts, metadata, c("PREOPE", "MET"), "All transcripts", dim_red = "umap")

plot_transcriptomic_data(umi_counts, metadata, c("PREOPE", "HC"), "All transcripts")
plot_transcriptomic_data(umi_counts, metadata, c("PREOPE", "HC"), "All transcripts", dim_red = "umap")

plot_transcriptomic_data(umi_counts, metadata, c("MET", "HC"), "All transcripts")
plot_transcriptomic_data(umi_counts, metadata, c("MET", "HC"), "All transcripts", dim_red = "umap")
