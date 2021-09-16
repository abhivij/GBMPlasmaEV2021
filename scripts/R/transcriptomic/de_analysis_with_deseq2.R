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

contrast_vec <- c("PREOPE", "MET")
get_comparison_results <- function(contrast_vec, comparison_num){
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

  

get_comparison_results(c("PREOPE", "MET"), 2)
get_comparison_results(c("PREOPE", "HC"), 2)
get_comparison_results(c("MET", "HC"), 2)




