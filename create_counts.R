library(tidyverse)
library(tximport)
library(biomaRt)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
quant_dir_name <- "quants_k21"

setwd(base_dir)

files <- list.files(quant_dir_name, recursive = TRUE, pattern = "*.sf")
samples <- sapply(list.files(quant_dir_name), FUN = 
                   function(x){
                     return(strsplit(x, fixed = TRUE, split = "_")[[1]][1])
                   }) 
names(files) <- samples

quant_dir <- paste(base_dir, quant_dir_name, sep = "/")
setwd(quant_dir)

print(paste("all file exists ?", all(file.exists(files))))


txi <- tximport(files, type = "salmon", txOut = TRUE)
counts <- txi$counts
counts <- as.data.frame(counts) %>%
  dplyr::filter(rowSums(counts) > 0)



#get miRNA-ensemblid mapping from ensembl biomart
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
miRNA <- getBM(attributes=c("mirbase_id","ensembl_transcript_id"),
              filters=c("with_mirbase"),values=list(TRUE), mart=ensembl)


counts_modified <- counts %>%
  mutate(id = rownames(counts)) %>%
  separate("id", c("ensembl_transcript_id", NA))

length(counts_modified$ensembl_transcript_id)
length(unique(counts_modified$ensembl_transcript_id))

miRNA_counts <- counts_modified %>%
  inner_join(miRNA)
length(miRNA_counts$ensembl_transcript_id)
length(miRNA_counts$mirbase_id)
length(unique(miRNA_counts$ensembl_transcript_id))
length(unique(miRNA_counts$mirbase_id))


counts_filtered_for_miRNA <- counts_modified %>%
  filter(ensembl_transcript_id %in% unique(miRNA_counts$ensembl_transcript_id))
rownames(counts_filtered_for_miRNA) <- counts_filtered_for_miRNA$ensembl_transcript_id
counts_filtered_for_miRNA <- counts_filtered_for_miRNA %>%
  select(-c(ensembl_transcript_id))

filename <- paste(quant_dir_name, "transcript_counts.csv", sep = "_")
filepath <- paste("..", filename, sep = "/")
write.table(x = counts, file = filepath, quote = FALSE, sep = ",", row.names=TRUE, col.names=NA)

filename <- paste(quant_dir_name, "mirna_counts.csv", sep = "_")
filepath <- paste("..", filename, sep = "/")
write.table(x = counts_filtered_for_miRNA, file = filepath, quote = FALSE, sep = ",", row.names=TRUE, col.names=NA)


miRNA_subset <- miRNA %>%
  filter(ensembl_transcript_id %in% unique(miRNA_counts$ensembl_transcript_id))
length(unique(miRNA_subset$mirbase_id))
length(unique(miRNA_subset$ensembl_transcript_id))

filename <- paste("mirna_name_mapping.csv", sep = "_")
write.table(x = miRNA_subset, file = filename, quote = FALSE, sep = ",", row.names = FALSE)
