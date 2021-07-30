library(tidyverse)
library(tximport)
library(biomaRt)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
quant_dir_name <- "quants_subset"

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
  filter(rowSums(counts) > 0)



#get miRNA-ensemblid mapping from ensembl biomart
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
miRNA <- getBM(attributes=c("mirbase_id","ensembl_transcript_id"),
              filters=c("with_mirbase"),values=list(TRUE), mart=ensembl)



miRNA_counts <- counts %>%
  mutate(id = rownames(counts)) %>%
  separate("id", c("ensembl_transcript_id", NA)) %>%
  inner_join(miRNA) %>%
  column_to_rownames("mirbase_id") %>%
  select(-c(ensembl_transcript_id))



filename <- paste(quant_dir_name, "transcript_counts.csv", sep = "_")
filepath <- paste("..", filename, sep = "/")
write.table(x = counts, file = filepath, quote = FALSE, sep = ",", row.names=TRUE, col.names=NA)

filename <- paste(quant_dir_name, "mirna_counts.csv", sep = "_")
filepath <- paste("..", filename, sep = "/")
write.table(x = miRNA_counts, file = filepath, quote = FALSE, sep = ",", row.names=TRUE, col.names=NA)
