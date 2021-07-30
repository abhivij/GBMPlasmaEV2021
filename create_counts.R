library(tidyverse)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME") 

tx2gene <- tx2gene %>%
  drop_na("GENEID")


library(tximport)
setwd("/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV")

files <- list.files("quants_subset", recursive = TRUE, pattern = "*.sf")
samples <- sapply(list.files("quants_subset"), FUN = 
                   function(x){
                     return(strsplit(x, fixed = TRUE, split = "_")[[1]][1])
                   }) 
names(files) <- samples

setwd("/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV/quants_subset")
all(file.exists(files))


txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)

counts <- txi$counts

counts <- as.data.frame(counts) %>%
  filter(rowSums(counts) > 0)


write.table(x = counts, file = "counts_subset.txt", quote = FALSE, sep = "\t", row.names=TRUE, col.names=NA)


txi <- tximport(files, type = "salmon", txOut = TRUE)
names(txi)

counts <- txi$counts

counts <- as.data.frame(counts) %>%
  filter(rowSums(counts) > 0)

rnames <- data.frame(rownames(counts))
colnames(rnames) <- "rnames"
  
rnames <- rnames %>%
  separate("rnames", c("n", NA))

rnames <- sapply(rownames(counts), FUN = 
                    function(x){
                      return(strsplit(x, fixed = TRUE, split = ".")[[1]][1])
                    }) 


write.table(x = counts, file = "transcript_counts_subset.txt", quote = FALSE, sep = "\t", row.names=TRUE, col.names=NA)


library(biomaRt)


ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
miRNA <- getBM( attributes=c("mirbase_id","ensembl_transcript_id"),
                filters=c("with_mirbase"),values=list(TRUE), mart=ensembl)

