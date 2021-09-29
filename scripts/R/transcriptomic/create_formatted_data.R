library(tidyverse)
library(readxl)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)


data <- read_excel("Data/RNA/158629.all_samples.summary.xlsx", sheet = "miRNA_piRNA")
mirna_data <- data[1:2505,]
pirna_data <- data[2507:2642,] %>%
  separate(miRNA, c("miRNA", NA, NA, NA), sep = "/")

dim(data)[1]
dim(mirna_data)[1] + dim(pirna_data)[1]

data <- rbind(mirna_data, pirna_data) %>%
  column_to_rownames("miRNA")

umi_counts <- data %>%
  select(ends_with("UMIs"))
colnames(umi_counts) <- gsub("-UMIs", "", colnames(umi_counts))
colnames(umi_counts) <- sapply(colnames(umi_counts), FUN = 
                                 function(x){
                                   strsplit(x, split = "_", fixed = TRUE)[[1]][1]
                                 }
)

write.csv(umi_counts, "Data/RNA/umi_counts.csv")
