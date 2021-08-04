library(tidyverse)
library(readxl)

setwd("/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV")

data <- read_excel("Data/RNA/158629.all_samples.summary.xlsx", sheet = "miRNA_piRNA")

data <- data %>%
  column_to_rownames("miRNA")

umi_counts <- data %>%
  select(ends_with("UMIs"))

read_counts <- data %>%
  select(ends_with("READs"))


split_function <- function(x){
  return (strsplit(x, split = "-", fixed = TRUE)[[1]][1])
}

length(colnames(read_counts))
length(colnames(umi_counts))
sum(sapply(colnames(umi_counts), split_function) == sapply(colnames(read_counts), split_function))


data_sub <- data[1:10, c(1:5, 122:126)]
new_names_prefix <- paste("Sample", c(1:5), sep = "")
new_names <- c(paste(new_names_prefix, "UMIs", sep = "-"),
               paste(new_names_prefix, "READs", sep = "-"))
colnames(data_sub) <- new_names

