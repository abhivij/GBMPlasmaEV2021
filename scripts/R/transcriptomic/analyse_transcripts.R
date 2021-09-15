library(tidyverse)
library(readxl)
library(edgeR)
library(caret)


base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

source("scripts/R/utils.R")
source("scripts/R/plot_data.R")


generated_mirna_counts <- read.table("mirna_counts.csv", header=TRUE, sep = ",", row.names = 1)

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
all_group_names <- sapply(colnames(umi_counts),
                          function(x){
                            prefix <- substr(x, 1, 2)
                            if (prefix == "ME") {
                              label <- "MET"
                            } else if (prefix == "HB") {
                              label <- "GBM"
                            } else {
                              label <- "HC"
                            }
                            return (label)
                          }
)

plot_transcriptomic_data <- function(data, file_name, title, groups = NA, dim_red = "pca"){
  if(length(groups) <= 1 && is.na(groups)){
    groups <- sapply(colnames(data),
                              function(x){
                                prefix <- substr(x, 1, 2)
                                if (prefix == "ME") {
                                  label <- "MET"
                                } else if (prefix == "HB") {
                                  label <- "GBM"
                                } else {
                                  label <- "HC"
                                }
                                return (label)
                              }
    )  
  }
  plot_data(t(data), file_name, title, groups = groups, dim_red = dim_red,
            colour_label = "Labels")
}

plot_transcriptomic_data(umi_counts, "umi_counts.jpg", "UMIs")
plot_transcriptomic_data(umi_counts, "umi_counts.jpg", "UMIs", dim_red = "tsne")
plot_transcriptomic_data(umi_counts, "umi_counts.jpg", "UMIs", dim_red = "umap")




keep <- filterByExpr(umi_counts, group = all_group_names)
umi_counts <- umi_counts[keep, ]

plot_transcriptomic_data(umi_counts, "filtered_umi_counts.jpg", "Filtered UMIs")
plot_transcriptomic_data(umi_counts, "filtered_umi_counts.jpg", "Filtered UMIs", dim_red = "tsne")
plot_transcriptomic_data(umi_counts, "filtered_umi_counts.jpg", "Filtered UMIs", dim_red = "umap")


umi_norm_data <- cpm(umi_counts, log = TRUE)
umi_norm_data <- scale(umi_norm_data)

umi_norm_data <- as.data.frame(t(as.matrix(umi_norm_data)))
umi_norm_data <- predict(preProcess(umi_norm_data), umi_norm_data)
umi_norm_data <- as.data.frame(t(as.matrix(umi_norm_data)))

plot_transcriptomic_data(umi_norm_data, "normlogcpm_umi_counts.jpg", "Norm Log CPM UMIs")
plot_transcriptomic_data(umi_norm_data, "normlogcpm_umi_counts.jpg", "Norm Log CPM UMIs", dim_red = "tsne")
plot_transcriptomic_data(umi_norm_data, "normlogcpm_umi_counts.jpg", "Norm Log CPM UMIs", dim_red = "umap")


all_group_names <- sapply(colnames(umi_norm_data),
                          function(x){
                            prefix <- substr(x, 1, 2)
                            if (prefix == "ME") {
                              label <- "MET"
                            } else if (prefix == "HB") {
                              label <- "GBM"
                            } else {
                              label <- "HC"
                            }
                            return (label)
                          }
)

ggplot() +
  geom_histogram(aes(x = umi_norm_data[,3]), stat = "bin", bins = 200) +
  xlab("UMIs") +
  ylab("Normalized counts")
ggsave("umi_hist.jpg", width = 10, height = 10)


ggplot() +
  geom_histogram(aes(x = umi_counts[,15]), stat = "bin", bins = 200) +
  xlab("UMIs") +
  ylab("Normalized counts")
ggsave("umi_hist.jpg", width = 10, height = 10)



###################

plot_data(generated_mirna_counts, "read_counts.jpg", "Read counts")


all_group_names <- sapply(colnames(generated_mirna_counts),
                          function(x){
                            prefix <- substr(x, 1, 2)
                            if (prefix == "ME") {
                              label <- "MET"
                            } else if (prefix == "HB") {
                              label <- "GBM"
                            } else {
                              label <- "HC"
                            }
                            return (label)
                          }
)
keep <- filterByExpr(generated_mirna_counts, group = all_group_names)
generated_mirna_counts <- generated_mirna_counts[keep, ]
plot_data(generated_mirna_counts, "filtered_read_counts.jpg", "Filtered Reads")


norm_data <- cpm(generated_mirna_counts, log = TRUE)
norm_data <- scale(norm_data)

norm_data <- as.data.frame(t(as.matrix(norm_data)))
norm_data <- predict(preProcess(norm_data), norm_data)
norm_data <- as.data.frame(t(as.matrix(norm_data)))

plot_data(norm_data, "normlogcpm_read_counts.jpg", "Norm Log CPM Reads")


