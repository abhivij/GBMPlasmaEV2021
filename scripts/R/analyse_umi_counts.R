library(tidyverse)
library(readxl)
library(edgeR)
library(caret)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)
filepath <- paste(base_dir, "mirna_counts.csv", sep = "/")
generated_mirna_counts <- read.table(filepath, header=TRUE, sep = ",", row.names = 1)

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


# read_counts <- data %>%
#   select(ends_with("READs"))
# 
# 
# split_function <- function(x){
#   return (strsplit(x, split = "-", fixed = TRUE)[[1]][1])
# }
# 
# length(colnames(read_counts))
# length(colnames(umi_counts))
# sum(sapply(colnames(umi_counts), split_function) == sapply(colnames(read_counts), split_function))
# 
# 
# data_sub <- data[1:10, c(1:5, 122:126)]
# new_names_prefix <- paste("Sample", c(1:5), sep = "")
# new_names <- c(paste(new_names_prefix, "UMIs", sep = "-"),
#                paste(new_names_prefix, "READs", sep = "-"))
# colnames(data_sub) <- new_names



plot_data <- function(norm_data, filename, title, width = 30, height = 30, perplexity = 30){
  all_group_names <- sapply(colnames(norm_data),
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
  all_sample_names <- colnames(norm_data)
  
  
  
  ne <- sapply(norm_data, function(x) as.numeric(x))  #causes rownames to be lost
  ne <- t(ne)
  
  pca_file_name <- paste("pca", filename, sep = "_")
  p <- prcomp(ne, scale. = TRUE, center = TRUE)
  pca_plotdata <- data.frame(p$x) %>%
    mutate(type = all_group_names) %>%
    mutate(name = all_sample_names)
  # 
  tsne_file_name <- paste("tsne", filename, sep = "_")
  set.seed(1)
  tsne_result <- Rtsne::Rtsne(ne, perplexity = perplexity)
  tsne_df <- data.frame(x = tsne_result$Y[,1], y = tsne_result$Y[,2], 
                        Colour = all_group_names,
                        name = all_sample_names)
  
  pca_plotdata %>%
    ggplot(aes(x = PC1, y = PC2, colour = type, label = name)) +
    geom_point(size = 2) +
    labs(title = title)
  ggsave(pca_file_name, width = 10, height = 10)
  
  tsne_plot <- tsne_df %>%
    ggplot(aes(x = x, y = y, colour = Colour, label = name)) +
    geom_point(size = 2) +
    ggplot2::labs(colour = "Groups", title = title) +
    ggplot2::xlab("tSNE 1") +
    ggplot2::ylab("tSNE 2")
  ggplot2::ggsave(tsne_file_name, tsne_plot)
  
}

plot_data(umi_counts, "umi_counts.jpg", "UMIs")


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

keep <- filterByExpr(umi_counts, group = all_group_names)
umi_counts <- umi_counts[keep, ]

plot_data(umi_counts, "filtered_umi_counts.jpg", "Filtered UMIs")

#unable to run below line
plot_data(cpm(umi_counts), "cpm_umi_counts.jpg", "CPM UMIs", perplexity = 2)


umi_norm_data <- cpm(umi_counts, log = TRUE)
umi_norm_data <- scale(umi_norm_data)

umi_norm_data <- as.data.frame(t(as.matrix(umi_norm_data)))
umi_norm_data <- predict(preProcess(umi_norm_data), umi_norm_data)
umi_norm_data <- as.data.frame(t(as.matrix(umi_norm_data)))

plot_data(umi_norm_data, "normlogcpm_umi_counts.jpg", "Norm Log CPM UMIs")


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
  geom_histogram(aes(x = umi_counts[,5]), stat = "bin", bins = 200) +
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



