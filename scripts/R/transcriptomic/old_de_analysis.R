library(tidyverse)
library(readxl)
library(edgeR)
library(caret)


base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

source("scripts/R/utils.R")
source("scripts/R/plot_data.R")

umi_counts <- read.csv("Data/RNA/umi_counts.csv", row.names = 1)
metadata <- read.csv("Data/metadata.csv")


ggplot(umi_counts) +
  geom_histogram(aes(x = MET21), stat = "bin", bins = 200) +
  xlab("UMI counts") +
  ylab("Number of transcripts")

ggplot(umi_counts) +
  geom_histogram(aes(x = HB50), stat = "bin", bins = 200) +
  xlab("UMI counts") +
  ylab("Number of transcripts")

ggplot(umi_counts) +
  geom_histogram(aes(x = HC5), stat = "bin", bins = 200) +
  xlab("UMI counts") +
  ylab("Number of transcripts")


mean_variance_plot <- function(samples){
  mean_counts <- apply(umi_counts[, samples], 1, mean)        
  variance_counts <- apply(umi_counts[, samples], 1, var)
  df <- data.frame(mean_counts, variance_counts)
  
  ggplot(df) +
    geom_point(aes(x=mean_counts, y=variance_counts)) + 
    scale_y_log10(limits = c(1,1e9)) +
    scale_x_log10(limits = c(1,1e9)) +
    geom_abline(intercept = 0, slope = 1, color="red")  
}

mean_variance_plot(c("HB10", "HB11", "HB12", "HB13"))
mean_variance_plot(c("MET4", "MET5", "MET6"))
mean_variance_plot(c("HC20", "HC21", "HC22"))

# data <- umi_counts
# metadata_col <- "GROUP_Q1to6"
# contrast <- "MET - HC"
# comparison_num <- 2

perform_de <- function(data, metadata, metadata_col, comparison_num, contrast){
  
  file_name <- paste("de", "comp", comparison_num, sub(" - ", "_", contrast), sep = "_")
  file_name <- paste(file_name, "png", sep = ".")
  
  metadata <- metadata %>%
    column_to_rownames("SUBJECT_ORIGINAL")
  
  group <- metadata[colnames(data), metadata_col]
  
  keep <- filterByExpr(data, group=group)
  data <- data[keep, ]
  
  dge_data <- DGEList(data)
  dge_data <- calcNormFactors(dge_data)
  
  model_matrix <- model.matrix(~0 + group)
  colnames(model_matrix) <- gsub("group", "", colnames(model_matrix))
  
  png(filename = append_path("plots/de/transcriptomic", paste("voom", file_name, sep = "_")))
  y <- voom(dge_data, model_matrix, plot = TRUE)
  dev.off()
  
  # y
  vfit <- lmFit(y, model_matrix)
  # head(coef(vfit))
  
  contr_matrix <- makeContrasts(contrasts = contrast,
                                levels = colnames(model_matrix))
  vfit <- contrasts.fit(vfit, contr_matrix)
  efit <- eBayes(vfit)
  
  png(filename = append_path("plots/de/transcriptomic", paste("sa", file_name, sep = "_")))
  plotSA(efit)
  dev.off()
  
  dt <- decideTests(efit)
  summary(dt)
  
  top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
    rownames_to_column("rna")
  
  result <- top.table %>%
    select(rna, logFC, adj.P.Val) %>%
    rename(Molecule = rna, adjPVal = adj.P.Val)
  

  title <- sub("-", "Vs", contrast)
  create_volcano_plot(result, 
                      title = title, 
                      file_name = paste("volcano", file_name, sep = "_"), 
                      dir_path = "plots/de/transcriptomic", logFC_cutoff = 1)
}


perform_de(umi_counts, metadata, "GROUP_Q1to6", 2, "PREOPE - MET")
perform_de(umi_counts, metadata, "GROUP_Q1to6", 2, "PREOPE - HC")
perform_de(umi_counts, metadata, "GROUP_Q1to6", 2, "MET - HC")



umi_norm_data <- cpm(umi_counts, log = TRUE)
umi_norm_data <- scale(umi_norm_data)

umi_norm_data <- as.data.frame(t(as.matrix(umi_norm_data)))
umi_norm_data <- predict(preProcess(umi_norm_data), umi_norm_data)
umi_norm_data <- as.data.frame(t(as.matrix(umi_norm_data)))


create_box_plot_transcriptomic(umi_counts, "Raw UMI Counts", "raw_umi_boxplot.png")
create_box_plot_transcriptomic(umi_norm_data, "Norm UMI Counts", "norm_umi_boxplot.png")

#cant have -ve counts
# perform_de(umi_norm_data, metadata, "GROUP_Q1to6", 2, "PREOPE - MET")
# perform_de(umi_norm_data, metadata, "GROUP_Q1to6", 2, "PREOPE - HC")
# perform_de(umi_norm_data, metadata, "GROUP_Q1to6", 2, "MET - HC")