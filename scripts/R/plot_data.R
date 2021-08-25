library(tidyverse)
library(viridis)
library(umap)
setwd("/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV/scripts/R")
source("utils.R")

setwd("/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV")

  
create_box_plot <- function(data, title, filename){
  data <- data %>%
    select(-c(GROUP_ORIGINAL)) %>%
    pivot_longer(!SUBJECT_ORIGINAL, names_to = "Protein", values_to = "LogIntensity")
  box_plot <- data %>%
    ggplot(aes(x=SUBJECT_ORIGINAL, y=LogIntensity)) +
    geom_boxplot() +
    ggtitle(title) +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.99, size=rel(0.5)))
  filename <- paste("plots/box_plot", filename, sep = "/")
  ggsave(filename, box_plot)  
}

norm_output1 <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_NA_equalizeMedians.csv")
create_box_plot(norm_output1, "Q1-6 Equalize median norm data", "q1to6_eqmed.png")

norm_output2 <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_NA_FALSE.csv")
create_box_plot(norm_output2, "Q1-6 without norm data", "q1to6_nonorm.png")


file_path <- "Data/Protein/Sumofnormalisedareas/Batch1-18_modified.csv"
sum_norm_area_data <- read.table(file_path, header = TRUE, sep = ",", 
                                 comment.char = "", na.strings = "#N/A",
                                 row.names = 1)
sum_norm_area_data <- data.frame(t(log2(sum_norm_area_data))) %>%
  rownames_to_column("SUBJECT_ORIGINAL") %>%
  mutate(GROUP_ORIGINAL = "dummy")
create_box_plot(sum_norm_area_data, "Log Sum norm area data", "log_sum_norm_area.png")

###################################

norm_data <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_NA_equalizeMedians.csv")
groups <- norm_data$GROUP_ORIGINAL
norm_data <- norm_data %>%
  select(-c(GROUP_ORIGINAL)) %>%
  column_to_rownames("SUBJECT_ORIGINAL")
norm_data <- data.frame(t(norm_data)) %>%
  drop_na()
norm_data <- data.frame(t(norm_data))
text <- NA  

plot_data <- function(norm_data, filename, title, colour_label,
                      perplexity = 30, 
                      groups, shownames = FALSE, text = NA, dim_red = "pca"){
  
  if(shownames && length(text) == 1 && is.na(text)){
    text <- rownames(norm_data)
  }
  set.seed(1)
  
  if(dim_red == "tsne"){
    result <- Rtsne::Rtsne(norm_data, perplexity = perplexity)
    dim_red_df <- data.frame(x = result$Y[,1], y = result$Y[,2], 
                             Colour = groups, 
                             Sample = text)    
    xlab <- "tSNE 1"
    ylab <- "tSNE 2"
  } else if(dim_red == "umap"){
    print(dim(norm_data)[1])
    # n_neighbors <- max(floor(dim(norm_data)[1] / 4), 2)
    # print(n_neighbors)
    # result <- umap(norm_data, n_neighbors = n_neighbors)
    result <- umap(norm_data)
    dim_red_df <- data.frame(x = result$layout[,1], y = result$layout[,2], 
                             Colour = groups, 
                             Sample = text)  
    xlab <- "UMAP 1"
    ylab <- "UMAP 2"
  } else if(dim_red == "pca"){
    result <- prcomp(norm_data)
    dim_red_df <- data.frame(x = result$x[,1], y = result$x[,2], 
                             Colour = groups, 
                             Sample = text)    
    xlab <- "PCA 1"
    ylab <- "PCA 2"    
  }
  
  
  if (shownames) {
    
    dim_red_plot <- ggplot2::ggplot(dim_red_df,
                                    ggplot2::aes(x = x, y = y, colour = Colour)) +
      ggplot2::geom_point() +
      geom_text_repel(aes(label=Sample)) +
      ggplot2::labs(title = title, colour = colour_label) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab)    
  } else {
    dim_red_plot <- ggplot2::ggplot(dim_red_df) +
      ggplot2::geom_point(ggplot2::aes(x = x, y = y, colour = Colour)) +
      ggplot2::labs(title = title, colour = colour_label) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab)
  }
  
  filename <- paste(dim_red, filename, sep = "_")
  filename <- paste("plots/dim_red", filename, sep = "/")
  ggplot2::ggsave(filename, dim_red_plot)  
}

plot_protein_norm_data <- function(norm_data, filename, 
                                   title, width = 30, height = 30, perplexity = 30,
                                   dim_red = "pca", colour_label = "Labels"){
  groups <- norm_data$GROUP_ORIGINAL
  norm_data <- norm_data %>%
    select(-c(GROUP_ORIGINAL)) %>%
    column_to_rownames("SUBJECT_ORIGINAL")
  norm_data <- data.frame(t(norm_data)) %>%
    drop_na()
  norm_data <- data.frame(t(norm_data))
  
  plot_data(norm_data, filename, title, colour_label = colour_label,
                        perplexity = 30, 
                        groups, shownames = FALSE, text = NA, dim_red)
}

norm_data <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_NA_equalizeMedians.csv")
plot_protein_norm_data(norm_data, "q1to6_eqmed.png", "Q1-6 Equalize median norm data", dim_red = "pca")
plot_protein_norm_data(norm_data, "q1to6_eqmed.png", "Q1-6 Equalize median norm data", dim_red = "tsne")
plot_protein_norm_data(norm_data, "q1to6_eqmed.png", "Q1-6 Equalize median norm data", dim_red = "umap")

norm_data <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_column_equalizeMedians.csv")
plot_protein_norm_data(norm_data, "column_q1to6_eqmed.png", "Column Q1-6 Equalize median norm data", dim_red = "pca")
plot_protein_norm_data(norm_data, "column_q1to6_eqmed.png", "Column Q1-6 Equalize median norm data", dim_red = "tsne")
plot_protein_norm_data(norm_data, "column_q1to6_eqmed.png", "Column Q1-6 Equalize median norm data", dim_red = "umap")


norm_data <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_column_FALSE.csv")
plot_protein_norm_data(norm_data, "column_q1to6_nonorm.png", "Column Q1-6 no norm data", dim_red = "pca")
plot_protein_norm_data(norm_data, "column_q1to6_nonorm.png", "Column Q1-6 no norm data", dim_red = "tsne")
plot_protein_norm_data(norm_data, "column_q1to6_nonorm.png", "Column Q1-6 no norm data", dim_red = "umap")

norm_data <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_disease_equalizeMedians.csv")
plot_protein_norm_data(norm_data, "disease_q1to6_eqmed.png", "Disease Q1-6 Equalize median norm data", dim_red = "pca")
plot_protein_norm_data(norm_data, "disease_q1to6_eqmed.png", "Disease Q1-6 Equalize median norm data", dim_red = "tsne")
plot_protein_norm_data(norm_data, "disease_q1to6_eqmed.png", "Disease Q1-6 Equalize median norm data", dim_red = "umap")

norm_data <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_disease_FALSE.csv")
plot_protein_norm_data(norm_data, "disease_q1to6_nonorm.png", "Disease Q1-6 no norm data", dim_red = "pca")
plot_protein_norm_data(norm_data, "disease_q1to6_nonorm.png", "Disease Q1-6 no norm data", dim_red = "tsne")
plot_protein_norm_data(norm_data, "disease_q1to6_nonorm.png", "Disease Q1-6 no norm data", dim_red = "umap")



norm_data1 <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_NA_equalizeMedians.csv") %>%
  select(c(SUBJECT_ORIGINAL, GROUP_ORIGINAL)) %>%
  rename("Q1to6Condition" = "GROUP_ORIGINAL")
norm_data5 <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_disease_equalizeMedians.csv") %>%
  select(c(SUBJECT_ORIGINAL, GROUP_ORIGINAL)) %>%
  rename("DiseaseType" = "GROUP_ORIGINAL")
norm_data3 <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_column_equalizeMedians.csv") %>%
  select(c(SUBJECT_ORIGINAL, GROUP_ORIGINAL)) %>%
  rename("Column" = "GROUP_ORIGINAL")

meta_data <- norm_data1 %>%
  inner_join(norm_data5) %>%
  full_join(norm_data3)


  
plot_protein_sum_of_area_data <- function(norm_data, meta_data, condition,
                                          l = FALSE, filename, 
                                          title, width = 30, height = 30, perplexity = 30,
                                          dim_red = "pca", colour_label = "Labels"){
  
  if(l == TRUE){
      norm_data <- log2(norm_data)
  }
  meta_data <- meta_data %>%
    select(c("SUBJECT_ORIGINAL", condition))
  norm_data <- data.frame(t(norm_data)) %>%
    rownames_to_column("SUBJECT_ORIGINAL") %>%
    inner_join(meta_data) %>%
    rename("GROUP_ORIGINAL" = condition)
  plot_protein_norm_data(norm_data, filename, 
                                     title, width = 30, height = 30, perplexity = 30,
                                     dim_red = dim_red, colour_label = colour_label)
}


file_path <- "Data/Protein/Sumofnormalisedareas/Batch1-18_modified.csv"
sum_norm_area_data <- read.table(file_path, header = TRUE, sep = ",", 
                                 comment.char = "", na.strings = "#N/A",
                                 row.names = 1)


plot_protein_sum_of_area_data(sum_norm_area_data, meta_data, "DiseaseType",
                              filename = "sumnormarea_disease.png", title = "Disease from Sum norm Area data", 
                              dim_red = "pca")
plot_protein_sum_of_area_data(sum_norm_area_data, meta_data, "DiseaseType",
                              filename = "sumnormarea_disease.png", title = "Disease from Sum norm Area data", 
                              dim_red = "umap")
plot_protein_sum_of_area_data(sum_norm_area_data, meta_data, "DiseaseType",
                              filename = "log_sumnormarea_disease.png", title = "Disease from Log Sum norm Area data", 
                              dim_red = "pca", l = TRUE)
plot_protein_sum_of_area_data(sum_norm_area_data, meta_data, "DiseaseType",
                              filename = "log_sumnormarea_disease.png", title = "Disease from Log Sum norm Area data", 
                              dim_red = "umap", l = TRUE)

plot_protein_sum_of_area_data(sum_norm_area_data, meta_data, "Column",
                              filename = "sumnormarea_column.png", title = "Column from Sum norm Area data", 
                              dim_red = "pca")
plot_protein_sum_of_area_data(sum_norm_area_data, meta_data, "Column",
                              filename = "sumnormarea_column.png", title = "Column from Sum norm Area data", 
                              dim_red = "umap")
plot_protein_sum_of_area_data(sum_norm_area_data, meta_data, "Column",
                              filename = "log_sumnormarea_column.png", title = "Column from Log Sum norm Area data", 
                              dim_red = "pca", l = TRUE)
plot_protein_sum_of_area_data(sum_norm_area_data, meta_data, "Column",
                              filename = "log_sumnormarea_column.png", title = "Column from Log Sum norm Area data", 
                              dim_red = "umap", l = TRUE)
