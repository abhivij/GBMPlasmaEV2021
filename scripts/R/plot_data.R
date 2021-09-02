library(tidyverse)
library(viridis)
library(umap)
library(ggrepel)
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


###################################

# norm_data <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_NA_equalizeMedians.csv")
# groups <- norm_data$GROUP_ORIGINAL
# norm_data <- norm_data %>%
#   select(-c(GROUP_ORIGINAL)) %>%
#   column_to_rownames("SUBJECT_ORIGINAL")
# norm_data <- data.frame(t(norm_data)) %>%
#   drop_na()
# norm_data <- data.frame(t(norm_data))
# text <- NA  

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

# result <- comparison_result_subset
# p_val_cutoff = 0.05
# logFC_cutoff = 0.5
# title = "PREOPE Vs MET"
# dir_path = "plots/de_protein/comp_2"
# file_name = "volcano_PREOPE_MET.png"
# logFC_cutoff = 0.5



create_volcano_plot <- function(result, title, file_name, dir_path = "",
                                p_val_cutoff = 0.05, logFC_cutoff = 5,
                                k = 5) {
  result$diffexpr <- "NO"
  result$diffexpr[result$logFC > logFC_cutoff & result$adjPVal < p_val_cutoff] <- "UP"
  result$diffexpr[result$logFC < -logFC_cutoff & result$adjPVal < p_val_cutoff] <- "DOWN"

  
  rownames(result) <- result$Molecule
  up_reg_result <- result %>%
    filter(logFC > 0 & adjPVal < p_val_cutoff)
  if(dim(up_reg_result)[1] == 0){
    up_reg_result <- result %>%
      filter(logFC > 0)    
  }
  up_reg_indices <- order(up_reg_result$logFC,decreasing = TRUE)[1:k]
  up_reg_indices <- rownames(up_reg_result[up_reg_indices,])
  
  down_reg_result <- result %>%
    filter(logFC < 0 & adjPVal < p_val_cutoff)
  if(dim(down_reg_result)[1] == 0){
    down_reg_result <- result %>%
      filter(logFC < 0) 
  }  
  down_reg_indices <- order(-down_reg_result$logFC, decreasing = TRUE)[1:k] 
  down_reg_indices <- rownames(down_reg_result[down_reg_indices,])
  
  result$de_mol <- NA
  result[up_reg_indices, "de_mol"] <- result[up_reg_indices, "Molecule"]
  result[down_reg_indices, "de_mol"] <- result[down_reg_indices, "Molecule"]

  
  colours <- c("red", "blue", "grey")
  names(colours) <- c("UP", "DOWN", "NO")
  
  vol_plot <- ggplot(data=result, aes(x=logFC, y=-log10(adjPVal), col=diffexpr, label=de_mol)) +
    geom_point() + 
    geom_text_repel(max.overlaps = 20, colour = "black") +
    geom_vline(xintercept=c(-logFC_cutoff, logFC_cutoff), col="green", linetype = "dashed") +
    geom_hline(yintercept=-log10(p_val_cutoff), col="green", linetype = "dashed") +
    scale_colour_manual(values = colours) +
    labs(title = title, colour = "Diff Expr")
  
  ggsave(append_path(dir_path, file_name), vol_plot) 
}
