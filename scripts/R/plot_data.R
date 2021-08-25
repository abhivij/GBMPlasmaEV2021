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