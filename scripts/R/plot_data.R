library(tidyverse)
library(viridis)
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
create_box_plot(norm_output1, "Q1-7 Equalize median norm data", "q1to7_eqmed.png")

norm_output2 <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_NA_FALSE.csv")
create_box_plot(norm_output2, "Q1-7 without norm data", "q1to7_nonorm.png")


file_path <- "Data/Protein/Sumofnormalisedareas/Batch1-18_modified.csv"
sum_norm_area_data <- read.table(file_path, header = TRUE, sep = ",", 
                                 comment.char = "", na.strings = "#N/A",
                                 row.names = 1)
sum_norm_area_data <- data.frame(t(log2(sum_norm_area_data))) %>%
  rownames_to_column("SUBJECT_ORIGINAL") %>%
  mutate(GROUP_ORIGINAL = "dummy")
create_box_plot(sum_norm_area_data, "Log Sum norm area data", "log_sum_norm_area.png")


plot_data <- function(norm_data, filename, title, colour_label,
                      perplexity = 30, 
                      groups, shownames = FALSE, text = NA, dim_red = "tsne"){
  
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
    n_neighbors <- max(floor(dim(norm_data)[1] / 4), 2)
    print(n_neighbors)
    result <- umap(norm_data, n_neighbors = n_neighbors)
    dim_red_df <- data.frame(x = result$layout[,1], y = result$layout[,2], 
                             Colour = groups, 
                             Sample = text)  
    xlab <- "UMAP 1"
    ylab <- "UMAP 2"
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
  filename <- paste("plots", filename, sep = "/")
  ggplot2::ggsave(filename, dim_red_plot)  
}

norm_data <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_NA_equalizeMedians.csv")
plot_protein_norm_data <- function(norm_data, dirname = "", filename, 
                                   title, width = 30, height = 30, perplexity = 30){
  groups <- norm_data$GROUP_ORIGINAL
  norm_data <- norm_data %>%
    select(-c(GROUP_ORIGINAL)) %>%
    column_to_rownames("SUBJECT_ORIGINAL")
  norm_data <- t(norm_data)
}



create_plot <- function(norm_data, dirname = "", filename, title, width = 30, height = 30, perplexity = 30){
  
  ne <- t(norm_data)
  
  # PCA on all samples and all proteins
  all_group_names <- ne["GROUP_ORIGINAL", ]
  all_sample_names <- ne["SUBJECT_ORIGINAL", ]
  colnames(ne) <- ne["SUBJECT_ORIGINAL", ]
  ne <- ne[-which(row.names(ne) %in% c("GROUP_ORIGINAL", "SUBJECT_ORIGINAL")), ]
  ne <- data.frame(ne)
  ne <- sapply(ne, function(x) as.numeric(x))  #causes rownames that is protein names to be lost
  ne <- t(ne)
  
  pca_file_name <- append_path(dirname, paste("pca", filename, sep = "_"))
  p <- prcomp(ne)
  pca_plotdata <- data.frame(p$x) %>%
    mutate(type = all_group_names) %>%
    mutate(name = all_sample_names)
  
  tsne_file_name <- append_path(dirname, paste("tsne", filename, sep = "_"))
  set.seed(1)
  tsne_result <- Rtsne::Rtsne(ne, perplexity = perplexity)
  tsne_df <- data.frame(x = tsne_result$Y[,1], y = tsne_result$Y[,2], 
                        Colour = all_group_names,
                        name = all_sample_names)
  
  if (dim(ne)[1] < 20) {
    pca_plotdata %>%
      ggplot(aes(x = PC1, y = PC2, colour = type, label = name)) +
      geom_point(size = 2) + 
      labs(title = title) +
      ggrepel::geom_text_repel(show.legend = FALSE)
    ggsave(pca_file_name, width = width, height = height, units = "cm")  
    
    tsne_plot <- tsne_df %>%
      ggplot(aes(x = x, y = y, colour = Colour, label = name)) +
      geom_point(size = 2) +
      ggplot2::labs(colour = "Groups", title = title) +
      ggplot2::xlab("tSNE 1") +
      ggplot2::ylab("tSNE 2") +
      ggrepel::geom_text_repel(show.legend = FALSE)  
  } else {
    pca_plotdata %>%
      ggplot(aes(x = PC1, y = PC2, colour = type, label = name)) +
      geom_point(size = 2) + 
      labs(title = title)
    ggsave(pca_file_name, width = width, height = height, units = "cm")  
    
    tsne_plot <- tsne_df %>%
      ggplot(aes(x = x, y = y, colour = Colour, label = name)) +
      geom_point(size = 2) +
      ggplot2::labs(colour = "Groups", title = title) +
      ggplot2::xlab("tSNE 1") +
      ggplot2::ylab("tSNE 2")
  }
  
  ggplot2::ggsave(tsne_file_name, tsne_plot)
  
}

data <- read.csv("Data/Protein/sample_columns.csv")

sample_column <- rbind(data.frame(Sample = data$Column1, Column = 'col1'),
                       data.frame(Sample = data$Column2, Column = 'col2'),
                       data.frame(Sample = data$Column3, Column = 'col3')) %>%
  filter(Sample != "")
sample_column_qc <- sample_column %>%
  filter(grepl("HB18|HC7", Sample))

write.csv(sample_column_qc, "Data/Protein/sample_column_qc.csv", row.names = FALSE, quote = FALSE)


norm_data <- read.csv(file = "Data/normQ1-6.csv")

sum(is.na(norm_data))
for(i in c(1:dim(norm_data)[1])){
  for(j in c(1:dim(norm_data)[2])){
    if(is.na(norm_data[i,j])){
      print(paste(i,j))
      break
    }
  }
}

#filtering out column with missing value
norm_data <- norm_data %>%
  select(-c(185))


create_plot(norm_data, dirname = "plots/all_q1-6", filename = "all_proteins.jpg", title = "All samples with all proteins")



norm_data_qc <- norm_data %>%
  filter(grepl("HB18|HC7", SUBJECT_ORIGINAL)) %>%
  inner_join(sample_column_qc, by = c("SUBJECT_ORIGINAL" = "Sample")) %>%
  select(-GROUP_ORIGINAL) %>%
  rename(GROUP_ORIGINAL = Column)

create_plot(norm_data = norm_data_qc, dirname = "plots/all_q1-6", filename = "qc_proteins.jpg", title = "QC samples with all proteins", 
            height = 10, width = 10, perplexity = 2)


column_data <- norm_data %>%
  inner_join(sample_column, by = c("SUBJECT_ORIGINAL" = "Sample")) %>%
  select(-GROUP_ORIGINAL) %>%
  rename(GROUP_ORIGINAL = Column)

create_plot(norm_data = column_data, dirname = "plots/all_q1-6", filename = "all_proteins_by_col.jpg", title = "All samples with all proteins based on columns", 
            height = 10, width = 10)
