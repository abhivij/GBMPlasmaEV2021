library(tidyverse)
library(readxl)
library(edgeR)
library(caret)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

source("scripts/R/utils.R")
source("scripts/R/plot_data.R")

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



##############
#analysing data with censoredInt = '0'

norm_output <- read.csv(file = "Data/Protein/norm_output/norm_annotatedQ1-6_NA_equalizeMedians.csv")

sum(is.na(norm_output))

create_box_plot(norm_output, "Q1-6 Equalize median norm data censored int = 0", "q1to6_eqmed_new.png")

plot_protein_norm_data(norm_data, "q1to6_eqmed_new.png", "Q1-6 Equalize median norm data ( data censored int = 0)", dim_red = "pca")
plot_protein_norm_data(norm_data, "q1to6_eqmed_new.png", "Q1-6 Equalize median norm data ( data censored int = 0)", dim_red = "tsne")
plot_protein_norm_data(norm_data, "q1to6_eqmed_new.png", "Q1-6 Equalize median norm data ( data censored int = 0)", dim_red = "umap")

norm_data %>%
  filter(GROUP_ORIGINAL %in% c("QC1", "QC2")) %>%
  select(SUBJECT_ORIGINAL, GROUP_ORIGINAL)

norm_data %>%
  filter(SUBJECT_ORIGINAL %in% c("HB18", "HC7")) %>%
  select(SUBJECT_ORIGINAL, GROUP_ORIGINAL)
