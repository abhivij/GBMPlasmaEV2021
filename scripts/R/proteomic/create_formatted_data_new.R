library(tidyverse)
library(readxl)

# base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
# setwd(base_dir)

source("scripts/R/utils.R")

if(!dir.exists("Data/Protein/formatted_data/")){
  dir.create("Data/Protein/formatted_data/", recursive = TRUE)
}

process_and_format_protein_data("Data/Protein/norm_output/norm__newcohort_processed_NA_FALSE.csv",
                                "Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil.csv",
                                impute = TRUE, filter_na_perc = 50)

##################################################################
### create data files with common proteins 

data <- read.csv("Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv", 
                 row.names=1)
validation_data <- read.csv("Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil.csv", 
                            row.names = 1)
length(intersect(rownames(data), rownames(validation_data)))
common_proteins <- intersect(rownames(data), rownames(validation_data))

initial_data_common <- data[common_proteins, ]
validation_data_common <- validation_data[common_proteins, ]

write.csv(initial_data_common, 
          "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil_common.csv")
write.csv(validation_data_common, 
          "Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil_common.csv")

##################################################################
