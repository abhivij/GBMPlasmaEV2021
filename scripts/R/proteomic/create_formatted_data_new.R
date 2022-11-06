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