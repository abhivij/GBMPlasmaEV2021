library(tidyverse)
library(readxl)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

if(!dir.exists("Data/Protein/formatted_data/")){
  dir.create("Data/Protein/formatted_data/")
}

input_file_path <- "Data/Protein/norm_output/norm_annotatedQ1-6_NA_FALSE.csv"
# max(formatted_data, na.rm = TRUE)
# 25.10626
# min(formatted_data, na.rm = TRUE)
# -2.664906

# input_file_path <- "Data/Protein/norm_output/norm_annotatedQ1-6_NA_equalizeMedians.csv"
# max(formatted_data, na.rm = TRUE)
# 33.42296
# min(formatted_data, na.rm = TRUE)
# -4.876739

# input_file_path <- "Data/Protein/norm_output/norm_annotatedQ1-6_NA_quantile.csv"
# max(formatted_data, na.rm = TRUE)
# 23.75979
# min(formatted_data, na.rm = TRUE)
# -8.416577

output_file_path <- "Data/Protein/formatted_data/Q1-6_nonorm_formatted.csv"
process_and_format_protein_data <- function(input_file_path, output_file_path){
  protein_data <- read.csv(file = input_file_path)
  
  formatted_data <-  protein_data %>%
    select(-c(GROUP_ORIGINAL)) %>%
    column_to_rownames("SUBJECT_ORIGINAL")
  
  max(formatted_data, na.rm = TRUE)
  min(formatted_data, na.rm = TRUE)
  quantiles <- quantile(formatted_data, na.rm = TRUE)
  
  na_repl_value <- quantiles["0%"] - quantiles["25%"]
  
  sum(is.na(formatted_data))
  formatted_data[is.na(formatted_data)] <- na_repl_value
  sum(is.na(formatted_data))
  
  formatted_data <- t(formatted_data)
  write.csv(formatted_data, output_file_path)
}

process_and_format_protein_data("Data/Protein/norm_output/norm_annotatedQ1-6_NA_FALSE.csv",
                                "Data/Protein/formatted_data/Q1-6_nonorm_formatted.csv")
process_and_format_protein_data("Data/Protein/norm_output/norm_annotatedQ7_NA_FALSE.csv",
                                "Data/Protein/formatted_data/Q7_nonorm_formatted.csv")

# data <- read.table(output_file_path, header=TRUE, sep=",", row.names=1, skip=0,
#                    nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")

