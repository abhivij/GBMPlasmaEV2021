library(MSstats)
library(tidyverse)
source("scripts/R/utils.R")
setwd("../..")

plot_de_results <- function(file_name){
  print(file_name)
  print(sapply(sys.call()[-1], deparse))
  file_name <- paste("comparison_result", comparison_num,
                     gsub(pattern = "data_process_output_", replacement = "", file_name),
                     sep = "_")
  comparison_result <- readRDS(append_path("Data/Protein/comparison_result", file_name))
}

test <- function(a, ...){
  print(a)
  print("others")
  print(sapply(sys.call()[-1], deparse))
}

