library(MSstats)
library(tidyverse)
source("scripts/R/utils.R")
setwd("../..")

file_name <- "data_process_output_annotatedQ1-6_NA_equalizeMedians.rds"
comparison_num <- "2"
plot_de_results <- function(file_name){
  print(file_name)
  print(sapply(sys.call()[-1], deparse))
  file_name <- paste("comparison_result", comparison_num,
                     gsub(pattern = "data_process_output_", replacement = "", file_name),
                     sep = "_")
  comparison_result <- readRDS(append_path("Data/Protein/comparison_result", file_name))
  names(comparison_result)
  head(comparison_result$ComparisonResult)
  levels(comparison_result$ComparisonResult$Label)
  
  dim(comparison_result$ComparisonResult)
  
  dim(comparison_result$ComparisonResult %>% filter(Label == "PREOPE_MET"))
  dim(comparison_result$ComparisonResult %>% filter(Label == "PREOPE_HC"))
  dim(comparison_result$ComparisonResult %>% filter(Label == "MET_HC")) 
  
  modified_comparison_result <- comparison_result$ComparisonResult %>%
    separate(Protein, c(NA, NA, "Protein"), sep = "\\|") %>%
    mutate(Protein = gsub(pattern = "_HUMAN", replacement = "", Protein))  
  
  dir_path <- paste("plots/de_protein/comp", comparison_num, sep = "_")
  if(!dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE)
  }
  
  modelBasedQCPlots(data=comparison_result, type="QQPlots",
                    width=5, height=5, address=append_path(dir_path, ""),
                    which.Protein = c(1:10))
  # residual plots
  modelBasedQCPlots(data=comparison_result, type="ResidualPlots",
                    width=5, height=5, address=append_path(dir_path, ""),
                    which.Protein = c(1:10))  
  
  groupComparisonPlots(data = modified_comparison_result, type = 'VolcanoPlot',
                       address=append_path(dir_path, ""))
  
  groupComparisonPlots(data = modified_comparison_result, type = 'Heatmap',
                       address = append_path(dir_path, ""))
  
}


