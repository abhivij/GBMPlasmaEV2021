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
  
  log_fc_cutoff <- 0.5
  # modelBasedQCPlots(data=comparison_result, type="QQPlots",
  #                   width=5, height=5, address=append_path(dir_path, ""),
  #                   which.Protein = c(1:10))
  # # residual plots
  # modelBasedQCPlots(data=comparison_result, type="ResidualPlots",
  #                   width=5, height=5, address=append_path(dir_path, ""),
  #                   which.Protein = c(1:10))  
  
  for(l in levels(modified_comparison_result$Label)){
    print(l)
    comparison_result_subset <- modified_comparison_result %>%
      filter(Label == l) %>%
      select(Protein, log2FC, adj.pvalue) %>%
      rename(Molecule = Protein, logFC = log2FC, adjPVal = adj.pvalue)
    create_volcano_plot(comparison_result_subset, 
                        title = gsub("_", " Vs ", l), 
                        file_name = paste("volcano", paste(l, ".png", sep = ""), sep = "_"), 
                        dir_path = dir_path, logFC_cutoff = log_fc_cutoff)
  }
  
  
                       
  
}


