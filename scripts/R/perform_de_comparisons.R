library(MSstats)
library(tidyverse)
source("utils.R")
setwd("../..")

#creates comparison object as shown below
# > levels(data_process_output$RunlevelData$GROUP_ORIGINAL)
# [1] "HC"        "MET"       "OUT"       "POSTOPE-P" "POSTOPE-T" "PREOPE"    "PREREC"    "QC1"      
# [9] "QC2"       "REC-P"     "REC-T" 

# comparison <- rbind(c(0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0), # PREOPE Vs MET
#                     c(-1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), # PREOPE Vs HC
#                     c(-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)) # MET Vs HC
# row.names(comparison) <- c("PREOPE_MET", 
#                            "PREOPE_HC",
#                            "MET_HC")

perform_de_comparisons <- function(file_name, comparison_list, comparison_num){
  file_path <- append_path("Data/Protein/data_process_output", file_name)
  
  data_process_output <- readRDS(file_path)
  
  all_groups <- levels(data_process_output$RunlevelData$GROUP_ORIGINAL)

  for(i in c(1:length(comparison_list))){
    comparison_row <- rep(0, length(all_groups))
    comparison_row[match(comparison_list[[i]][1], all_groups)] <- 1
    comparison_row[match(comparison_list[[i]][2], all_groups)] <- -1  
    if(i == 1){
      comparison <- rbind(comparison_row)
    } else {
      comparison <- rbind(comparison, comparison_row)
    }
  }
  row.names(comparison) <- sapply(comparison_list, function(x){
    return(paste(x[1], x[2], sep = "_"))
  })  
  
  comparison_result <- groupComparison(contrast.matrix = comparison,
                                       data = data_process_output)

  file_name <- paste("comparison_result", comparison_num,
        gsub(pattern = "data_process_output_", replacement = "", file_name),
        sep = "_")
  output_dir <- "Data/Protein/comparison_result"
  saveRDS(comparison_result, file = append_path(output_dir, file_name))
  
}

