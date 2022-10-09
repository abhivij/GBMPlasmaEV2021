library(MSstats)
library(tidyverse)
source("utils.R")
setwd("../..")

modify_condition <- function(data, condition_type = NA, sample_column = NA) {
  allowed_condition_types <- c("column", "disease")
  if (condition_type %in% allowed_condition_types) {
    if(condition_type == "column"){
      if(is.na(sample_column)){
        print("specify value for sample_column")
      }
      else{
        a <- unique(sample_column$Sample)
        b <- unique(data$BioReplicate)
        print(a[!a %in% b])
        print(b[!b %in% a])
        
        data <- data %>%
          inner_join(sample_column, by = c("BioReplicate" = "Sample")) %>%
          mutate(Condition = Column) %>%
          select(-c(Column))
      }
    }
    else if(condition_type == "disease") {
      data <- data %>%
        mutate(Condition = ifelse(grepl("HB", BioReplicate), 'GBM', 
                                  ifelse(grepl("MET", BioReplicate), 'MET', 'HC')))  
    }
  }
  return (data)
}

get_data_from_file <- function(file_name, condition_type = NA, sample_column = NA){
  data <- read.csv(file_name)
  print(paste(file_name, condition_type))
  return (modify_condition(data, condition_type, sample_column))
}



process_protein_data <- function(data_dir, condition_type, norm, remove_50_missing = FALSE){
  print(paste("Data dir :", data_dir))
  data_path <- paste("Data/Protein/MSstatsinput", data_dir, sep = "")
  file_names <- list.files(path = data_path, full.names = TRUE)
  
  if(!is.na(condition_type) && condition_type == "column"){
    sample_column <- read.csv("Data/Protein/sample_columns.csv", na.strings = "") %>%
      pivot_longer(cols = everything(), names_to = "Column", values_to = "Sample", values_drop_na = TRUE)
    data <- do.call(rbind,
                    lapply(file_names, get_data_from_file, "column", sample_column))  
  } else {
    data <- do.call(rbind,
                    lapply(file_names, get_data_from_file, condition_type))
  }
  
  
  print(paste("Normalization :", norm))
  
  data <- SkylinetoMSstatsFormat(data)
  data_process_output <- dataProcess(data, logTrans = 2, normalization = norm,
                                     censoredInt = '0', remove50missing = remove_50_missing)
  
  file_name_substring <- condition_type
  if(remove_50_missing){
    file_name_substring <- paste(file_name_substring, "remove_50_missing", sep = "_")
  }
  
  
  file_name <- paste(paste("data_process_output", data_dir, file_name_substring, norm, sep = "_"), "rds", sep = ".")
  output_dir <- "Data/Protein/data_process_output"
  saveRDS(data_process_output, file = append_path(output_dir, file_name))
  
  normed <- data_process_output$RunlevelData %>%
    select(ProteinName, LogIntensities, GROUP_ORIGINAL, SUBJECT_ORIGINAL) %>%
    separate(ProteinName, c(NA, "Protein", NA), sep = "\\|") %>% 
    pivot_wider(names_from = Protein, values_from = LogIntensities)
  
  file_name <- paste(paste("norm", data_dir, file_name_substring, norm, sep = "_"), "csv", sep = ".")
  output_dir <- "Data/Protein/norm_output"
  write.table(normed, append_path(output_dir, file_name), 
              quote = FALSE, sep = ",", row.names = FALSE)
  
}








