library(MSstats)
library(tidyverse)
# library(stringr)   for str_detect()




args = commandArgs(trailingOnly = TRUE)

if(length(args) >= 2){
  data_dir <- args[2]
} else{
  data_dir <- "annotatedQ1-6"
}
if(length(args) >= 3){
  condition_type <- args[3]
} else {
  condition_type <- NA
}
if(length(args) >= 4){
  norm <- args[4]
} else {
  norm <- "equalizeMedians"
}



setwd("~/UNSW/VafaeeLab/GBMPlasmaEV")
# setwd("/srv/scratch/vafaeelab/AbhishekVijayan/GBMPlasmaEV")


print(paste("Data dir :", data_dir))
data_path <- paste("Data/Protein/MSstatsinput", data_dir, sep = "")
file_names <- list.files(path = data_path, full.names = TRUE)

# a <- unique(sample_column$Sample)  
# b <- unique(data$BioReplicate)
# a[!a %in% b]
# b[!b %in% a]

modify_condition <- function(data, condition_type = NA, sample_column = NA) {
  allowed_condition_types <- c("column", "disease")
  if (condition_type %in% allowed_condition_types) {
    if(condition_type == "column"){
      if(is.na(sample_column)){
        print("specify value for sample_column")
      }
      else{
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

get_data_from_file <- function(filename, condition_type = NA, sample_column = NA){
  data <- read.csv(filename)
  return (modify_condition(data, condition_type, sample_column))
}

file_names <- file_names[1:3]
print(file_names)

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
data_process_output <- dataProcess(data, logTrans = "2", normalization = norm)
normed <- data_process_output$RunlevelData %>%
  select(Protein, LogIntensities, GROUP_ORIGINAL, SUBJECT_ORIGINAL)

prots <- strsplit(as.character(normed$Protein), split="\\|")
prots <- unlist(lapply(prots, function(x) x[2]))
normed$Protein <- prots
normed <- normed %>%
  pivot_wider(names_from = Protein, values_from = LogIntensities)

filename <- paste(paste("norm", data_dir, condition_type, norm, sep = "_"), "csv", sep = ".")
print(filename)
write.table(normed, filename, quote = FALSE, sep = ",", row.names = FALSE)