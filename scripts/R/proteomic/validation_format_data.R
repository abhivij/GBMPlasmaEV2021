library(tidyverse)

data_path <- "Data/Protein/MSstatsinput_newcohort_Oct28"
validation_metadata <- read.csv("Data/metadata_glionet.csv", stringsAsFactors = FALSE)
output_dir_path <- "Data/Protein/MSstatsinput_newcohort_processed"

file_names <- list.files(path = data_path, full.names = TRUE)
for(f in file_names){
  
  # f <- file_names[1]
  data <- read.csv(f)
  
  # data_copy <- data
  # data <- data_copy
  # length(unique(data$Replicate))
  
  data <- data %>%
    mutate(BioReplicate = sub("20220714_Eclipse-OT_SH_|20220719_Eclipse-OT_SH_", "", File.Name)) %>%
    mutate(BioReplicate = sub(".raw", "", BioReplicate, fixed = TRUE)) %>%
    mutate(BioReplicate = sub("-", "", BioReplicate, fixed = TRUE))

  data <- data %>%
    left_join(validation_metadata %>% select(sample_id, category_old_name), by = c("BioReplicate" = "sample_id")) %>%
    mutate(Condition = case_when(is.na(category_old_name) ~ 'QC',
                                 TRUE ~ category_old_name), .before = BioReplicate) %>%
    select(-c(category_old_name))
  
  # length(unique(data$BioReplicate))

  f_actual_name <- sub(paste0(data_path,"/"), "", f, fixed = TRUE)
  if(!dir.exists(output_dir_path)){
    dir.create(output_dir_path, recursive = TRUE)
  }
  write.csv(data, paste(output_dir_path, f_actual_name, sep = "/"), row.names = FALSE)
}
