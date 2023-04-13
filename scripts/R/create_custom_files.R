library(tidyverse)
library(xlsx)


#write best biomarkers - combined + combat into excel file

best_features <- read.csv("Data/selected_features/best_features_with_add_col.csv") %>%
  filter(is_best == 1) %>%
  filter(grepl(pattern = "GBM_combined_", x = dataset_id)) %>%
  mutate(dataset_id = sub("GBM_combined_", "", dataset_id, fixed = TRUE)) %>%
  mutate(dataset_id = sub("_common_combat", "", dataset_id, fixed = TRUE))
all_protein_names <- read.csv("Data/Protein/formatted_data/all_protein_names.csv")

file_name <- "combined_common_combat_best_features.xlsx"


for(i in c(1:nrow(best_features))){
  dataset_id <- best_features[i, "dataset_id"]
  biomarkers <- strsplit(best_features[i, "biomarkers"], split = "|", fixed = TRUE)[[1]]
  
  biomarkers_df <- data.frame(biomarkers = biomarkers)
  if(grepl("transcriptomic", dataset_id)){
    #replace _ in mirna names to - : that's how it is represnted in mirbase
    #but pirnabank has _ in names
    biomarkers_df <- biomarkers_df %>%
      mutate(biomarkers = ifelse(!grepl("piR", biomarkers), 
                                 gsub("_", "-", biomarkers, fixed = TRUE),
                                 biomarkers))
  } else{
    #this case is proteomics
    #add protein name
    biomarkers_df <- biomarkers_df %>%
      left_join(all_protein_names, by = c("biomarkers" = "from_id")) %>%
      relocate(primary_gene_id, .before = protein_name)
  }
  
  #print(biomarkers_df)
  # 
  print(dim(biomarkers_df))
  write.xlsx(biomarkers_df, file = paste0("Data/selected_features/", file_name), 
             append = TRUE, sheetName = dataset_id, row.names = FALSE)
}
