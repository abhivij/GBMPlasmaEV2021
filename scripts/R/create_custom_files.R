library(tidyverse)
library(xlsx)

create_biomarker_file <- function(best_features, file_name, dir_path = "Data/selected_features/"){
  for(i in c(1:nrow(best_features))){
    dataset_id <- paste(best_features[i, "dataset_id"], i, sep = "_")
    is_best <- ifelse(best_features[i, "is_best"] == 1, "best", "")
    details <- paste(best_features[i, "description"], best_features[i, "min_iter_feature_presence"], sep = "_")
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

    biomarkers_df <- biomarkers_df %>%
      mutate(method = details, .before = biomarkers)
    
    print(dim(biomarkers_df))
    dataset_id <- sub("transcriptomic", "tr_", dataset_id)
    dataset_id <- sub("proteomic", "pr_", dataset_id)
    sheet_name <- paste(dataset_id, is_best, sep = "_")
    print(sheet_name)
    write.xlsx(biomarkers_df, file = paste0(dir_path, file_name), 
               append = TRUE, sheetName = sheet_name, 
               row.names = FALSE)
  }
}

best_features <- read.csv("Data/selected_features/best_features_with_add_col.csv") %>%
  filter(is_best == 1) %>%
  filter(grepl(pattern = "GBM_combined_", x = dataset_id)) %>%
  filter(grepl(pattern = "_common_combat", x = dataset_id)) %>%
  filter(!grepl(pattern = "_common_combat_mod", x = dataset_id)) %>%
  mutate(dataset_id = sub("GBM_combined_", "", dataset_id, fixed = TRUE)) %>%
  mutate(dataset_id = sub("_common_combat", "", dataset_id, fixed = TRUE))
all_protein_names <- read.csv("Data/Protein/formatted_data/all_protein_names.csv")

create_biomarker_file(best_features, file_name = "combined_common_combat_best_features_withnewsheetname.xlsx")





best_features <- read.csv("Data/selected_features/best_features_with_add_col.csv") %>%
  filter(is_best > 0) %>%
  filter(grepl(pattern = "GBM_combined_", x = dataset_id)) %>%
  filter(grepl(pattern = "_combat_compset2_", x = dataset_id)) %>%
  mutate(dataset_id = sub("GBM_combined_", "", dataset_id, fixed = TRUE)) %>%
  mutate(dataset_id = sub("_combat_compset2_", "", dataset_id, fixed = TRUE))
all_protein_names <- read.csv("Data/Protein/formatted_data/all_protein_names.csv")

create_biomarker_file(best_features, file_name = "PREOPE_MET_HC_best_features.xlsx")
