library(tidyverse)
library(sva)
library(readxl)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

source("scripts/R/prediction_pipeline/cm_logistic_regression.R")
source("scripts/R/prediction_pipeline/cm_svm.R")
source("scripts/R/prediction_pipeline/cm_rf.R")
source("scripts/R/integration/run_all_models.R")

data_file_path_prot <- "Data/Protein/combined_data.combat.POSTOPE_TPVsREC_TP.csv"
data_file_path_tra <- "Data/RNA/combined_data.combat.POSTOPE_TPVsREC_TP.csv"
comparison = "POSTOPE_TPVsREC_TP"
conditions = c("REC_TP", "POSTOPE_TP")


base_model_pipeline <- function(comparison, conditions, 
                                data_file_path_prot,
                                data_file_path_tra){

  phenotype_file_path_prot <- "Data/proteomic_phenotype_combined.txt"
  phenotype_file_path_tra <- "Data/transcriptomic_phenotype_combined.txt"
  
  best_features_file_path <- "Data/selected_features/best_features_with_add_col.csv"
  dataset_replace_string_prot = "GBM_combined_proteomic_common_combat_"
  dataset_replace_string_tra = "GBM_combined_transcriptomic_common_combat_"
  
  best_features <- read.csv(best_features_file_path)  
  
  best_features_prot <- best_features %>%
    mutate(dataset_id = gsub(dataset_replace_string_prot, "", dataset_id)) %>%
    filter(is_best == 1, dataset_id == comparison)
  best_features_prot <- strsplit(best_features_prot$biomarkers, split = "|", fixed = TRUE)[[1]] 
  best_features_tra <- best_features %>%
    mutate(dataset_id = gsub(dataset_replace_string_tra, "", dataset_id)) %>%
    filter(is_best == 1, dataset_id == comparison)
  best_features_tra <- strsplit(best_features_tra$biomarkers, split = "|", fixed = TRUE)[[1]] 
  
  data_prot <- read.table(data_file_path_prot, header=TRUE, sep=",", row.names=1, skip=0,
                          nrows=-1, comment.char="", fill=TRUE, na.strings = "NA") 
  data_prot <- data_prot[best_features_prot, ]
  data_tra <- read.table(data_file_path_tra, header=TRUE, sep=",", row.names=1, skip=0,
                          nrows=-1, comment.char="", fill=TRUE, na.strings = "NA") 
  data_tra <- data_tra[best_features_tra, ]
  
  label_prot <- read.table(phenotype_file_path_prot, header=TRUE, sep="\t") %>%
    rename("Label" = comparison) %>%
    filter(Label %in% conditions) %>%
    dplyr::select(Sample, Label, data_cohort) %>%
    arrange(Sample)
  label_tra <- read.table(phenotype_file_path_tra, header=TRUE, sep="\t") %>%
    rename("Label" = comparison) %>%
    filter(Label %in% conditions) %>%
    dplyr::select(Sample, Label, data_cohort) %>%
    arrange(Sample)
  
  data_prot <- data_prot[, label_prot$Sample]
  label_prot$Sample <- sub("HB0", "HB", label_prot$Sample) 
  colnames(data_prot) <- label_prot$Sample
  
  data_tra <- data_tra[, label_tra$Sample]
  
  #all_labels is necessary because there are some unique samples in prot and tra
  all_labels <- rbind(label_prot, label_tra) %>%
    distinct() %>%
    arrange(Sample) %>%
    mutate(Label_data_cohort = paste(Label, data_cohort, sep = "_"))
  
  all_labels$Sample[!all_labels$Sample %in% label_prot$Sample]
  all_labels$Sample[!all_labels$Sample %in% label_tra$Sample]
  
  data_prot <- as.data.frame(t(data_prot))
  data_tra <- as.data.frame(t(data_tra))
  
  set.seed(1000)
  k_val <- 5
  times_val <- 6
  k_prod_times <- k_val * times_val
  train_index <- caret::createMultiFolds(y = all_labels$Label_data_cohort, k = k_val, times = times_val)
  
  for (i in c(1:k_prod_times)){
    # i <- 1
    all_labels.train <- all_labels[train_index[[i]], , drop = FALSE]
    all_labels.test <- all_labels[-train_index[[i]], , drop = FALSE]
    
    label_prot.train <- label_prot[label_prot$Sample %in% all_labels.train$Sample, ]
    data_prot.train <- data_prot[label_prot.train$Sample, ]
    label_prot.test <- label_prot[label_prot$Sample %in% all_labels.test$Sample, ]
    data_prot.test <- data_prot[label_prot.test$Sample, ]
    
    label_tra.train <- label_tra[label_tra$Sample %in% all_labels.train$Sample, ]
    data_tra.train <- data_tra[label_tra.train$Sample, ]
    label_tra.test <- label_tra[label_tra$Sample %in% all_labels.test$Sample, ]
    data_tra.test <- data_tra[label_tra.test$Sample, ]
    
    result_df_prot <- cbind(omics_type = "prot", 
                            run_all_models(data_prot.train, label_prot.train,
                                     data_prot.test, label_prot.test,
                                     conditions))
    result_df_tra <- cbind(omics_type = "tra", 
                           run_all_models(data_tra.train, label_tra.train,
                                          data_tra.test, label_tra.test,
                                          conditions))
    result_df <- cbind(iter = i,
                       rbind(result_df_prot,
                             result_df_tra))
    
    if(i == 1){
      result_df_all <- result_df
    } else{
      result_df_all <- rbind(result_df_all, result_df)
    }
  }
  
  result_file_dir_path <- "Data/prediction_result/integration/"
  result_file_name <- paste0(comparison, ".csv")
  if(!dir.exists(result_file_dir_path)){
    dir.create(result_file_dir_path, recursive = TRUE)
  }
  write.csv(format(result_df_all, digits = 3), paste0(result_file_dir_path, result_file_name), 
            row.names = FALSE)
}
