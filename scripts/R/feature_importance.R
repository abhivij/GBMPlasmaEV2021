#used to obtain feature importance of transcripts/proteins of interest, by passing them through Random Forest classifier
#particularly created to get feature importances of DE transcripts

library(tidyverse)
library(sva)
library(readxl)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

source("scripts/R/prediction_pipeline/cm_rf.R")



#returns feature importance of all features in the data
#so ideal to pass as argument only the data subset with features of interest

# data_file_path = "Data/RNA_all/newquant_Nov2023_umi_counts_PREOPE_MET_HC_filter90_de_PREOPEVsMET.csv"
# phenotype_file_path = "Data/transcriptomic_phenotype_PREOPE_MET_HC.txt"
#
# comparison = "PREOPEVsMET"
# classes = c("MET", "PREOPE")
# perform_filter = TRUE
# norm = "log_cpm"
# 
# logfc_file_path <- "plots_RNA_all/PREOPE_MET_HC/volcano/volcano_1_PREOPEVsMET.csv"
# result_file_dir_path <- "plots_RNA_all/PREOPE_MET_HC/volcano/"
# 
# file_name_prefix = "DE_"


rf_feature_importance <- function(data_file_path,
                                  phenotype_file_path,
                                  logfc_file_path,
                                  result_file_dir_path,
                                  comparison, classes,
                                  perform_filter = FALSE, norm = "none", file_name_prefix = ""){

  data <- read.table(data_file_path, header=TRUE, sep=",", row.names=1, skip=0,
                     nrows=-1, comment.char="", fill=TRUE, na.strings = "NA") 
  
  label <- read.table(phenotype_file_path, header=TRUE, sep="\t") %>%
    dplyr::rename("Label" = comparison) %>%
    dplyr::filter(Label %in% classes) %>%
    dplyr::select(Sample, Label, data_cohort) %>%
    arrange(Sample) %>%
    mutate(Label_data_cohort = paste(Label, data_cohort, sep = "_"))

  data <- data[, label$Sample]
  data <- as.data.frame(t(data))
  
  set.seed(1000)
  k_val <- 5
  times_val <- 6
  k_prod_times <- k_val * times_val
  train_index <- caret::createMultiFolds(y = label$Label_data_cohort, k = k_val, times = times_val)
  
  for (i in c(1:k_prod_times)){
    label.train <- label[train_index[[i]], , drop = FALSE]
    label.test <- label[-train_index[[i]], , drop = FALSE]

    data.train <- data[label.train$Sample, ]
    data.test <- data[label.test$Sample, ]
    
    result_list <- rf_model(data.train, label.train, data.test, label.test, 
                               classes, classifier_feature_imp = TRUE)
    result_df <- result_list[[1]]
    feature_importance <- result_list[[2]]
    
    result_df <- result_df %>%
      mutate(iter = i)
    feature_importance <- feature_importance %>%
      arrange(desc(MeanDecreaseGini)) %>%
      mutate(iter = i)
    
    if(i == 1){
      result_df_all <- result_df
      feature_importance_df_all <- feature_importance
    } else{
      result_df_all <- rbind(result_df_all, result_df)
      feature_importance_df_all <- rbind(feature_importance_df_all, feature_importance)
    }
  }
  
  feature_importance_df_all <- feature_importance_df_all %>%
    group_by(feature) %>%
    summarize(MeanDecreaseGini_mean = mean(MeanDecreaseGini))
  
  logFC_data <- read.csv(logfc_file_path) %>%
    dplyr::select(c(Molecule, logFC))
  
  feature_importance_df_all <- feature_importance_df_all %>%
    inner_join(logFC_data, by = c("feature" = "Molecule")) %>%
    arrange(desc(MeanDecreaseGini_mean))
  
  result_file_name <- paste0(file_name_prefix, "feature_importance_", comparison, ".csv")
  if(!dir.exists(result_file_dir_path)){
    dir.create(result_file_dir_path, recursive = TRUE)
  }
  write.csv(format(as.data.frame(feature_importance_df_all), digits = 3), paste0(result_file_dir_path, result_file_name), 
            row.names = FALSE)
}




rf_feature_importance(data_file_path = "Data/RNA_all/newquant_Nov2023_umi_counts_PREOPE_MET_HC_filter90_de_PREOPEVsMET.csv",
                      phenotype_file_path = "Data/transcriptomic_phenotype_PREOPE_MET_HC.txt",
                      logfc_file_path <- "plots_RNA_all/PREOPE_MET_HC/volcano/volcano_1_PREOPEVsMET.csv",
                      result_file_dir_path <- "plots_RNA_all/PREOPE_MET_HC/volcano/",
                      comparison = "PREOPEVsMET", classes = c("MET", "PREOPE"),
                      perform_filter = TRUE, norm = "log_cpm", file_name_prefix = "DE_")

rf_feature_importance(data_file_path = "Data/RNA_all/newquant_Nov2023_umi_counts_PREOPE_MET_HC_filter90_de_PREOPEVsHC.csv",
                      phenotype_file_path = "Data/transcriptomic_phenotype_PREOPE_MET_HC.txt",
                      logfc_file_path <- "plots_RNA_all/PREOPE_MET_HC/volcano/volcano_2_PREOPEVsHC.csv",
                      result_file_dir_path <- "plots_RNA_all/PREOPE_MET_HC/volcano/",
                      comparison = "PREOPEVsHC", classes = c("HC", "PREOPE"),
                      perform_filter = TRUE, norm = "log_cpm", file_name_prefix = "DE_")
