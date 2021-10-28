setwd("~/UNSW/VafaeeLab/GBMPlasmaEV/")
source("scripts/R/dataset_pipeline_arguments.R")
source("scripts/R/utils.R")
library(tidyverse)
library(viridis)
library(ComplexHeatmap)
library(UpSetR)


dataset_vec <- c("GBMPlasmaEV_transcriptomic_PREOPEVsMET",
                 "GBMPlasmaEV_transcriptomic_PREOPEVsHC",
                 "GBMPlasmaEV_transcriptomic_METVsHC",
                 "GBMPlasmaEV_proteomic_norm_quantile_PREOPEVsMET",
                 "GBMPlasmaEV_proteomic_norm_quantile_PREOPEVsHC",
                 "GBMPlasmaEV_proteomic_norm_quantile_METVsHC",
                 "GBMPlasmaEV_proteomic_impute75fil_norm_quantile_PREOPEVsMET",
                 "GBMPlasmaEV_proteomic_impute75fil_norm_quantile_PREOPEVsHC",
                 "GBMPlasmaEV_proteomic_impute75fil_norm_quantile_METVsHC",
                 "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsMET",
                 "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsHC",
                 "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_METVsHC"                 
)
dataset_vec <- dataset_vec[c(1:3, 10:12)]


fsm_vec <- c("all", 
             "t-test", "t-test_BH",
             "t-test_pval_0.025", "t-test_pval_0.01", "t-test_pval_0.005",
             "wilcoxontest", "wilcoxontest_BH",
             "wilcoxontest_pval_0.025", "wilcoxontest_pval_0.001", "wilcoxontest_pval_0.005",
             "ranger_impu_cor", 
             "mrmr10", "mrmr20",
             "mrmr30", "mrmr50", 
             "mrmr75", "mrmr100",
             "RF_RFE", "ga_rf")
# fsm_vector <- fsm_vec[c(1:2, 4, 6:12)]




# best_fsm <- "mrmr100"

best_fsm_vec <- c("ranger_impu_cor",
                  "t-test", "wilcoxontest",
                  "mrmr30", "mrmr100")

min_iter_feature_presence = 29
create_common_feature_plots_and_datasubsets <- function(dataset_id,
                                                        best_fsm_vec, 
                                                        min_iter_feature_presence){
  
  features_file <- paste(dataset_id, "features.csv", sep = "_")
  features_file <- paste("fem_pipeline_results", features_file, sep = "/")
  print(features_file)
  
  selected_features <- list()
  for(best_fsm in best_fsm_vec){
    features_info <- read.table(features_file, sep = ',', header = TRUE)
    
    features_info <- features_info %>%
      filter(FSM == best_fsm) %>%
      select(-c(FSM, Iter))
    print(best_fsm)
    print(dim(features_info[,colSums(features_info) >= min_iter_feature_presence]))
    
    selected_features[[best_fsm]] <- colnames(features_info[,colSums(features_info) >= min_iter_feature_presence])   
  }
  ###########################################
  #upsetplot
  file_name <- paste(dataset_id, min_iter_feature_presence, 
                     "upset.png", sep = "_")
  
  dir_path <- "plots/FEMPipeline/common_features_upset"
  if(!dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE)
  }
  
  print("creating upset plot...")
  print(paste(dir_path, file_name, sep = "/"))
  png(filename = paste(dir_path, file_name, sep = "/"), 
      units = "cm", width = 20, height = 15, res = 1200)
  upset(fromList(selected_features))
  dev.off()
  
  #create new data subsets
  
  data <- read.table("Data/RNA/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                     nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
  for(best_fsm in best_fsm_vec){
    data_sub <- data[gsub(".", "-", selected_features[[best_fsm]], fixed = TRUE),]
    print(dim(data_sub))
    print(sum(is.na(data_sub)))
    subset_filename <- paste0("Data/RNA/", dataset_id, "_", best_fsm, 
                              "_", min_iter_feature_presence,
                              "_umi_counts.csv")
    print(subset_filename)
    write.csv(data_sub, subset_filename)
  }
  
}

ds <- dataset_pipeline_arguments[[1]]
dataset_id <- paste(ds$dataset_id, ds$classification_criteria, sep = "_")
print(dataset_id) 

create_common_feature_plots_and_datasubsets(dataset_id = dataset_id,
                                            best_fsm_vec = best_fsm_vec,
                                            min_iter_feature_presence = 28)
create_common_feature_plots_and_datasubsets(dataset_id = dataset_id,
                                            best_fsm_vec = best_fsm_vec,
                                            min_iter_feature_presence = 29)
