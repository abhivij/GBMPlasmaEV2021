setwd("~/UNSW/VafaeeLab/GBMPlasmaEV/")
source("scripts/R/dataset_pipeline_arguments.R")
source("scripts/R/utils.R")
setwd("~/UNSW/VafaeeLab/GBMPlasmaEV/fem_pipeline_results/")
library(tidyverse)
library(viridis)
library(ComplexHeatmap)
library(UpSetR)


dataset_vec <- c("GBMPlasmaEV_transcriptomic_PREOPEVsMET",
                 "GBMPlasmaEV_transcriptomic_PREOPEVsHC",
                 "GBMPlasmaEV_transcriptomic_METVsHC",
                 "GBMPlasmaEV_proteomic_norm_quantile_PREOPEVsMET",
                 "GBMPlasmaEV_proteomic_norm_quantile_PREOPEVsHC",
                 "GBMPlasmaEV_proteomic_norm_quantile_METVsHC"
)


min_iter_feature_presence <- 28

fsm_vec <- c("all", "t-test", "t-test_BH",
             "wilcoxontest", "wilcoxontest_BH",
             "ranger_impu_cor", 
             "mrmr30", "mrmr50", "mrmr75", "mrmr100", "RF_RFE", "ga_rf")
fsm_vector <- fsm_vec[c(1:2, 4, 6:12)]

ds <- dataset_pipeline_arguments[[1]]
dataset_id <- paste(ds$dataset_id, ds$classification_criteria, sep = "_")
print(dataset_id) 

features_file <- paste(dataset_id, "features.csv", sep = "_")

# best_fsm <- "mrmr100"

best_fsm_vec <- c("mrmr30", "mrmr100", "wilcoxontest", "ranger_impu_cor")

selected_features <- list()
for(best_fsm in best_fsm_vec){
  features_info <- read.table(features_file, sep = ',', header = TRUE)
  
  features_info <- features_info %>%
    filter(FSM == best_fsm) %>%
    select(-c(FSM, Iter))
  
  selected_features[[best_fsm]] <- colnames(features_info[,colSums(features_info) >= 28])   
}
###########################################
#upsetplot

upset(fromList(selected_features))

VennDiagram::venn.diagram(selected_features, "plots/venn_selected_features.png")

###########################################
setwd("~/UNSW/VafaeeLab/GBMPlasmaEV/")
data <- read.table("Data/RNA/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                   nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
for(best_fsm in best_fsm_vec){
  data_sub <- data[gsub(".", "-", selected_features[[best_fsm]], fixed = TRUE),]
  print(dim(data_sub))
  print(sum(is.na(data_sub)))
  subset_filename <- paste0("Data/RNA/", dataset_id, "_", best_fsm, "_umi_counts.csv")
  print(subset_filename)
  write.csv(data_sub, subset_filename)
}


