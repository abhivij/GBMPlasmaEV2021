library(tidyverse)
library(ggvenn)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

source("scripts/R/prediction_pipeline/cm_logistic_regression.R")
source("scripts/R/prediction_pipeline/cm_svm.R")
source("scripts/R/prediction_pipeline/cm_rf.R")

######################################################################
  
# comparison = "POSTOPE_TPVsREC_TP"
# omics_type = "transcriptomic"
# 
# #conditions : c(pos_class, neg_class, validation_class)
# conditions = c("POSTOPE_TP", "REC_TP", "PREREC")
# phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP"
# best_features_file_path = "Data/selected_features/best_features_with_add_col.csv"
# train_index = NA

comparison = "PREOPEVsREC_TP"
omics_type = "transcriptomic"
classes = c("REC_TP", "PREOPE")
best_features_file_path = "Data/selected_features/best_features_with_add_col.csv"
use_common_biomarkers = TRUE
result_file_name <- "Data/validation_prediction_result/PREOPEVsREC_TP.csv"
result_file_name <- "Data/validation_prediction_result/PREOPEVsREC_TP_commonbiomarkerswithtest.csv"

pipeline <- function(comparison, omics_type, classes, 
                     best_features_file_path, 
                     result_file_name = "Data/validation_prediction_result/PREOPEVsREC_TP.csv",
                     use_common_biomarkers = FALSE){

  best_features <- read.csv(best_features_file_path)  
  
  if(omics_type == "transcriptomic"){
    dataset_id <- paste0("GBM_tr_initial_",
                         comparison)
  }else{
    # dataset_id <- paste0("GBMPlasmaEV_proteomic_impute50fil_quantile_",
    #                      comparison)
  }
  best_features_sub <- best_features %>%
    filter(dataset_id == !!dataset_id,
           is_best == 1) 
  biomarkers <- strsplit(best_features_sub$biomarkers, split = "|", fixed = TRUE)[[1]]  
  
  if(omics_type == "transcriptomic"){
    data.train <- read.table("Data/RNA/umi_counts_initial_cohort.csv", header=TRUE, sep=",", row.names=1, skip=0,
                       nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")  
    data.test <- read.table("Data/RNA/umi_counts_validation_cohort.csv", header=TRUE, sep=",", row.names=1, skip=0,
                            nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
    colnames(data.test) <- paste0("S", colnames(data.test))
    norm <- "norm_log_cpm_simple"
    perform_filter <- TRUE
    split_str <- "simple_norm_"
    phenotype <- read.table("Data/transcriptomic_phenotype.txt", header=TRUE, sep="\t")
    test_metadata <- read.table("Data/RNA_validation/metadata_glionet.csv", header = TRUE, sep = ",")
    lim = c(-3, 5)
    breaks = seq(-3, 5, 1)
  } else {
    # norm <- "quantile"
    # perform_filter <- FALSE
    # split_str <- "quantile_"
    # phenotype <- read.table("Data/proteomic_phenotype.txt", header=TRUE, sep="\t")
    # if(grepl(pattern = "REC-TP", x = dataset_id, fixed = TRUE)){
    #   #currently this case will cause issues while reading phenotype 
    #   # and requiring multiple columns in phenotype
    #   data <- read.table("Data/Protein/formatted_data/Q7_nonorm_formatted_impute50fil.csv", header=TRUE, sep=",", row.names=1, skip=0,
    #                      nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
    # } else{
    #   data <- read.table("Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv", header=TRUE, sep=",", row.names=1, skip=0,
    #                      nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")  
    # }
    # lim = c(5, 17.5)
    # breaks = seq(5, 17.5, 2.5)
  } 
  
  output_labels <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes) %>%
    dplyr::select(Sample, Label)
  
  output_labels_test <- test_metadata %>%
    rename("Label" = "category_old_name") %>%
    filter(Label %in% c(classes, "UNK")) %>%
    rename("Sample" = "sample_id") %>%
    dplyr::select(Sample, Label)
  
  # if(omics_type == "proteomic"){
  #   #converting the proteomic labels to be same as proteomic label
  #   #NOTE : currently this is specific to POSTOPE_TP, REC_TP, PREREC samples
  #   output_labels <- output_labels %>%
  #     mutate(Sample = gsub("HB0", "HB", Sample, fixed = TRUE)) %>%
  #     filter(Sample != "HB6")
  #   
  #   colnames(data) <- gsub("HB0", "HB", colnames(data), fixed = TRUE)
  # }
  
  #currently data.train, data.test format : (transcripts x samples)
  data.train <- data.train %>% dplyr::select(output_labels$Sample)
  data.train <- as.data.frame(t(as.matrix(data.train)))
  
  data.test <- data.test %>% dplyr::select(output_labels_test$Sample)
  data.test <- as.data.frame(t(as.matrix(data.test)))
  
  #now data.train, data.test format : (samples x transcripts)
  
  
  ggvenn(list("selected biomarkers" = biomarkers, 
              "Validation Cohort before filter" = colnames(data.test)),
         fill_color = c("green", "blue"),
         stroke_size = 0.1,
         set_name_size = 5,
         text_size = 3)
  ggsave(paste0("plots/", comparison, "_beforefilter.png"))
  
  
  if(perform_filter){
    data.train <- as.data.frame(t(as.matrix(data.train)))
    data.test <- as.data.frame(t(as.matrix(data.test)))
    
    keep <- edgeR::filterByExpr(data.train, group = output_labels$Label)
    data.train <- data.train[keep, ]
    
    keep <- edgeR::filterByExpr(data.test, group = output_labels_test$Label)
    data.test <- data.test[keep, ]  
  }
  
  
  if(norm == "norm_log_cpm_simple"){
    #calculating norm log cpm
    print("norm_log_cpm_simple")
    
    data.train <- edgeR::cpm(data.train, log=TRUE)
    data.test <- edgeR::cpm(data.test, log=TRUE)
    
    data.train <- as.data.frame(t(as.matrix(data.train)))
    data.test <- as.data.frame(t(as.matrix(data.test)))
    
    #normalizing the data
    normparam <- caret::preProcess(data.train) 
    data.train <- predict(normparam, data.train)
    
    normparam <- caret::preProcess(data.test)
    data.test <- predict(normparam, data.test) 
    
  } else if(norm == "quantile"){
    
    # norm_data <- preprocessCore::normalize.quantiles(as.matrix(data.train))
    # norm_data <- data.frame(norm_data, row.names = rownames(data.train))
    # colnames(norm_data) <- colnames(data.train)
    # data.train <- norm_data
    # 
    # norm_data <- preprocessCore::normalize.quantiles(as.matrix(data.test))
    # norm_data <- data.frame(norm_data, row.names = rownames(data.test))
    # colnames(norm_data) <- colnames(data.test)
    # data.test <- norm_data
    # 
    # norm_data <- preprocessCore::normalize.quantiles(as.matrix(data.test2))
    # norm_data <- data.frame(norm_data, row.names = rownames(data.test2))
    # colnames(norm_data) <- colnames(data.test2)
    # data.test2 <- norm_data
  }
  #now data, data.test2 format : (samples x transcripts)
  
  ggvenn(list("selected biomarkers" = biomarkers, 
              "Validation Cohort after filter" = colnames(data.test)),
         fill_color = c("green", "orange"),
         stroke_size = 0.1,
         set_name_size = 5,
         text_size = 3)
  ggsave(paste0("plots/", comparison, "_afterfilter.png"))
  
  
  if(use_common_biomarkers){
    #get best biomarkers present in test_data
    common_biomarkers <- intersect(biomarkers, colnames(data.test))
    data.train <- data.train[, common_biomarkers]
    data.test <- data.test[, common_biomarkers]
  } else{
    #get all best biomarkers only
    data.train <- data.train[, biomarkers]
    
    missing_biomarkers <- setdiff(biomarkers, colnames(data.test))
    
    for(mb in missing_biomarkers){
      data.test[[mb]] <- 0  
    }
    
    data.test <- data.test[, biomarkers]    
  }

  
  if(omics_type == "transcriptomic"){
    #for now using only rf
    result_df <- rf_model(data.train, output_labels, data.test, output_labels_test, classes)

    
  } else if(omics_type == "proteomic"){
    # result_df <- svm_model(data.train, label.train, data.test, label.test, 
    #           data.test2, label.test2,
    #           classes, kernel = "sigmoid")
    # result_file_name <- "Data/prediction_result/proteomics.csv"
  }
  write.csv(format(result_df, digits = 3), result_file_name, row.names = FALSE)
  
}  

# train_index <- pipeline(comparison = "POSTOPE_TPVsREC_TP", omics_type = "transcriptomic", 
#                         conditions = c("POSTOPE_TP", "REC_TP", "PREREC"),
#                         phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
#                         best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
#                         train_index = NA)  
# pipeline(comparison = "POSTOPE_TPVsREC_TP", omics_type = "proteomic", 
#          conditions = c("POSTOPE_TP", "REC_TP", "PREREC"),
#          phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
#          best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
#          train_index = train_index)


pipeline(comparison = "PREOPEVsREC_TP", 
         omics_type = "transcriptomic", 
         classes = c("REC_TP", "PREOPE"),
         best_features_file_path = "Data/selected_features/best_features_with_add_col.csv")


pipeline(comparison = "PREOPEVsREC_TP", 
         omics_type = "transcriptomic", 
         classes = c("REC_TP", "PREOPE"),
         best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
         use_common_biomarkers = TRUE,
         result_file_name <- "Data/validation_prediction_result/PREOPEVsREC_TP_commonbiomarkerswithtest.csv")


result_df <- read.csv("Data/validation_prediction_result/PREOPEVsREC_TP.csv")
train_results <- result_df %>%
  filter(Type == "train")
acc.train <- sum(train_results$TrueLabel == train_results$PredictedLabel)/dim(train_results)[1]

pr <- ROCR::prediction(train_results$Probability, train_results$TrueLabel, label.ordering = classes)
auc <- ROCR::performance(pr, measure = "auc")@y.values[[1]]


test_results <- result_df %>%
  filter(Type == "test", TrueLabel != "UNK")
acc.test <- sum(test_results$TrueLabel == test_results$PredictedLabel) / dim(test_results)[1]

pr <- ROCR::prediction(test_results$Probability, test_results$TrueLabel, label.ordering = classes)
auc <- ROCR::performance(pr, measure = "auc")@y.values[[1]]



result_df <- read.csv("Data/validation_prediction_result/PREOPEVsREC_TP_commonbiomarkerswithtest.csv")
train_results <- result_df %>%
  filter(Type == "train")
acc.train <- sum(train_results$TrueLabel == train_results$PredictedLabel)/dim(train_results)[1]

pr <- ROCR::prediction(train_results$Probability, train_results$TrueLabel, label.ordering = classes)
auc <- ROCR::performance(pr, measure = "auc")@y.values[[1]]


test_results <- result_df %>%
  filter(Type == "test", TrueLabel != "UNK")
acc.test <- sum(test_results$TrueLabel == test_results$PredictedLabel) / dim(test_results)[1]

pr <- ROCR::prediction(test_results$Probability, test_results$TrueLabel, label.ordering = classes)
auc <- ROCR::performance(pr, measure = "auc")@y.values[[1]]
