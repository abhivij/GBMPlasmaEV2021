library(tidyverse)
library(readxl)
base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

source("scripts/R/prediction_pipeline/cm_logistic_regression.R")
source("scripts/R/prediction_pipeline/cm_svm.R")
source("scripts/R/prediction_pipeline/cm_rf.R")

######################################################################

comparison = "PREOPEVsPOSTOPE_TP"
# comparison = "POSTOPE_TPVsREC_TP"
# omics_type = "transcriptomic"
omics_type = "proteomic"

#conditions : c(pos_class, neg_class, validation_class)
# conditions = c("POSTOPE_TP", "REC_TP", "PREREC")

# conditions = c("POSTOPE_TP", "REC_TP")
conditions = c("POSTOPE_TP", "PREOPE")

phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP"

best_features_file_path = "Data/selected_features/made_old_2022June19/best_features_with_add_col.csv"

# train_index = NA
# 
# result_file_name <- "Data/prediction_result/transcriptomics_with_new_validation_data.csv"

result_file_name <- "Data/prediction_result/proteomics_with_new_validation_data.csv"
perform_filter = FALSE


comparison = "PREOPEVsPOSTOPE_TP"
omics_type = "proteomic"
conditions = c("POSTOPE_TP", "PREOPE")
phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP"
best_features_file_path = "Data/selected_features/best_features_with_add_col.csv"
result_file_name = "Data/prediction_result/proteomics_no_norm/PREOPEVsPOSTOPE_TP.csv"
model = "l2_log_reg"
norm = "no_norm"

pipeline_with_validation_data <- function(comparison, omics_type, conditions,
                                          phenotype_column, best_features_file_path,
                                          result_file_name, model, norm = NA,
                                          dataset_replace_str = NA,
                                          data_file_path = NA,
                                          validation_data_file_path = NA){
  
  classes = conditions
  best_features <- read.csv(best_features_file_path)  
  
  categories <- strsplit(comparison, split = "Vs", fixed = TRUE)[[1]]
  
  if(is.na(dataset_replace_str)){
    if(omics_type == "transcriptomic"){
      dataset_id <- paste0("GBMPlasmaEV_transcriptomic_simple_norm_",
                           comparison)
    }else{
      dataset_id <- paste0("GBM_initial_proteomic_impute50fil_", norm, "_",
                           comparison)
    }
    best_features_sub <- best_features %>%
      filter(dataset_id == !!dataset_id,
             is_best == 1)
  } else{
    best_features_sub <- best_features %>%
      mutate(dataset_id = gsub(dataset_replace_str, "", dataset_id)) %>%
      filter(is_best == 1, dataset_id == comparison)
  }
  biomarkers <- strsplit(best_features_sub$biomarkers, split = "|", fixed = TRUE)[[1]]  
  
  if(omics_type == "transcriptomic"){
    data <- read.table("Data/RNA/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                       nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")  
    norm <- "norm_log_cpm_simple"
    perform_filter <- TRUE
    phenotype <- read.table("Data/transcriptomic_phenotype.txt", header=TRUE, sep="\t")
    lim = c(-3, 5)
    breaks = seq(-3, 5, 1)
    
    validation_data <- read.table("Data/RNA_validation/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                                  nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
    validation_metadata <- read.csv("Data/RNA_validation/metadata_glionet.csv") %>%
      mutate(sample_category = factor(sample_category)) %>%
      mutate(sample_category = recode_factor(sample_category, "PRE-OP" = "PREOPE",
                                       "POST-OP" = "POSTOPE_TP",
                                       "RECURRENCE" = "REC_TP"))
  } else {
    if(is.na(norm)){
      norm <- "quantile"
    }
    perform_filter <- FALSE
    phenotype <- read.table("Data/proteomic_phenotype.txt", header=TRUE, sep="\t")
    if(is.na(data_file_path)){
      data_file_path <- "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv" 
    }
    if(is.na(validation_data_file_path)){
      validation_data_file_path <- "Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil.csv"
    }
    
    data <- read.table(data_file_path, header=TRUE, sep=",", row.names=1, skip=0,
                         nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")  
    lim = c(5, 17.5)
    breaks = seq(5, 17.5, 2.5)
    validation_data <- read.csv(validation_data_file_path, row.names = 1)
    validation_metadata <- read.csv("Data/RNA_validation/metadata_glionet.csv") %>%
      mutate(sample_category = category_old_name) %>%
      select(-c(category_old_name))
    
    colnames(validation_data)[colnames(validation_data) == "SB12_01"] = "SB12"
    
    #use SB22.02
    colnames(validation_data)[colnames(validation_data) == "SB22.02"] = "SBtobeused22"
    colnames(validation_data)[colnames(validation_data) == "SB22"] = "SB22_dont_include"
    colnames(validation_data)[colnames(validation_data) == "SBtobeused22"] = "SB22"
    
    validation_metadata <- validation_metadata %>%
      filter(sample_id != "SB7")
  } 
  
  output_labels <- phenotype %>%
    rename("Label" = phenotype_column) %>%
    filter(Label %in% conditions) %>%
    dplyr::select(Sample, Label)
  
  if(omics_type == "proteomic"){
    #converting the proteomic labels to be same as proteomic label
    #NOTE : currently this is specific to POSTOPE_TP, REC_TP, PREREC samples
    output_labels <- output_labels %>%
      mutate(Sample = gsub("HB0", "HB", Sample, fixed = TRUE)) %>%
      filter(Sample != "HB6")
    
    colnames(data) <- gsub("HB0", "HB", colnames(data), fixed = TRUE)
  }
  #currently data format : (transcripts x samples)

  data <- data %>% dplyr::select(output_labels$Sample)
  data <- as.data.frame(t(as.matrix(data)))
  
  #now data format : (samples x transcripts)

  data.train <- data
  label.train <- output_labels
  
  #TODO : get new data and labels as data.test
  label.test <- validation_metadata %>%
    select(sample_id, sample_category) %>%
    filter(sample_category %in% c(conditions, "UNK")) %>%
    arrange(sample_category)
  colnames(label.test) <- colnames(label.train)
  
  data.test <- validation_data %>% dplyr::select(label.test$Sample)
  data.test <- as.data.frame(t(as.matrix(data.test)))
  
  data.train <- as.data.frame(t(as.matrix(data.train)))
  data.test <- as.data.frame(t(as.matrix(data.test)))
  if(perform_filter){
    
    keep <- edgeR::filterByExpr(data.train, group = label.train$Label)
    data.train <- data.train[keep, ]
    
    keep <- edgeR::filterByExpr(data.test, group = label.test$Label)
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
    data.test <- predict(normparam, data.test) #normalizing test data using params from train data   
    
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
  } else if(norm == "quantile_train_param"){
    #adapted from https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
    data.train.rank <- apply(data.train, 2, rank, ties.method="average")
    data.train.sorted <- data.frame(apply(data.train, 2, sort))
    data.train.mean <- apply(data.train.sorted, 1, mean)
    index_to_mean <- function(index, data_mean){
      #index can be int or int+0.5
      #if int+0.5, take average of the numbers in those positions
      int.result <- data_mean[index]
      index.int <- floor(index)
      #some of the values in point5.result might be NA
      #but they won't be chosen
      point5.result <- (data_mean[index.int] + data_mean[index.int+1])/2
      point5.indices <- index%%1 != 0
      result <- int.result
      result[point5.indices] <- point5.result[point5.indices]
      return (result)
    }
    data.train.norm <- apply(data.train.rank, 2, index_to_mean, data_mean = data.train.mean)
    rownames(data.train.norm) <- rownames(data.train)
    data.train <- data.train.norm
    
    data.test.rank <- apply(data.test, 2, rank, ties.method="average")
    #use params i.e. mean values of rows, from training data
    data.test.norm <- apply(data.test.rank, 2, index_to_mean, data_mean = data.train.mean)
    rownames(data.test.norm) <- rownames(data.test)
    data.test <- data.test.norm
  }
  data.train <- as.data.frame(t(as.matrix(data.train)))
  data.test <- as.data.frame(t(as.matrix(data.test)))
  #now data, data.test2 format : (samples x transcripts)
  
  #get best biomarkers only
  data.train <- data.frame(data.train)[, biomarkers]  #data.frame() replaces - in colnames to .
  
  colnames(data.test) <- gsub("-", ".", colnames(data.test))
  
  sum(biomarkers %in% colnames(data.test))
  
  available_biomarkers <- c(biomarkers[biomarkers %in% colnames(data.test)])
  print(available_biomarkers)
  non_available_biomarkers <- biomarkers[!biomarkers %in% colnames(data.test)] 
  print(non_available_biomarkers)
  
  print(length(available_biomarkers))
  print(length(non_available_biomarkers))
  
  data.test <- data.test %>%
    select(available_biomarkers)
  for(nab in non_available_biomarkers){
    data.test[[nab]] <- 0
  }
  
  # data.test2 = NA
  # label.test2 = NA
  # regularize = 'l2'
  if(model == "rf"){
    result_df <- rf_model(data.train, label.train, data.test, label.test, classes)
  } else if(model == "radial_svm"){
    result_df <- svm_model(data.train, label.train,
                           data.test, label.test,
                           classes, kernel = "radial")
  } else if(model == "sigmoid_svm"){
    result_df <- svm_model(data.train, label.train,
                           data.test, label.test,
                           classes, kernel = "sigmoid")
  } else if(model == "l2_log_reg"){
    result_df <- log_reg_model(data.train, label.train, data.test, label.test,
                               classes, regularize = "l2")
  } else if(model == "l1_log_reg"){
    result_df <- log_reg_model(data.train, label.train, data.test, label.test,
                               classes, regularize = "l1")
  }
  write.csv(format(result_df, digits = 3), result_file_name, row.names = FALSE)
}  


comparison = "PREOPEVsREC_TP"
classes = c("REC_TP", "PREOPE")
result_file_path = "Data/prediction_result/proteomics_no_norm/PREOPEVsREC_TP.csv"
metric_output_file_path = "Data/prediction_result/proteomics_no_norm/metrics.csv"

show_metrics <- function(comparison, classes, result_file_path,
                         metric_output_file_path = "data/prediction_result/metrics.csv"){
  print(comparison)
  result_df <- read.csv(result_file_path) %>%
    filter(TrueLabel %in% classes)
  
  results.train <- result_df %>%
    filter(Type == "train")
  results.test <- result_df %>%
    filter(Type == "test")
  
  acc.train <- mean(results.train$TrueLabel == results.train$PredictedLabel)
  acc.test <- mean(results.test$TrueLabel == results.test$PredictedLabel)
  
  pr <- ROCR::prediction(results.train$Pred_prob, results.train$TrueLabel, label.ordering = classes)
  auc.train <- ROCR::performance(pr, measure = "auc")@y.values[[1]]
  
  pr <- ROCR::prediction(results.test$Pred_prob, results.test$TrueLabel, label.ordering = classes)
  auc.test <- ROCR::performance(pr, measure = "auc")@y.values[[1]]
  
  #classes[2] assigned to class1 deliberately
  #classes are provided as c(neg_class, pos_class)
  class1 <- classes[2]
  class2 <- classes[1]
  class1_count <- results.train %>% dplyr::filter(TrueLabel == class1) %>% nrow()
  class2_count <- results.train %>% dplyr::filter(TrueLabel == class2) %>% nrow() 
  train_count <- paste0(class1, "=", class1_count, " | ", class2, "=", class2_count)
  class1_count <- results.test %>% dplyr::filter(TrueLabel == class1) %>% nrow()
  class2_count <- results.test %>% dplyr::filter(TrueLabel == class2) %>% nrow() 
  test_count <- paste0(class1, "=", class1_count, " | ", class2, "=", class2_count)
  metrics <- data.frame(Comparison = comparison,
                        TrainCount = train_count,
                        TestCount = test_count,
                        TrainAccuracy = acc.train,
                        TrainAUC = auc.train,
                        TestAccuracy = acc.test,
                        TestAUC = auc.test)
  
  write.table(x = metrics, file = metric_output_file_path, append = TRUE, 
              col.names = !file.exists(metric_output_file_path), sep = ",",
              row.names = FALSE)
}

pipeline_with_validation_data(comparison = "POSTOPE_TPVsREC_TP", omics_type = "transcriptomic", 
                        conditions = c("POSTOPE_TP", "REC_TP"),
                        phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                        best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                        result_file_name = "Data/prediction_result/transcriptomics_POSTOPE_TPVsREC_TP_with_new_validation_data.csv") 






#################
#check number of transcripts


best_features <- read.csv("Data/selected_features/best_features_with_add_col.csv")  
comparison = "POSTOPE_TPVsREC_TP"

categories <- strsplit(comparison, split = "Vs", fixed = TRUE)[[1]]
if(omics_type == "transcriptomic"){
  dataset_id <- paste0("GBMPlasmaEV_transcriptomic_simple_norm_",
                       comparison)
}else{
  dataset_id <- paste0("GBMPlasmaEV_proteomic_impute50fil_quantile_",
                       comparison)
}
best_features_sub <- best_features %>%
  filter(dataset_id == !!dataset_id,
         is_best == 1) 
biomarkers <- strsplit(best_features_sub$biomarkers, split = "|", fixed = TRUE)[[1]]  


validation_data <- read.table("Data/RNA_validation/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                              nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
rownames(validation_data)[1:10]
transcripts_identified <- gsub("-", ".", rownames(validation_data), fixed = TRUE)
biomarkers[biomarkers %in% transcripts_identified]

validation_data[gsub(".", "-", c("hsa.miR.1307.3p", "hsa.miR.328.3p", "hsa.miR.4513"), fixed = TRUE),]

data <- read_excel("Data/qiagen_results_test/GBMPlasmaEV_initial_cohort_RNAData/219155.all_samples.summary.xlsx", 
                   sheet = "miRNA_piRNA")
mirna_data <- data[1:324,]
pirna_data <- data[326:461,] %>%
  separate(miRNA, c("miRNA", NA, NA, NA), sep = "/")

dim(data)[1]
dim(mirna_data)[1] + dim(pirna_data)[1]

data <- rbind(mirna_data, pirna_data) %>%
  column_to_rownames("miRNA")

umi_counts_a <- data %>%
  select(ends_with("UMIs"))
colnames(umi_counts_a) <- gsub("-UMIs", "", colnames(umi_counts_a))
colnames(umi_counts_a) <- sapply(colnames(umi_counts_a), FUN = 
                                   function(x){
                                     strsplit(x, split = "_", fixed = TRUE)[[1]][1]
                                   }
)

transcripts_identified <- gsub("-", ".", rownames(umi_counts_a), fixed = TRUE)
biomarkers[biomarkers %in% transcripts_identified]

umi_counts_a[gsub(".", "-", c("hsa.miR.1307.3p", "hsa.miR.328.3p", "hsa.miR.4513"), fixed = TRUE),]



###############################
#proteomics

pipeline_with_validation_data(comparison = "PREOPEVsPOSTOPE_TP", 
                              omics_type = "proteomic", 
                              conditions = c("POSTOPE_TP", "PREOPE"),
                              phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                              best_features_file_path = "Data/selected_features/made_old_2022June19/best_features_with_add_col.csv",
                              result_file_name = "Data/prediction_result/proteomics/PREOPEVsPOSTOPE_TP.csv")

pipeline_with_validation_data(comparison = "PREOPEVsREC_TP", 
                              omics_type = "proteomic", 
                              conditions = c("REC_TP", "PREOPE"),
                              phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                              best_features_file_path = "Data/selected_features/made_old_2022June19/best_features_with_add_col.csv",
                              result_file_name = "Data/prediction_result/proteomics/PREOPEVsREC_TP.csv")

pipeline_with_validation_data(comparison = "POSTOPE_TPVsREC_TP", 
                              omics_type = "proteomic", 
                              conditions = c("REC_TP", "POSTOPE_TP"),
                              phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                              best_features_file_path = "Data/selected_features/made_old_2022June19/best_features_with_add_col.csv",
                              result_file_name = "Data/prediction_result/proteomics/POSTOPE_TPVsREC_TP.csv")


show_metrics(comparison = "PREOPEVsPOSTOPE_TP", classes = c("POSTOPE_TP", "PREOPE"), 
             result_file_path = "Data/prediction_result/proteomics/PREOPEVsPOSTOPE_TP.csv",
             metric_output_file_path = "Data/prediction_result/proteomics/metrics.csv")
show_metrics(comparison = "PREOPEVsREC_TP", classes = c("REC_TP", "PREOPE"), 
             result_file_path = "Data/prediction_result/proteomics/PREOPEVsREC_TP.csv",
             metric_output_file_path = "Data/prediction_result/proteomics/metrics.csv")
show_metrics(comparison = "POSTOPE_TPVsREC_TP", classes = c("REC_TP", "POSTOPE_TP"), 
             result_file_path = "Data/prediction_result/proteomics/POSTOPE_TPVsREC_TP.csv",
             metric_output_file_path = "Data/prediction_result/proteomics/metrics.csv")



#######


pipeline_with_validation_data(comparison = "PREOPEVsPOSTOPE_TP", 
                              omics_type = "proteomic", 
                              conditions = c("POSTOPE_TP", "PREOPE"),
                              phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                              best_features_file_path = "Data/selected_features/made_old_2022June19/best_features_with_add_col.csv",
                              result_file_name = "Data/prediction_result/proteomics_nonorm/PREOPEVsPOSTOPE_TP.csv")

pipeline_with_validation_data(comparison = "PREOPEVsREC_TP", 
                              omics_type = "proteomic", 
                              conditions = c("REC_TP", "PREOPE"),
                              phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                              best_features_file_path = "Data/selected_features/made_old_2022June19/best_features_with_add_col.csv",
                              result_file_name = "Data/prediction_result/proteomics_nonorm/PREOPEVsREC_TP.csv")

pipeline_with_validation_data(comparison = "POSTOPE_TPVsREC_TP", 
                              omics_type = "proteomic", 
                              conditions = c("REC_TP", "POSTOPE_TP"),
                              phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                              best_features_file_path = "Data/selected_features/made_old_2022June19/best_features_with_add_col.csv",
                              result_file_name = "Data/prediction_result/proteomics_nonorm/POSTOPE_TPVsREC_TP.csv")


show_metrics(comparison = "PREOPEVsPOSTOPE_TP", classes = c("POSTOPE_TP", "PREOPE"), 
             result_file_path = "Data/prediction_result/proteomics_nonorm/PREOPEVsPOSTOPE_TP.csv",
             metric_output_file_path = "Data/prediction_result/proteomics_nonorm/metrics.csv")
show_metrics(comparison = "PREOPEVsREC_TP", classes = c("REC_TP", "PREOPE"), 
             result_file_path = "Data/prediction_result/proteomics_nonorm/PREOPEVsREC_TP.csv",
             metric_output_file_path = "Data/prediction_result/proteomics_nonorm/metrics.csv")
show_metrics(comparison = "POSTOPE_TPVsREC_TP", classes = c("REC_TP", "POSTOPE_TP"), 
             result_file_path = "Data/prediction_result/proteomics_nonorm/POSTOPE_TPVsREC_TP.csv",
             metric_output_file_path = "Data/prediction_result/proteomics_nonorm/metrics.csv")



###########



pipeline_with_validation_data(comparison = "PREOPEVsPOSTOPE_TP", 
                              omics_type = "proteomic", 
                              conditions = c("POSTOPE_TP", "PREOPE"),
                              phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                              best_features_file_path = "Data/selected_features/made_old_2022June19/best_features_with_add_col.csv",
                              result_file_name = "Data/prediction_result/proteomics_trainparamnorm/PREOPEVsPOSTOPE_TP.csv")

pipeline_with_validation_data(comparison = "PREOPEVsREC_TP", 
                              omics_type = "proteomic", 
                              conditions = c("REC_TP", "PREOPE"),
                              phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                              best_features_file_path = "Data/selected_features/made_old_2022June19/best_features_with_add_col.csv",
                              result_file_name = "Data/prediction_result/proteomics_trainparamnorm/PREOPEVsREC_TP.csv")

pipeline_with_validation_data(comparison = "POSTOPE_TPVsREC_TP", 
                              omics_type = "proteomic", 
                              conditions = c("REC_TP", "POSTOPE_TP"),
                              phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                              best_features_file_path = "Data/selected_features/made_old_2022June19/best_features_with_add_col.csv",
                              result_file_name = "Data/prediction_result/proteomics_trainparamnorm/POSTOPE_TPVsREC_TP.csv")


show_metrics(comparison = "PREOPEVsPOSTOPE_TP", classes = c("POSTOPE_TP", "PREOPE"), 
             result_file_path = "Data/prediction_result/proteomics_trainparamnorm/PREOPEVsPOSTOPE_TP.csv",
             metric_output_file_path = "Data/prediction_result/proteomics_trainparamnorm/metrics.csv")
show_metrics(comparison = "PREOPEVsREC_TP", classes = c("REC_TP", "PREOPE"), 
             result_file_path = "Data/prediction_result/proteomics_trainparamnorm/PREOPEVsREC_TP.csv",
             metric_output_file_path = "Data/prediction_result/proteomics_trainparamnorm/metrics.csv")
show_metrics(comparison = "POSTOPE_TPVsREC_TP", classes = c("REC_TP", "POSTOPE_TP"), 
             result_file_path = "Data/prediction_result/proteomics_trainparamnorm/POSTOPE_TPVsREC_TP.csv",
             metric_output_file_path = "Data/prediction_result/proteomics_trainparamnorm/metrics.csv")



###################
pipeline_with_validation_data(comparison = "PREOPEVsPOSTOPE_TP", 
                              omics_type = "proteomic", 
                              conditions = c("POSTOPE_TP", "PREOPE"),
                              phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                              best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                              result_file_name = "Data/prediction_result/proteomics_no_norm/PREOPEVsPOSTOPE_TP.csv",
                              model = "l2_log_reg", norm = "no_norm")

pipeline_with_validation_data(comparison = "POSTOPE_TPVsREC_TP", 
                              omics_type = "proteomic", 
                              conditions = c("REC_TP", "POSTOPE_TP"),
                              phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                              best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                              result_file_name = "Data/prediction_result/proteomics_no_norm/POSTOPE_TPVsREC_TP.csv",
                              model = "l2_log_reg", norm = "no_norm")

pipeline_with_validation_data(comparison = "PREOPEVsREC_TP", 
                              omics_type = "proteomic", 
                              conditions = c("REC_TP", "PREOPE"),
                              phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                              best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                              result_file_name = "Data/prediction_result/proteomics_no_norm/PREOPEVsREC_TP.csv",
                              model = "radial_svm", norm = "no_norm")

show_metrics(comparison = "PREOPEVsPOSTOPE_TP", classes = c("POSTOPE_TP", "PREOPE"), 
             result_file_path = "Data/prediction_result/proteomics_no_norm/PREOPEVsPOSTOPE_TP.csv",
             metric_output_file_path = "Data/prediction_result/proteomics_no_norm/metrics.csv")
show_metrics(comparison = "POSTOPE_TPVsREC_TP", classes = c("REC_TP", "POSTOPE_TP"), 
             result_file_path = "Data/prediction_result/proteomics_no_norm/POSTOPE_TPVsREC_TP.csv",
             metric_output_file_path = "Data/prediction_result/proteomics_no_norm/metrics.csv")
show_metrics(comparison = "PREOPEVsREC_TP", classes = c("REC_TP", "PREOPE"), 
             result_file_path = "Data/prediction_result/proteomics_no_norm/PREOPEVsREC_TP.csv",
             metric_output_file_path = "Data/prediction_result/proteomics_no_norm/metrics.csv")

pipeline_with_validation_data(comparison = "PREOPEVsPOSTOPE_TP", 
                              omics_type = "proteomic", 
                              conditions = c("POSTOPE_TP", "PREOPE"),
                              phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                              best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                              result_file_name = "Data/prediction_result/proteomics_quantile_train_param/PREOPEVsPOSTOPE_TP.csv",
                              model = "l2_log_reg", norm = "quantile_train_param")

pipeline_with_validation_data(comparison = "POSTOPE_TPVsREC_TP", 
                              omics_type = "proteomic", 
                              conditions = c("REC_TP", "POSTOPE_TP"),
                              phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                              best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                              result_file_name = "Data/prediction_result/proteomics_quantile_train_param/POSTOPE_TPVsREC_TP.csv",
                              model = "radial_svm", norm = "quantile_train_param")

pipeline_with_validation_data(comparison = "PREOPEVsREC_TP", 
                              omics_type = "proteomic", 
                              conditions = c("REC_TP", "PREOPE"),
                              phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                              best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                              result_file_name = "Data/prediction_result/proteomics_quantile_train_param/PREOPEVsREC_TP.csv",
                              model = "l2_log_reg", norm = "quantile_train_param")


show_metrics(comparison = "PREOPEVsPOSTOPE_TP", classes = c("POSTOPE_TP", "PREOPE"), 
             result_file_path = "Data/prediction_result/proteomics_quantile_train_param/PREOPEVsPOSTOPE_TP.csv",
             metric_output_file_path = "Data/prediction_result/proteomics_quantile_train_param/metrics.csv")
show_metrics(comparison = "POSTOPE_TPVsREC_TP", classes = c("REC_TP", "POSTOPE_TP"), 
             result_file_path = "Data/prediction_result/proteomics_quantile_train_param/POSTOPE_TPVsREC_TP.csv",
             metric_output_file_path = "Data/prediction_result/proteomics_quantile_train_param/metrics.csv")
show_metrics(comparison = "PREOPEVsREC_TP", classes = c("REC_TP", "PREOPE"), 
             result_file_path = "Data/prediction_result/proteomics_quantile_train_param/PREOPEVsREC_TP.csv",
             metric_output_file_path = "Data/prediction_result/proteomics_quantile_train_param/metrics.csv")


####################
#common proteins no norm
###################
pipeline_with_validation_data(comparison = "POSTOPE_TPVsREC_TP", 
                              omics_type = "proteomic", 
                              conditions = c("REC_TP", "POSTOPE_TP"),
                              phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                              best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                              result_file_name = "Data/prediction_result/proteomics_common_no_norm/POSTOPE_TPVsREC_TP.csv",
                              model = "sigmoid_svm", norm = "no_norm",
                              dataset_replace_str = "GBM_initial_proteomic_impute50fil_common_no_norm_",
                              data_file_path = "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil_common.csv", 
                              validation_data_file_path = "Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil_common.csv")
pipeline_with_validation_data(comparison = "PREOPEVsPOSTOPE_TP", 
                              omics_type = "proteomic", 
                              conditions = c("POSTOPE_TP", "PREOPE"),
                              phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                              best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                              result_file_name = "Data/prediction_result/proteomics_common_no_norm/PREOPEVsPOSTOPE_TP.csv",
                              model = "radial_svm", norm = "no_norm",
                              dataset_replace_str = "GBM_initial_proteomic_impute50fil_common_no_norm_",
                              data_file_path = "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil_common.csv", 
                              validation_data_file_path = "Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil_common.csv")
pipeline_with_validation_data(comparison = "PREOPEVsREC_TP", 
                              omics_type = "proteomic", 
                              conditions = c("REC_TP", "PREOPE"),
                              phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                              best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                              result_file_name = "Data/prediction_result/proteomics_common_no_norm/PREOPEVsREC_TP.csv",
                              model = "rf", norm = "no_norm",
                              dataset_replace_str = "GBM_initial_proteomic_impute50fil_common_no_norm_",
                              data_file_path = "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil_common.csv", 
                              validation_data_file_path = "Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil_common.csv")

show_metrics(comparison = "POSTOPE_TPVsREC_TP", classes = c("REC_TP", "POSTOPE_TP"), 
             result_file_path = "Data/prediction_result/proteomics_common_no_norm/POSTOPE_TPVsREC_TP.csv",
             metric_output_file_path = "Data/prediction_result/proteomics_common_no_norm/metrics.csv")
show_metrics(comparison = "PREOPEVsPOSTOPE_TP", classes = c("POSTOPE_TP", "PREOPE"), 
             result_file_path = "Data/prediction_result/proteomics_common_no_norm/PREOPEVsPOSTOPE_TP.csv",
             metric_output_file_path = "Data/prediction_result/proteomics_common_no_norm/metrics.csv")
show_metrics(comparison = "PREOPEVsREC_TP", classes = c("REC_TP", "PREOPE"), 
             result_file_path = "Data/prediction_result/proteomics_common_no_norm/PREOPEVsREC_TP.csv",
             metric_output_file_path = "Data/prediction_result/proteomics_common_no_norm/metrics.csv")


#########################
#common proteins quantile train param norm
pipeline_with_validation_data(comparison = "POSTOPE_TPVsREC_TP", 
                              omics_type = "proteomic", 
                              conditions = c("REC_TP", "POSTOPE_TP"),
                              phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                              best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                              result_file_name = "Data/prediction_result/proteomics_common_quantile_train_param/POSTOPE_TPVsREC_TP.csv",
                              model = "sigmoid_svm", norm = "quantile_train_param",
                              dataset_replace_str = "GBM_initial_proteomic_impute50fil_common_quantile_train_param_",
                              data_file_path = "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil_common.csv", 
                              validation_data_file_path = "Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil_common.csv")
pipeline_with_validation_data(comparison = "PREOPEVsPOSTOPE_TP", 
                              omics_type = "proteomic", 
                              conditions = c("POSTOPE_TP", "PREOPE"),
                              phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                              best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                              result_file_name = "Data/prediction_result/proteomics_common_quantile_train_param/PREOPEVsPOSTOPE_TP.csv",
                              model = "radial_svm", norm = "quantile_train_param",
                              dataset_replace_str = "GBM_initial_proteomic_impute50fil_common_quantile_train_param_",
                              data_file_path = "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil_common.csv", 
                              validation_data_file_path = "Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil_common.csv")
pipeline_with_validation_data(comparison = "PREOPEVsREC_TP", 
                              omics_type = "proteomic", 
                              conditions = c("REC_TP", "PREOPE"),
                              phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP", 
                              best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                              result_file_name = "Data/prediction_result/proteomics_common_quantile_train_param/PREOPEVsREC_TP.csv",
                              model = "radial_svm", norm = "quantile_train_param",
                              dataset_replace_str = "GBM_initial_proteomic_impute50fil_common_quantile_train_param_",
                              data_file_path = "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil_common.csv", 
                              validation_data_file_path = "Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil_common.csv")

show_metrics(comparison = "POSTOPE_TPVsREC_TP", classes = c("REC_TP", "POSTOPE_TP"), 
             result_file_path = "Data/prediction_result/proteomics_common_quantile_train_param/POSTOPE_TPVsREC_TP.csv",
             metric_output_file_path = "Data/prediction_result/proteomics_common_quantile_train_param/metrics.csv")
show_metrics(comparison = "PREOPEVsPOSTOPE_TP", classes = c("POSTOPE_TP", "PREOPE"), 
             result_file_path = "Data/prediction_result/proteomics_common_quantile_train_param/PREOPEVsPOSTOPE_TP.csv",
             metric_output_file_path = "Data/prediction_result/proteomics_common_quantile_train_param/metrics.csv")
show_metrics(comparison = "PREOPEVsREC_TP", classes = c("REC_TP", "PREOPE"), 
             result_file_path = "Data/prediction_result/proteomics_common_quantile_train_param/PREOPEVsREC_TP.csv",
             metric_output_file_path = "Data/prediction_result/proteomics_common_quantile_train_param/metrics.csv")
