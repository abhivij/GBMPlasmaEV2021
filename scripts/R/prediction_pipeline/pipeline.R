library(tidyverse)
base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

source("scripts/R/prediction_pipeline/cm_logistic_regression.R")
source("scripts/R/prediction_pipeline/cm_svm.R")

  ######################################################################
  
  comparison = "POSTOPE_TPVsREC_TP"
  omics_type = "transcriptomic"
  omics_type = "proteomic"
  
  #conditions : c(pos_class, neg_class, validation_class)
  conditions = c("POSTOPE_TP", "REC_TP", "PREREC")
  classes = conditions[1:2]
  
  phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP"
  
  best_features_file_path = "Data/selected_features/best_features_with_add_col.csv"
  
  train_index = NA
  
  
pipeline <- function(comparison, omics_type, conditions, classes, 
                     phenotype_column, best_features_file_path,
                     train_index = NA){
  
  best_features <- read.csv(best_features_file_path)  
  
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
  
  if(omics_type == "transcriptomic"){
    data <- read.table("Data/RNA/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                       nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")  
    norm <- "norm_log_cpm_simple"
    perform_filter <- TRUE
    split_str <- "simple_norm_"
    phenotype <- read.table("Data/transcriptomic_phenotype.txt", header=TRUE, sep="\t")
    lim = c(-3, 5)
    breaks = seq(-3, 5, 1)
  } else {
    norm <- "quantile"
    perform_filter <- FALSE
    split_str <- "quantile_"
    phenotype <- read.table("Data/proteomic_phenotype.txt", header=TRUE, sep="\t")
    if(grepl(pattern = "REC-TP", x = dataset_id, fixed = TRUE)){
      #currently this case will cause issues while reading phenotype 
      # and requiring multiple columns in phenotype
      data <- read.table("Data/Protein/formatted_data/Q7_nonorm_formatted_impute50fil.csv", header=TRUE, sep=",", row.names=1, skip=0,
                         nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
    } else{
      data <- read.table("Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv", header=TRUE, sep=",", row.names=1, skip=0,
                         nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")  
    }
    lim = c(5, 17.5)
    breaks = seq(5, 17.5, 2.5)
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
  
  #create 3 groups : train, test, test2
  
  #currently data format : (transcripts x samples)
  test2_class <- conditions[3]
  label.test2 <- output_labels %>%
    filter(Label == test2_class)
  data.test2 <- data %>% dplyr::select(label.test2$Sample)  
  data.test2 <- as.data.frame(t(as.matrix(data.test2)))
  
  output_labels <- output_labels %>%
    filter(Label != test2_class)
  data <- data %>% dplyr::select(output_labels$Sample)
  data <- as.data.frame(t(as.matrix(data)))
  
  #now data, data.test2 format : (samples x transcripts)

  if(is.na(train_index)){
    print("create data partition")
    random_seed = 1000
    set.seed(random_seed)
    train_index <- caret::createDataPartition(output_labels$Label, p = .8, 
                                              list = FALSE, 
                                              times = 1)    
  } else{
    print("using existing partition")
  }
  
  data.train <- data[train_index, ]
  label.train <- output_labels[train_index, ]
  
  data.test <- data[-train_index, ]
  label.test <- output_labels[-train_index, ]
  
  if(perform_filter){
    data.train <- as.data.frame(t(as.matrix(data.train)))
    data.test <- as.data.frame(t(as.matrix(data.test)))
    data.test2 <- as.data.frame(t(as.matrix(data.test2)))
    
    keep <- edgeR::filterByExpr(data.train, group = label.train$Label)
    data.train <- data.train[keep, ]
    data.test <- data.test[keep, ]  
    data.test2 <- data.test2[keep, ]
  }
  
  
  if(norm == "norm_log_cpm_simple"){
    #calculating norm log cpm
    print("norm_log_cpm_simple")
    
    data.train <- edgeR::cpm(data.train, log=TRUE)
    data.test <- edgeR::cpm(data.test, log=TRUE)
    data.test2 <- edgeR::cpm(data.test2, log=TRUE) 
    
    data.train <- as.data.frame(t(as.matrix(data.train)))
    data.test <- as.data.frame(t(as.matrix(data.test)))
    data.test2 <- as.data.frame(t(as.matrix(data.test2)))
    
    #normalizing the data
    normparam <- caret::preProcess(data.train) 
    data.train <- predict(normparam, data.train)
    data.test <- predict(normparam, data.test) #normalizing test data using params from train data   
    data.test2 <- predict(normparam, data.test2)
    
  } else if(norm == "quantile"){
    
    norm_data <- preprocessCore::normalize.quantiles(as.matrix(data.train))
    norm_data <- data.frame(norm_data, row.names = rownames(data.train))
    colnames(norm_data) <- colnames(data.train)
    data.train <- norm_data
    
    norm_data <- preprocessCore::normalize.quantiles(as.matrix(data.test))
    norm_data <- data.frame(norm_data, row.names = rownames(data.test))
    colnames(norm_data) <- colnames(data.test)
    data.test <- norm_data
    
    norm_data <- preprocessCore::normalize.quantiles(as.matrix(data.test2))
    norm_data <- data.frame(norm_data, row.names = rownames(data.test2))
    colnames(norm_data) <- colnames(data.test2)
    data.test2 <- norm_data
  }
  #now data, data.test2 format : (samples x transcripts)
  
  #get best biomarkers only
  data.train <- data.frame(data.train)[, biomarkers]  #data.frame() replaces - in colnames to .
  data.test <- data.frame(data.test)[, biomarkers]
  data.test2 <- data.frame(data.test2)[, biomarkers]
  
  if(omics_type == "transcriptomic"){
    result_df <- logistic_regression(data.train, label.train, data.test, label.test, 
                        data.test2, label.test2,
                        classes = classes, 
                        regularize = 'l2')    
    result_file_name <- "Data/prediction_result/transcriptomics.csv"
  } else if(omics_type == "proteomic"){
    result_df <- svm_model(data.train, label.train, data.test, label.test, 
              data.test2, label.test2,
              classes, kernel = "sigmoid")
    result_file_name <- "Data/prediction_result/proteomics.csv"
  }
  write.csv(format(result_df, digits = 3), result_file_name)
  
  return (train_index)
  
}  


pr_op_label <- output_labels
pr_op_label <- pr_op_label %>%
  mutate(Sample = gsub("HB0", "HB", Sample, fixed = TRUE)) %>%
  filter(Sample != "HB6")
all.equal(tr_op_label, pr_op_label)
