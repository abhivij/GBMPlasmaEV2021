library(tidyverse)
library(readxl)
base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

source("scripts/R/prediction_pipeline/cm_logistic_regression.R")
source("scripts/R/prediction_pipeline/cm_svm.R")

######################################################################

comparison = "POSTOPE_TPVsREC_TP"
omics_type = "transcriptomic"
# omics_type = "proteomic"

#conditions : c(pos_class, neg_class, validation_class)
# conditions = c("POSTOPE_TP", "REC_TP", "PREREC")

conditions = c("POSTOPE_TP", "REC_TP")

phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP"

best_features_file_path = "Data/selected_features/best_features_with_add_col.csv"

train_index = NA

result_file_name <- "Data/prediction_result/transcriptomics_with_new_validation_data.csv"

result_file_name <- "Data/prediction_result/proteomics_with_new_validation_data.csv"


pipeline_with_validation_data <- function(comparison, omics_type, conditions,
                     phenotype_column, best_features_file_path,
                     result_file_name){
  
  classes = conditions
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
    
    validation_data <- read.table("Data/RNA_validation/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                                  nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
    validation_metadata <- read.csv("Data/RNA_validation/metadata_glionet.csv") %>%
      mutate(sample_category = factor(sample_category)) %>%
      mutate(sample_category = recode_factor(sample_category, "PRE-OP" = "PREOPE",
                                       "POST-OP" = "POSTOPE_TP",
                                       "RECURRENCE" = "REC_TP"))
    
    
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
  
  
  #currently data format : (transcripts x samples)

  data <- data %>% dplyr::select(output_labels$Sample)
  data <- as.data.frame(t(as.matrix(data)))
  
  #now data, data.test2 format : (samples x transcripts)

  data.train <- data
  label.train <- output_labels
  
  #TODO : get new data and labels as data.test
  label.test <- validation_metadata %>%
    select(sample_id, sample_category) %>%
    filter(sample_category %in% c(conditions, NA)) %>%
    arrange(sample_category)
  colnames(label.test) <- colnames(label.train)
  
  data.test <- validation_data %>% dplyr::select(label.test$Sample)
  data.test <- as.data.frame(t(as.matrix(data.test)))
  
  if(perform_filter){
    data.train <- as.data.frame(t(as.matrix(data.train)))
    data.test <- as.data.frame(t(as.matrix(data.test)))
    
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
    
    norm_data <- preprocessCore::normalize.quantiles(as.matrix(data.train))
    norm_data <- data.frame(norm_data, row.names = rownames(data.train))
    colnames(norm_data) <- colnames(data.train)
    data.train <- norm_data
    
    norm_data <- preprocessCore::normalize.quantiles(as.matrix(data.test))
    norm_data <- data.frame(norm_data, row.names = rownames(data.test))
    colnames(norm_data) <- colnames(data.test)
    data.test <- norm_data
  }
  #now data, data.test2 format : (samples x transcripts)
  
  #get best biomarkers only
  data.train <- data.frame(data.train)[, biomarkers]  #data.frame() replaces - in colnames to .
  
  colnames(data.test) <- gsub("-", ".", colnames(data.test))
  
  sum(biomarkers %in% colnames(data.test))
  
  available_biomarkers <- c(biomarkers[biomarkers %in% colnames(data.test)])
  print(available_biomarkers)
  non_available_biomarkers <- biomarkers[!biomarkers %in% colnames(data.test)] 
  print(non_available_biomarkers)
  
  data.test <- data.test %>%
    select(available_biomarkers)
  for(nab in non_available_biomarkers){
    data.test[[nab]] <- 0
  }
  
  
  # data.test2 = NA
  # label.test2 = NA
  # regularize = 'l2'
  if(omics_type == "transcriptomic"){
    result_df <- logistic_regression(data.train, label.train, data.test, label.test, 
                                     classes = classes, 
                                     regularize = 'l2')    
  } else if(omics_type == "proteomic"){
    #TODO : remove test2
    result_df <- svm_model(data.train, label.train, data.test, label.test, 
                           data.test2, label.test2,
                           classes, kernel = "sigmoid")
  }
  write.csv(format(result_df, digits = 3), result_file_name, row.names = FALSE)
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
