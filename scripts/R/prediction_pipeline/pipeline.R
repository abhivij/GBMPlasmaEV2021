# source("R/data_extraction/extract.R")
# source("R/feature_extraction_arguments.R")
# source("R/run_fsm_and_models.R")
# source("R/helper.R")

#' FEM pipeline
#' 
#' Runs the Feature Extraction Method Comparsion pipeline for the given data
#' @param phenotype_file_name Name of the file containing phenotype info - class of each sample.
#' 
#' - should be a tab separated file 
#' 
#' - should have one field named "Sample" and contain sample names as in the read count file
#' 
#' 
#' example phenotype file contents given below :
#' 
#' Sample Age Gender Technology    GBMVsControl     GBMVsGlioma
#' 
#' s1     60    M    Microarray    GBM     GBM
#' 
#' s2     65    M   RNASeq        GBM     GBM
#' 
#' s3     59    F    RNASeq        Control NA
#' 
#' s4     55    F    Microarray    NA        Glioma
#' 
#' s5     50    M     RNASeq        Control NA
#' 
#' @param read_count_dir_path directory path for read count file
#' @param read_count_file_name filename of the read count file. Should be in (transcripts x samples) format.
#' Other biomarkers can also be used.
#' In general case : (biomarkers x samples) : biomarkers along rows, samples along columns
#'
#' @param sep field separator character in read count file
#' @param classification_criteria Column in the phenotype file to perform classification on
#' Eg : GBMVsControl 
#' @param filter_expression Filtering to be done on the samples based on a column in the phenotype file 
#' Eg: Age > 55
#' @param classes Classes to be compared : c("negativeclassname", "positiveclassname")
#' @param extracted_count_file_name name of the file to output the read count file after filtering for classification_criteria and filter_expression.
#' File data format : (transcripts x samples)
#' @param output_label_file_name name of the file with labels for filtered samples as in extracted_count_file.
#' File data format : (2 columns : Sample, Label)
#' @param dataset_id An ID for the data to be written in results
#' @param cores Number of cores to be used for parallel processing
#' @param fems_to_run Vector of names of FEMs to run  Eg: c("t-test", "mrmr30", "mrmr50"). Empty vector value runs pipeline on all FEMs
#' @param fems_to_ignore Vector of names of FEMs to not run from the list of all allowed FEMs  Eg: c("t-test_holm", "ga_rf")
#' 
#'  To see all allowed FEMs use show_allowed_fems()
#' 
#' @param perform_filter Should filtering for low expressed transcripts be performed
#' @param norm Normalization method to be used
#' @param classifier_feature_imp Should feature importance by the classifier be written to a file. 
#' Applicable only for Random Forest
#' @param random_seed random seed value set before train test k-fold split
#' @export 
execute_pipeline <- function(phenotype_file_name, 
                             read_count_dir_path, 
                             read_count_file_name,
                             sep = "",
                             skip_row_count = 0, 
                             row_count = -1,
                             na_strings = "NA",
                             classification_criteria, 
                             filter_expression = expression(TRUE), 
                             classes,
                             extracted_count_file_name = "read_counts.txt",
                             output_label_file_name = "output_labels.txt",
                             dataset_id, 
                             cores = 3,
                             results_dir_path = "results",
                             fems_to_run = c(),
                             fems_to_ignore = c(),
                             perform_filter = TRUE,
                             norm = c("norm_log_cpm", "norm_log_cpm_simple",
                                      "quantile", "norm_quantile", 
                                      "vsn", FALSE),
                             classifier_feature_imp = FALSE,
                             random_seed = 1000){
  start_time <- Sys.time()
  print(paste("Pipeline Execution on", dataset_id, classification_criteria))
  
  ######################################################################
  
  comparison = "POSTOPE_TPVsREC_TP"
  omics_type = "transcriptomic"
  omics_type = "proteomic"
  
  #conditions : c(pos_class, neg_class, validation_class)
  conditions = c("POSTOPE_TP", "REC_TP", "PREREC")
  classes = conditions[1:2]
  
  phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP"
  
  best_features_file_path = "Data/selected_features/best_features_with_add_col.csv"
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
  
  random_seed = 1000
  set.seed(random_seed)
  train_index <- caret::createDataPartition(output_labels$Label, p = .8, 
                                            list = FALSE, 
                                            times = 1)
  
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
    logistic_regression(data.train, label.train, data.test, label.test, 
                        data.test2, label.test2,
                        classes = classes, 
                        regularize = 'l2',
                        result_file_name = "Data/prediction_result/transcriptomics.csv")    
  } else if(omics_type == "proteomic"){
    svm_model(data.train, label.train, data.test, label.test, 
                          data.test2, label.test2,
                          classes, kernel = "sigmoid", 
                          result_file_name = "Data/prediction_result/proteomics.csv")
  }

  
  end_time <- Sys.time()
  print(end_time - start_time)
}
