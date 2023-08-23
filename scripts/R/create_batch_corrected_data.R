library(tidyverse)
library(sva)

comparison = "POSTOPE_TPVsREC_TP"
classes = c("POSTOPE_TP", "REC_TP")
omics_type = "proteomics"
perform_filter = TRUE
norm = "quantile_train_param"
batch_effect_correction = "combat"

comparison = "PREOPEVsREC_TP"
classes = c("PREOPE", "REC_TP")
omics_type = "proteomics"
norm = "quantile_train_param"
filter_before_common = TRUE

create_batch_corrected_data <- function(comparison, classes, omics_type,
                                        norm,
                                        perform_filter = TRUE, filter_before_common = FALSE,
                                        batch_effect_correction = "combat", 
                                        include_unknown = FALSE){
  if(omics_type == "proteomics"){
    data_file_path <- "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv"
    phenotype_file_path <- "Data/proteomic_phenotype.txt"
    validation_data_file_path <- "Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil.csv"
    validation_metadata_file_path <- "Data/RNA_validation/metadata_glionet.csv"
    output_dir_path <- "Data/Protein"
    
  } else if(omics_type == "transcriptomics"){
    data_file_path <- "Data/RNA/umi_counts_initial_cohort.csv"
    phenotype_file_path <- "Data/transcriptomic_phenotype.txt"
    validation_data_file_path <- "Data/RNA/umi_counts_validation_cohort.csv"      
    validation_metadata_file_path <- "Data/RNA_validation/metadata_glionet.csv"
    output_dir_path <- "Data/RNA"
    
  } else{
    return("Invalid omics_type")
  }
  
  data <- read.csv(data_file_path, row.names = 1)
  phenotype <- read.table(phenotype_file_path, header=TRUE, sep="\t")
  
  validation_data <- read.csv(validation_data_file_path, row.names = 1)
  validation_metadata <- read.csv(validation_metadata_file_path) %>%
    rename("Label" = "category_old_name", "Sample" = "sample_id") %>%
    select(Sample, age, gender, Label)
  
  if(omics_type == "proteomics"){
    colnames(validation_data)[colnames(validation_data) == "SB12_01"] = "SB12"
    #use SB22.02
    colnames(validation_data)[colnames(validation_data) == "SB22.02"] = "SBtobeused22"
    colnames(validation_data)[colnames(validation_data) == "SB22"] = "SB22_dont_include"
    colnames(validation_data)[colnames(validation_data) == "SBtobeused22"] = "SB22"
    
    validation_metadata <- validation_metadata %>%
      filter(Sample != "SB7")    
  } else if(omics_type == "transcriptomics"){
    colnames(validation_data) <- paste0("S", colnames(validation_data))
  }
  
  phenotype <- phenotype %>%
    mutate(cohort = "initial")
  validation_metadata <- validation_metadata %>%
    mutate(cohort = "validation")
  
  output_labels.train <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes) %>%
    dplyr::select(Sample, Label, cohort)
  
  if(!include_unknown){
    output_labels.test <- validation_metadata %>%
      filter(Label %in% classes) %>%
      dplyr::select(Sample, Label, cohort)    
    file_name_prefix = ""
  } else{
    output_labels.test <- validation_metadata %>%
      filter(Label %in% classes | Sample %in% c('SB10', 'SB14', 'SB15', 'SB18', 'SB19', 
                                                'SB23', 'SB26', 'SB27', 'SB30', 'SB34', 
                                                'SB41', 'SB44', 'SB47', 'SB51', 'SB54', 
                                                'SB55', 'SB6')) %>%
      dplyr::select(Sample, Label, cohort)
    file_name_prefix = "include_unknown"
  }
  
  #currently data format : (transcripts x samples)
  
  data.train <- data %>% dplyr::select(output_labels.train$Sample)
  data.test <- validation_data %>% dplyr::select(output_labels.test$Sample)
  
  if(filter_before_common){
    if(perform_filter){
      keep <- edgeR::filterByExpr(data.train, group = output_labels.train$Label)
      data.train <- data.train[keep, ]
      keep <- edgeR::filterByExpr(data.test, group = output_labels.test$Label)
      data.test <- data.test[keep, ]  
    }    
    common <- intersect(rownames(data.train), rownames(data.test))  
    data.train <- data.train[common, ]
    data.test <- data.test[common, ]
    
    file_name_prefix = paste0(file_name_prefix, "filter_before.")
  } else{
    common <- intersect(rownames(data.train), rownames(data.test))  
    data.train <- data.train[common, ]
    data.test <- data.test[common, ]
    if(perform_filter){
      keep <- edgeR::filterByExpr(data.train, group = output_labels.train$Label)
      data.train <- data.train[keep, ]
      data.test <- data.test[keep, ]  
    }        
  }
  
  if(norm == "quantile_train_param"){
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
  } else if(norm == "log_cpm"){
    #calculating norm log cpm
    data.train <- edgeR::cpm(data.train, log=TRUE)
    data.test <- edgeR::cpm(data.test, log=TRUE)
  }
  data.train <- as.data.frame(t(as.matrix(data.train)))
  data.test <- as.data.frame(t(as.matrix(data.test)))  
  
  data_of_interest <- rbind(data.train, data.test)
  output_labels <- rbind(output_labels.train, output_labels.test)
  
  all.equal(rownames(data_of_interest), output_labels$Sample)
  data <- data_of_interest
  
  if(batch_effect_correction == "combat"){
    data <- as.data.frame(t(as.matrix(data)))
    data.combat = ComBat(dat=data, batch=output_labels$cohort)
    data.combat <- as.data.frame(data.combat)
  } else if(batch_effect_correction == "combat_with_mod"){
    data <- as.data.frame(t(as.matrix(data)))
    mod <- model.matrix(~as.factor(output_labels$Label), data = data)
    data.combat = ComBat(dat=data, batch=output_labels$cohort, mod=mod)
    data.combat <- as.data.frame(data.combat)
    file_name_prefix <- paste0(file_name_prefix, "mod_")
  } else if(batch_effect_correction == "combat_with_mod_testsame"){
    data <- as.data.frame(t(as.matrix(data)))
    output_labels.modified <- output_labels %>%
      mutate(Label = ifelse(cohort == "validation", "unk", Label))
    mod <- model.matrix(~as.factor(output_labels.modified$Label), data = data)
    data.combat = ComBat(dat=data, batch=output_labels.modified$cohort, mod=mod)
    data.combat <- as.data.frame(data.combat)    
    file_name_prefix <- paste0(file_name_prefix, "mod_testsame_")
  }
  data.train <- data.combat[, output_labels.train$Sample]
  data.test <- data.combat[, output_labels.test$Sample]

  print(dim(data.combat))
  print(dim(data.train))
  print(dim(data.test))
  
  write.csv(data.combat, file = paste0(output_dir_path, "/combined_data.combat.", 
                                      file_name_prefix,
                                      comparison, ".csv"))
  write.csv(data.train, file = paste0(output_dir_path, "/initial_data.combat.", 
                                      file_name_prefix,
                                      comparison, ".csv"))
  write.csv(data.test, file = paste0(output_dir_path, "/validation_data.combat.", 
                                     file_name_prefix,
                                     comparison, ".csv"))
}


create_batch_corrected_data(comparison = "POSTOPE_TPVsREC_TP",
                            classes = c("POSTOPE_TP", "REC_TP"),
                            omics_type = "proteomics",
                            norm = "quantile_train_param")
create_batch_corrected_data(comparison = "PREOPEVsPOSTOPE_TP",
                            classes = c("PREOPE", "POSTOPE_TP"),
                            omics_type = "proteomics",
                            norm = "quantile_train_param")
create_batch_corrected_data(comparison = "PREOPEVsREC_TP",
                            classes = c("PREOPE", "REC_TP"),
                            omics_type = "proteomics",
                            norm = "quantile_train_param")

# d1 <- read.csv("Data/Protein/initial_data.combat.PREOPEVsREC_TP.csv", row.names = 1)
# d2 <- read.csv("Data/Protein/validation_data.combat.PREOPEVsREC_TP.csv", row.names = 1)

create_batch_corrected_data(comparison = "POSTOPE_TPVsREC_TP",
                            classes = c("POSTOPE_TP", "REC_TP"),
                            omics_type = "transcriptomics",
                            norm = "log_cpm")
create_batch_corrected_data(comparison = "PREOPEVsPOSTOPE_TP",
                            classes = c("PREOPE", "POSTOPE_TP"),
                            omics_type = "transcriptomics",
                            norm = "log_cpm")
create_batch_corrected_data(comparison = "PREOPEVsREC_TP",
                            classes = c("PREOPE", "REC_TP"),
                            omics_type = "transcriptomics",
                            norm = "log_cpm")

# d1 <- read.csv("Data/RNA/initial_data.combat.PREOPEVsPOSTOPE_TP.csv", row.names = 1)
# d2 <- read.csv("Data/RNA/validation_data.combat.PREOPEVsPOSTOPE_TP.csv", row.names = 1)

#note : though files with filter_before_common = TRUE were created as shown below,
#it has not been run in the pipeline because proteins were similar or lesser in count 
#compared to previous
#didn't compare transcriptomics data

create_batch_corrected_data(comparison = "POSTOPE_TPVsREC_TP",
                            classes = c("POSTOPE_TP", "REC_TP"),
                            omics_type = "proteomics",
                            norm = "quantile_train_param",
                            filter_before_common = TRUE)
create_batch_corrected_data(comparison = "PREOPEVsPOSTOPE_TP",
                            classes = c("PREOPE", "POSTOPE_TP"),
                            omics_type = "proteomics",
                            norm = "quantile_train_param",
                            filter_before_common = TRUE)
create_batch_corrected_data(comparison = "PREOPEVsREC_TP",
                            classes = c("PREOPE", "REC_TP"),
                            omics_type = "proteomics",
                            norm = "quantile_train_param",
                            filter_before_common = TRUE)

create_batch_corrected_data(comparison = "POSTOPE_TPVsREC_TP",
                            classes = c("POSTOPE_TP", "REC_TP"),
                            omics_type = "transcriptomics",
                            norm = "log_cpm",
                            filter_before_common = TRUE)
create_batch_corrected_data(comparison = "PREOPEVsPOSTOPE_TP",
                            classes = c("PREOPE", "POSTOPE_TP"),
                            omics_type = "transcriptomics",
                            norm = "log_cpm",
                            filter_before_common = TRUE)
create_batch_corrected_data(comparison = "PREOPEVsREC_TP",
                            classes = c("PREOPE", "REC_TP"),
                            omics_type = "transcriptomics",
                            norm = "log_cpm",
                            filter_before_common = TRUE)

#################################
create_batch_corrected_data(comparison = "POSTOPE_TPVsREC_TP",
                            classes = c("POSTOPE_TP", "REC_TP"),
                            omics_type = "proteomics",
                            norm = "quantile_train_param",
                            batch_effect_correction = "combat_with_mod")
create_batch_corrected_data(comparison = "PREOPEVsPOSTOPE_TP",
                            classes = c("PREOPE", "POSTOPE_TP"),
                            omics_type = "proteomics",
                            norm = "quantile_train_param",
                            batch_effect_correction = "combat_with_mod")
create_batch_corrected_data(comparison = "PREOPEVsREC_TP",
                            classes = c("PREOPE", "REC_TP"),
                            omics_type = "proteomics",
                            norm = "quantile_train_param",
                            batch_effect_correction = "combat_with_mod")

create_batch_corrected_data(comparison = "POSTOPE_TPVsREC_TP",
                            classes = c("POSTOPE_TP", "REC_TP"),
                            omics_type = "transcriptomics",
                            norm = "log_cpm",
                            batch_effect_correction = "combat_with_mod")
create_batch_corrected_data(comparison = "PREOPEVsPOSTOPE_TP",
                            classes = c("PREOPE", "POSTOPE_TP"),
                            omics_type = "transcriptomics",
                            norm = "log_cpm",
                            batch_effect_correction = "combat_with_mod")
create_batch_corrected_data(comparison = "PREOPEVsREC_TP",
                            classes = c("PREOPE", "REC_TP"),
                            omics_type = "transcriptomics",
                            norm = "log_cpm",
                            batch_effect_correction = "combat_with_mod")


############################################################

#can't run the below ones
# gives
# Found2batches
# Adjusting for2covariate(s) or covariate level(s)
# Error in ComBat(dat = data, batch = output_labels.modified$cohort, mod = mod) : 
#   At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat



# create_batch_corrected_data(comparison = "POSTOPE_TPVsREC_TP",
#                             classes = c("POSTOPE_TP", "REC_TP"),
#                             omics_type = "proteomics",
#                             norm = "quantile_train_param",
#                             batch_effect_correction = "combat_with_mod_testsame")
# create_batch_corrected_data(comparison = "PREOPEVsPOSTOPE_TP",
#                             classes = c("PREOPE", "POSTOPE_TP"),
#                             omics_type = "proteomics",
#                             norm = "quantile_train_param",
#                             batch_effect_correction = "combat_with_mod_testsame")
# create_batch_corrected_data(comparison = "PREOPEVsREC_TP",
#                             classes = c("PREOPE", "REC_TP"),
#                             omics_type = "proteomics",
#                             norm = "quantile_train_param",
#                             batch_effect_correction = "combat_with_mod_testsame")
# 
# create_batch_corrected_data(comparison = "POSTOPE_TPVsREC_TP",
#                             classes = c("POSTOPE_TP", "REC_TP"),
#                             omics_type = "transcriptomics",
#                             norm = "log_cpm",
#                             batch_effect_correction = "combat_with_mod_testsame")
# create_batch_corrected_data(comparison = "PREOPEVsPOSTOPE_TP",
#                             classes = c("PREOPE", "POSTOPE_TP"),
#                             omics_type = "transcriptomics",
#                             norm = "log_cpm",
#                             batch_effect_correction = "combat_with_mod_testsame")
# create_batch_corrected_data(comparison = "PREOPEVsREC_TP",
#                             classes = c("PREOPE", "REC_TP"),
#                             omics_type = "transcriptomics",
#                             norm = "log_cpm",
#                             batch_effect_correction = "combat_with_mod_testsame")


#################

create_batch_corrected_data(comparison = "POSTOPE_TPVsREC_TP",
                            classes = c("POSTOPE_TP", "REC_TP"),
                            omics_type = "proteomics",
                            norm = "quantile_train_param", include_unknown = TRUE)
create_batch_corrected_data(comparison = "PREOPEVsPOSTOPE_TP",
                            classes = c("PREOPE", "POSTOPE_TP"),
                            omics_type = "proteomics",
                            norm = "quantile_train_param", include_unknown = TRUE)
create_batch_corrected_data(comparison = "PREOPEVsREC_TP",
                            classes = c("PREOPE", "REC_TP"),
                            omics_type = "proteomics",
                            norm = "quantile_train_param", include_unknown = TRUE)

create_batch_corrected_data(comparison = "POSTOPE_TPVsREC_TP",
                            classes = c("POSTOPE_TP", "REC_TP"),
                            omics_type = "transcriptomics",
                            norm = "log_cpm", include_unknown = TRUE)
create_batch_corrected_data(comparison = "PREOPEVsPOSTOPE_TP",
                            classes = c("PREOPE", "POSTOPE_TP"),
                            omics_type = "transcriptomics",
                            norm = "log_cpm", include_unknown = TRUE)
create_batch_corrected_data(comparison = "PREOPEVsREC_TP",
                            classes = c("PREOPE", "REC_TP"),
                            omics_type = "transcriptomics",
                            norm = "log_cpm", include_unknown = TRUE)



####################

comparison = "PREOPEVsMET"
classes = c("MET", "PREOPE")
omics_type = "proteomics"
norm = "quantile_train_param"
perform_filter = FALSE
batch_effect_correction = "combat"

create_batch_corrected_data_PMH <- function(comparison, classes, omics_type,
                                            norm,
                                            perform_filter = TRUE,
                                            batch_effect_correction = "combat"){
  if(omics_type == "proteomics"){
    data_file_path <- "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv"
    validation_data_file_path <- "Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil.csv" 
    phenotype_file_path <- "Data/proteomic_phenotype_PREOPE_MET_HC.txt"
    output_dir_path <- "Data/Protein"
  } else if(omics_type == "transcriptomics"){
    data_file_path <- "Data/RNA/umi_counts_initial_cohort.csv"
    validation_data_file_path <- "Data/RNA/umi_counts_validation_cohort.csv"      
    phenotype_file_path <- "Data/transcriptomic_phenotype_PREOPE_MET_HC.txt"
    output_dir_path <- "Data/RNA"
  } else{
    return("Invalid omics_type")
  }
  
  data <- read.csv(data_file_path, row.names = 1)
  validation_data <- read.csv(validation_data_file_path, row.names = 1)
  
  phenotype <- read.table(phenotype_file_path, header=TRUE, sep="\t")
  
  if(omics_type == "proteomics"){
    colnames(validation_data)[colnames(validation_data) == "SB12_01"] = "SB12"
    #use SB22.02
    colnames(validation_data)[colnames(validation_data) == "SB22.02"] = "SBtobeused22"
    colnames(validation_data)[colnames(validation_data) == "SB22"] = "SB22_dont_include"
    colnames(validation_data)[colnames(validation_data) == "SBtobeused22"] = "SB22"
   
  } else if(omics_type == "transcriptomics"){
    colnames(validation_data) <- paste0("S", colnames(validation_data))
  }
  
  output_labels.cohort1 <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes) %>%
    dplyr::select(Sample, Label, data_cohort, Subgroup, Sex, Age) %>%
    filter(data_cohort == "initial")
  output_labels.cohort2 <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes) %>%
    dplyr::select(Sample, Label, data_cohort, Subgroup, Sex, Age) %>%
    filter(data_cohort == "validation")
  
  #currently data format : (transcripts x samples)
  
  data.cohort1 <- data %>% dplyr::select(output_labels.cohort1$Sample)
  data.cohort2 <- validation_data %>% dplyr::select(output_labels.cohort2$Sample)
  
  #take common only if both cohorts contain samples
  if(ncol(data.cohort1) > 0 & ncol(data.cohort2) > 0){
    common <- intersect(rownames(data.cohort1), rownames(data.cohort2))  
    data.cohort1 <- data.cohort1[common, ]
    data.cohort2 <- data.cohort2[common, ]
    data <- cbind(data.cohort1, data.cohort2)    
  } else if(ncol(data.cohort1) > 0){
    data <- data.cohort1
  } else if(ncol(data.cohort2) > 0){
    data <- data.cohort2
  } else{
    print("no samples")
  }
  
  output_labels <- rbind(output_labels.cohort1, output_labels.cohort2)

  if(perform_filter){
    keep <- edgeR::filterByExpr(data, group = output_labels$Label)
    data_left_out <- data[!keep, ]
    data <- data[keep, ]
  }
  
  if(norm == "quantile_train_param"){
    #adapted from https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
    data.rank <- apply(data, 2, rank, ties.method="average")
    data.sorted <- data.frame(apply(data, 2, sort))
    data.mean <- apply(data.sorted, 1, mean)
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
    data.norm <- apply(data.rank, 2, index_to_mean, data_mean = data.mean)
    rownames(data.norm) <- rownames(data)
    data <- data.norm
    
  } else if(norm == "log_cpm"){
    #calculating norm log cpm
    data <- edgeR::cpm(data, log=TRUE)
  }
  
  data <- as.data.frame(t(as.matrix(data)))  
  
  all.equal(rownames(data), output_labels$Sample)
  
  if(batch_effect_correction == "combat"){
    data <- as.data.frame(t(as.matrix(data)))
    data.bec = ComBat(dat=data, batch=output_labels$data_cohort)
    data.bec <- as.data.frame(data.bec)
  } else if(batch_effect_correction == "combat_with_mod"){
    data <- as.data.frame(t(as.matrix(data)))
    mod <- model.matrix(~as.factor(output_labels$Label), data = data)
    data.bec = ComBat(dat=data, batch=output_labels$data_cohort, mod=mod)
    data.bec <- as.data.frame(data.bec)
    file_name_prefix <- paste0(file_name_prefix, "mod_")
  } else{
    data.bec <- as.data.frame(t(as.matrix(data)))
  }
  
  print(dim(data.bec))
  
  write.csv(data.bec, file = paste0(output_dir_path, "/combined_data.", 
                                    batch_effect_correction, ".", 
                                    comparison, ".csv"))
}

create_batch_corrected_data_PMH(comparison = "PREOPEVsMET",
                                classes = c("MET", "PREOPE"),
                                omics_type = "proteomics",
                                norm = "quantile_train_param",
                                perform_filter = FALSE,
                                batch_effect_correction = "combat")
create_batch_corrected_data_PMH(comparison = "PREOPEVsHC",
                                classes = c("HC", "PREOPE"),
                                omics_type = "proteomics",
                                norm = "quantile_train_param",
                                perform_filter = FALSE,
                                batch_effect_correction = "combat")
#running the below to create normalized data file for Agota/Susannah
create_batch_corrected_data_PMH(comparison = "METVsHC",
                                classes = c("HC", "MET"),
                                omics_type = "proteomics",
                                norm = "quantile_train_param",
                                perform_filter = FALSE,
                                batch_effect_correction = "")

create_batch_corrected_data_PMH(comparison = "PREOPEVsMET",
                                classes = c("MET", "PREOPE"),
                                omics_type = "transcriptomics",
                                norm = "log_cpm",
                                perform_filter = TRUE,
                                batch_effect_correction = "combat")
create_batch_corrected_data_PMH(comparison = "PREOPEVsHC",
                                classes = c("HC", "PREOPE"),
                                omics_type = "transcriptomics",
                                norm = "log_cpm",
                                perform_filter = TRUE,
                                batch_effect_correction = "combat")
#running the below to create normalized data file for Agota/Susannah
create_batch_corrected_data_PMH(comparison = "METVsHC",
                                classes = c("HC", "MET"),
                                omics_type = "transcriptomics",
                                norm = "log_cpm",
                                perform_filter = TRUE,
                                batch_effect_correction = "")
