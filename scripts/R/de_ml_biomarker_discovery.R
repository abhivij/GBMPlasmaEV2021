library(tidyverse)
library(edgeR)
library(tidyverse)
library(sva)
library(readxl)
library(ranger)

source('scripts/R/prediction_pipeline/cm_rf.R')

#create 30 train-test splits 
#in each split, use train data for de analysis
#get logFC, p value for all
#create a list of siginificant tra/prot in each iter, append biomarkers identified to this list, append prev study markers to this list
#for each element of this list, get avg logFC and avg logP

#note - for some reason the functiond efinition execution stops due to some extra brackets issue
#so read all args, and execute the function code manually


#################################################################################################
data_file_path <- "Data/RNA_all/newquant_Nov2023_umi_counts_PREOPE_MET_HC_filter90.csv"
phenotype_file_path <- "Data/transcriptomic_phenotype_PREOPE_MET_HC_withaddicolumn.txt"
omics_type = "transcriptomics"
comparison = "PREOPEVsHC"
conditions = c("HC", "PREOPE")
pval_cutoff <- 0.05
best_features_file_path <- "Data/selected_features/best_features_with_add_col.csv"
dataset_replace_string <- "GBM_combined_transcriptomic_combat_compset2_new_quant_"
previous_study_biomarkers <- read_excel("Data/selected_features/features_of_interest.xlsx", sheet = "SerumExosome_GBMVsHC")
colnames(previous_study_biomarkers) <- "Molecule"
previous_study_biomarkers <- previous_study_biomarkers %>%
  mutate(source = "SerumExosome_GBMVsHC")
#note : doing it this way to support proteomics biomarkers from 2 sources
prev_data_available = TRUE
result_file_path = "DE_results_2024/tra_result_PREOPEVsHC_agg.csv"

#################################################################################################
data_file_path <- "Data/RNA_all/newquant_Nov2023_umi_counts_PREOPE_MET_HC_filter90.csv"
phenotype_file_path <- "Data/transcriptomic_phenotype_PREOPE_MET_HC_withaddicolumn.txt"
omics_type = "transcriptomics"
comparison = "PREOPEVsMET"
conditions = c("MET", "PREOPE")
pval_cutoff <- 0.05
best_features_file_path <- "Data/selected_features/best_features_with_add_col.csv"
dataset_replace_string <- "GBM_combined_transcriptomic_combat_compset2_new_quant_"
prev_data_available = FALSE
previous_study_biomarkers = NA
result_file_path = "DE_results_2024/tra_result_PREOPEVsMET_agg.csv"

#################################################################################################
data_file_path <- "Data/RNA_all/newquant_Nov2023_umi_counts_PREOPE_MET_HC_filter90.csv"
phenotype_file_path <- "Data/transcriptomic_phenotype_PREOPE_MET_HC_withaddicolumn.txt"
omics_type = "transcriptomics"
comparison = "METVsHC"
conditions = c("HC", "MET")
pval_cutoff <- 0.05
best_features_file_path <- "Data/selected_features/best_features_with_add_col.csv"
dataset_replace_string <- "GBM_combined_transcriptomic_new_quant_compset2_"
prev_data_available = FALSE
previous_study_biomarkers = NA
result_file_path = "DE_results_2024/tra_result_METVsHC_agg.csv"



#################################################################################################
data_file_path <- "Data/Protein/formatted_data/PREOPE_MET_HC_data.csv"
phenotype_file_path <- "Data/proteomic_phenotype_PREOPE_MET_HC_withaddicolumn.txt"
omics_type = "proteomics"
comparison = "PREOPEVsHC"
conditions = c("HC", "PREOPE")
pval_cutoff <- 0.05
best_features_file_path <- "Data/selected_features/best_features_with_add_col.csv"
dataset_replace_string <- "GBM_combined_proteomic_combat_compset2_"  

previous_study_biomarkers <- read_excel("Data/selected_features/features_of_interest.xlsx", sheet = "PlasmaEV_GBMVsHC")
colnames(previous_study_biomarkers) <- "Molecule"
previous_study_biomarkers <- previous_study_biomarkers %>%
  mutate(source1 = "PlasmaEV_GBMVsHC")

previous_study_biomarkers2 <- read_excel("Data/selected_features/features_of_interest.xlsx", sheet = "UrineEV_GBMVsHC")
colnames(previous_study_biomarkers2) <- "Molecule"
previous_study_biomarkers2 <- previous_study_biomarkers2 %>%
  mutate(source2 = "UrineEV_GBMVsHC")

previous_study_biomarkers <- previous_study_biomarkers %>%
  full_join(previous_study_biomarkers2, by = "Molecule")

prev_data_available = TRUE
result_file_path = "DE_results_2024/prot_result_PREOPEVsHC_agg.csv"

#################################################################################################
data_file_path <- "Data/Protein/formatted_data/PREOPE_MET_HC_data.csv"
phenotype_file_path <- "Data/proteomic_phenotype_PREOPE_MET_HC_withaddicolumn.txt"
omics_type = "proteomics"
comparison = "PREOPEVsMET"
conditions = c("MET", "PREOPE")
pval_cutoff <- 0.05
best_features_file_path <- "Data/selected_features/best_features_with_add_col.csv"
dataset_replace_string <- "GBM_combined_proteomic_combat_compset2_"  
previous_study_biomarkers <- NA
prev_data_available = FALSE
result_file_path = "DE_results_2024/prot_result_PREOPEVsMET_agg.csv"

#################################################################################################

data_file_path <- "Data/Protein/formatted_data/PREOPE_MET_HC_data.csv"
phenotype_file_path <- "Data/proteomic_phenotype_PREOPE_MET_HC_withaddicolumn.txt"
omics_type = "proteomics"
comparison = "METVsHC"
conditions = c("HC", "MET")
pval_cutoff <- 0.05
best_features_file_path <- "Data/selected_features/best_features_with_add_col.csv"
dataset_replace_string <- "GBM_combined_proteomic_combat_compset2_"  
previous_study_biomarkers <- NA
prev_data_available = FALSE
result_file_path = "DE_results_2024/prot_result_METVsHC_agg.csv"

#################################################################################################

de_ml_bimarker_discovery.compute_feature_rank <- function(data_file_path, phenotype_file_path, omics_type,
                                                          comparison, conditions, pval_cutoff,
                                                          best_features_file_path, dataset_replace_string,
                                                          prev_data_available = FALSE,
                                                          previous_study_biomarkers = NA,
                                                          result_file_path){
  data <- read.csv(data_file_path, row.names = 1)
  phenotype <- read.table(phenotype_file_path, header=TRUE, sep="\t")
  
  label <- phenotype %>%
    dplyr::rename("Label" = comparison) %>%
    dplyr::select(all_of(c("Sample", "Label", "data_cohort", "PREOPE_MET_HC"))) %>%
    mutate(condition_cohort = paste(PREOPE_MET_HC, data_cohort, sep = "_"))
  
  data <- data[, label$Sample]
  all.equal(label$Sample, colnames(data))
  
  #############################################
  
  best_features <- read.csv(best_features_file_path)
  best_features <- best_features %>%
    mutate(dataset_id = gsub(dataset_replace_string, "", dataset_id)) %>%
    filter(is_best == 1, dataset_id == comparison)
  best_features <- strsplit(best_features$biomarkers, split = "|", fixed = TRUE)[[1]]
  
  best_features <- gsub(".", "-", best_features, fixed = TRUE)
  
  #############################################
  
  if(omics_type == "transcriptomics"){
    data <- preprocess_data(data, label, perform_filter = TRUE, norm = "log_cpm")
  } else if(omics_type == "proteomics"){
    data <- preprocess_data(data, label, perform_filter = FALSE, norm = "quantile")
  } else{
    print("invalid omics_type")
    return
  }
  data.combat <- ComBat(dat = data, batch = label$data_cohort)
  data <- as.data.frame(data.combat)
  
  #############################################
  
  set.seed(1000)
  k_val <- 5
  times_val <- 6
  k_prod_times <- k_val * times_val
  train_index <- caret::createMultiFolds(y = label$condition_cohort, k = k_val, times = times_val)
  
  for (i in c(1:k_prod_times)) {
    
    # i <- 1
    label.train <- label[train_index[[i]], , drop = FALSE]
    label.test <- label[-train_index[[i]], , drop = FALSE]
    # print(head(label.train$Sample))
    # print(head(label.test$Sample))
    
    # print(summary(factor(label.train$PREOPE_MET_HC)))
    # print(summary(factor(label.test$PREOPE_MET_HC)))
    
    data.train <- data[, label.train$Sample]
    data.test <- data[, label.test$Sample]
    
    result.de <- perform_de(data, label, conditions) %>%
      mutate("iter" = i) %>%
      mutate("significant" = case_when((PVal < pval_cutoff) ~ 1,
                                       TRUE ~ 0))
    
    ##########################################
    label.train <- label.train %>%
      filter(!is.na(Label))
    label.test <- label.test %>%
      filter(!is.na(Label))
    data.train <- data[, label.train$Sample]
    data.test <- data[, label.test$Sample]
    
    result_list <- rf_model(as.data.frame(t(data.train)), label.train,
                            as.data.frame(t(data.test)), label.test,
                            conditions, classifier_feature_imp = TRUE)
    feature_importance <- result_list[[2]]
    feature_importance <- feature_importance %>%
      arrange(desc(MeanDecreaseGini)) %>%
      dplyr::rename(c("Molecule" = "feature"))
    
    result <- result.de %>%
      inner_join(feature_importance)
    ##########################################
    
    #obtain features selected by ranger to count number of times feature is chosen - to be used as stability score
    
    ranger_model <- ranger::ranger(x = as.data.frame(t(data.train)),
                                   y = factor(label.train$Label), importance = "impurity_corrected")
    features <- which(ranger_model$variable.importance > 0)
    
    feature_chosen <- data.frame(Molecule = rownames(data.train))
    feature_chosen <- feature_chosen %>%
      mutate(ranger_chosen = ifelse(Molecule %in% names(features), 1, 0))
    
    result <- result %>%
      inner_join(feature_chosen)
    
    if(i == 1){
      result_all <- result
    } else{
      result_all <- rbind(result_all, result)
    }
    
  }
  
  result.agg <- result_all %>%
    group_by(Molecule) %>%
    summarize(sig_count = sum(significant),
              fischers_combined_p = -2*sum(log(PVal)),
              mean_logFC = mean(logFC),
              MeanDecreaseGini_mean = mean(MeanDecreaseGini),
              stability = sum(ranger_chosen)) %>%
    mutate(identified_biomarker = case_when(Molecule %in% best_features ~ 1,
                                            TRUE ~ 0))
  
  if(prev_data_available){
    result.agg <- result.agg %>%
      left_join(previous_study_biomarkers)
  }
  
  #################################################################
  
  #obtain ranking and combined_score
  result.agg <- result.agg %>%
    mutate(rank.fischers_combined_p = rank(fischers_combined_p)) %>%
    mutate(rank.abs_mean_logFC = rank(abs(mean_logFC))) %>%
    mutate(rank.MeanDecreaseGini_mean = rank(MeanDecreaseGini_mean)) %>%
    mutate(rank.stability = rank(stability)) %>%
    mutate(score.identified_biomarker = 2 * nrow(result.agg) * identified_biomarker) #twice the weightage for selected biomarkers
  
  result.agg <- result.agg %>%
    mutate(combined_score = rank.fischers_combined_p + rank.abs_mean_logFC +
             rank.MeanDecreaseGini_mean + rank.stability +
             score.identified_biomarker,
           .after = Molecule) %>%
    arrange(desc(combined_score))
  
  #################################################################
  
  if(omics_type == "proteomics"){
    all_protein_names <- read.csv("Data/Protein/formatted_data/all_protein_names.csv") %>% 
      dplyr::select(c(from_id, primary_gene_id))
    
    result.agg <- result.agg %>%
      left_join(all_protein_names, 
                by = c("Molecule" = "from_id")) %>%
      relocate(primary_gene_id, .before = "Molecule")
  }
  
  write.csv(result.agg, result_file_path, row.names = FALSE)
  
}



preprocess_data <- function(data, output_labels, perform_filter, norm){
  #data format : (features x samples)
  
  if(perform_filter){
    keep <- edgeR::filterByExpr(data, group = output_labels$PREOPE_MET_HC)
    data <- data[keep, ]
  }
  
  if(norm == "log_cpm"){
    data <- edgeR::cpm(data, log = TRUE)
  } else if(norm == "quantile"){
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
  }
}


perform_de <- function(data, label, conditions){
  model_matrix <- model.matrix(~ 0 + label$PREOPE_MET_HC)
  colnames(model_matrix) <- sub("label$PREOPE_MET_HC", "", colnames(model_matrix), fixed = TRUE)
  
  contr_matrix <- makeContrasts(contrasts = paste0(conditions[2], " - ", conditions[1]),
                                levels = colnames(model_matrix))
  fit <- lmFit(data, model_matrix)
  # head(coef(fit))
  fit <- contrasts.fit(fit, contr_matrix)
  # head(coef(fit))
  efit <- eBayes(fit)
  top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
    rownames_to_column("molecule")
  result <- top.table %>%
    dplyr::select(molecule, logFC, P.Value, adj.P.Val) %>%
    dplyr::rename(Molecule = molecule, adjPVal = adj.P.Val, PVal = P.Value) %>%
    arrange(logFC)
  
  return(result)
}
