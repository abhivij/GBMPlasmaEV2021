library(tidyverse)
library(readxl)
library(ROCR)

library(ggvenn)
library(ComplexHeatmap)


base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV/"
setwd(base_dir)

#create heatmap of average prediction probabilities

best_features <- read.csv("Data/selected_features/best_features_with_add_col.csv") %>%
  filter(is_best == 1, 
         grepl("_combined_", dataset_id),
         grepl('_common_combat_', dataset_id)) %>%
  mutate(comparison = gsub("GBM_combined_transcriptomic_common_combat_|GBM_combined_proteomic_common_combat_",
                           "", dataset_id)) %>%
  mutate(dataset_id = paste(dataset_id, description, 
                            min_iter_feature_presence, comparison,
                            sep = "_")) %>%
  dplyr::select(c(dataset_id, comparison))

sample_wise_results.proteomics <- read.csv("fem_pipeline_results_combined_pr_common_combat_subset/all_samplewise_result_df.csv") %>%
  inner_join(best_features %>% filter(grepl("proteomic", dataset_id)), 
             by = c("DataSetId" = "dataset_id")) %>%
  filter(FSM == "all") %>%
  select(-c(FSM))
sample_wise_results.transcriptomics <- read.csv("fem_pipeline_results_combined_tr_common_combat_subset/all_samplewise_result_df.csv") %>%
  inner_join(best_features %>% filter(grepl("transcriptomic", dataset_id)), 
             by = c("DataSetId" = "dataset_id")) %>%
  filter(FSM == "all") %>%
  select(-c(FSM))

sample_type <- "test"
comparison_of_interest <- "POSTOPE_TPVsREC_TP"
classes <- c("REC_TP", "POSTOPE_TP")
create_mean_prob_heatmap <- function(sample_wise_results.transcriptomics, 
                                     sample_wise_results.proteomics,
                                     comparison_of_interest,
                                     classes,
                                     sample_type){
  results.tra <- sample_wise_results.transcriptomics %>%
    filter(Type == sample_type, comparison == comparison_of_interest) %>%
    dplyr::select(-c(DataSetId, comparison, Type)) %>%
    group_by(Sample, TrueLabel, Model) %>%
    summarize(mean_prob = mean(PredProb))
  results.prot <- sample_wise_results.proteomics %>%
    filter(Type == sample_type, comparison == comparison_of_interest) %>%
    dplyr::select(-c(DataSetId, comparison, Type)) %>%
    group_by(Sample, TrueLabel, Model) %>%
    summarize(mean_prob = mean(PredProb))
  
  sample_label.tra <- results.tra %>%
    dplyr::select(c(Sample, TrueLabel)) %>%
    distinct() %>%
    arrange(Sample)
  sample_label.prot <- results.prot %>%
    dplyr::select(c(Sample, TrueLabel)) %>%
    distinct() %>%
    mutate(Sample = sub("HB0", "HB", Sample)) %>%
    arrange(Sample)
  
  setdiff(sample_label.prot$Sample, sample_label.tra$Sample)
  setdiff(sample_label.tra$Sample, sample_label.prot$Sample)
  
  results.prot <- results.prot %>%
    mutate(Sample = sub("HB0", "HB", Sample)) 
  
  results <- rbind(results.tra %>% mutate(omics_type = "transcriptomics"),
                   results.prot %>% mutate(omics_type = "proteomics")) %>%
    ungroup() %>%
    mutate(model_omics_type = paste(Model, omics_type, sep = "_"))
  
  data_to_plot <- results %>%
    dplyr::select(c(Sample, model_omics_type, mean_prob)) %>%
    pivot_wider(names_from = model_omics_type, values_from = mean_prob) %>%
    arrange(Sample) %>%
    column_to_rownames("Sample")
  data_to_plot <- data.matrix(data_to_plot)
  
  meta_data.row <- rbind(sample_label.tra, sample_label.prot) %>%
    distinct()
  meta_data.col <- data.frame(model_omics_type = colnames(data_to_plot)) %>%
    separate(model_omics_type, into = c("model", "omics_type"), sep = "_", remove = FALSE)
  
  Heatmap(data_to_plot, name = "Mean Prediction probability",
          col = viridis(5),
          rect_gp = gpar(col = "white", lwd = 1),
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          row_title = "Samples",
          row_names_side = "left",
          row_split = meta_data.row$TrueLabel,
          column_split = meta_data.col$omics_type,
          column_title = NULL,
          show_column_names = FALSE,
          bottom_annotation = HeatmapAnnotation(
            "Omics type" = meta_data.col$omics_type,
            "Model" = meta_data.col$model,
            annotation_name_side = "left"
          )) + 
    HeatmapAnnotation("TrueLabel" = meta_data.row$TrueLabel, 
                                 which = "row", 
                      col = list("TrueLabel" = c("POSTOPE_TP" = "#440154", "REC_TP" = "#fde725")))
  
  
  Heatmap(data_to_plot, name = "Mean Prediction probability",
          col = magma(5),
          rect_gp = gpar(col = "white", lwd = 1),
          cluster_columns = FALSE,
          column_order = replaced_dataset_id_vec[replaced_dataset_id_vec %in% colnames(data_to_plot)],
          column_names_rot = 60,
          column_names_gp = gpar(fontsize = 10),
          row_title = "Samples",
          row_names_side = "left",
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.3f", data_to_plot[i, j]), x, y, gp = gpar(fontsize = 7, col = "slateblue3"))
          })
  
  
}


sample_wise_results.train <- sample_wise_results %>%
  filter(Type == "train")
sample_wise_results.test <- sample_wise_results %>%
  filter(Type == "test")

current_comparison <- "CFRDVsNGT"
classes <- c("NGT", "CFRD")
for(i in c(1:30)){
  #i <- 1
  results_subset <- sample_wise_results.test %>%
    filter(Iter == i, comparison == current_comparison)
  
  pr <- ROCR::prediction(results_subset$PredProb, 
                         results_subset$TrueLabel, label.ordering = classes)
  #compute ROC curve, and AUC  
  prf <- ROCR::performance(pr, measure = "tpr", x.measure = "fpr")
  if(i == 1){
    plot(prf, col = i)  
  } else{
    plot(prf, col = i, add = TRUE)
  }
  
  print(ROCR::performance(pr, measure = "auc")@y.values[[1]])
}
