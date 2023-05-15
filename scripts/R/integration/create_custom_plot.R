library(tidyverse)
library(readxl)
library(ROCR)

library(ggvenn)
library(ComplexHeatmap)
library(viridis)
library(RColorBrewer)

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
  dplyr::select(-c(FSM))
sample_wise_results.transcriptomics <- read.csv("fem_pipeline_results_combined_tr_common_combat_subset/all_samplewise_result_df.csv") %>%
  inner_join(best_features %>% filter(grepl("transcriptomic", dataset_id)), 
             by = c("DataSetId" = "dataset_id")) %>%
  filter(FSM == "all") %>%
  dplyr::select(-c(FSM))
# 
sample_type <- "test"
comparison_of_interest <- "POSTOPE_TPVsREC_TP"
classes <- c("REC_TP", "POSTOPE_TP")




#classes assumed to be given as (neg_class, pos_class)
create_mean_prob_heatmap <- function(sample_wise_results.transcriptomics, 
                                     sample_wise_results.proteomics,
                                     comparison_of_interest,
                                     classes,
                                     sample_type,
                                     row_col){
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
    mutate("Sample_char_part" = str_extract(Sample, "[A-Z]+"),
           "Sample_num_part" = as.numeric(str_extract(Sample, "[0-9]+"))) %>%
    arrange(Sample_char_part, Sample_num_part) %>%
    dplyr::select(-c(Sample_char_part, Sample_num_part)) %>%
    column_to_rownames("Sample")
  data_to_plot <- data.matrix(data_to_plot)
  
  model_names <- c("Simple logistic regression",
                   "Radial Kernel SVM",
                   "Sigmoid Kernel SVM",
                   "Random Forest",
                   "L1 Regularized logistic regression",
                   "L2 Regularized logistic regression",
                   "Elastic net logistic regression")
  
  meta_data.row <- data.frame(Sample = rownames(data_to_plot)) %>% 
    inner_join(rbind(sample_label.tra, sample_label.prot) %>%
                 distinct()) %>%
    mutate(TrueLabel = factor(TrueLabel, levels = rev(classes)))
  meta_data.col <- data.frame(model_omics_type = colnames(data_to_plot)) %>%
    separate(model_omics_type, into = c("model", "omics_type"), sep = "_", remove = FALSE) %>%
    mutate(model = factor(model, levels = model_names))
  
  row_col <- list()
  row_col[["TrueLabel"]] <- c("#440154", "#fde725")
  names(row_col[["TrueLabel"]]) <- classes
  column_col <- list("Omics type" = c("proteomics" = "skyblue1",
                                      "transcriptomics" = "indianred1"))
  column_col[["Model"]] <- brewer.pal(n = 7, name = "Paired")
  names(column_col[["Model"]]) <- model_names
  
  ht <- Heatmap(data_to_plot, name = "Mean Prediction probability",
                col = viridis(5),
                rect_gp = gpar(col = "white", lwd = 1),
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                row_title = "Samples",
                row_names_side = "left",
                row_split = meta_data.row$TrueLabel,
                column_split = meta_data.col$model,
                column_title = NULL,
                show_column_names = FALSE,
                row_names_gp = gpar(fontsize = 5),
                bottom_annotation = HeatmapAnnotation(
                  "Omics type" = meta_data.col$omics_type,
                  "Model" = meta_data.col$model,
                  annotation_name_side = "left",
                  col = column_col
                )) + 
    HeatmapAnnotation("TrueLabel" = meta_data.row$TrueLabel, 
                      which = "row", 
                      col = row_col)
  
  png(paste0("plots/custom_heatmap/", comparison_of_interest, 
             sample_type,
             ".png"), units = "cm", width = 20, height = 15, res = 1200)  
  draw(ht, column_title = paste(sample_type, sub("Vs", " Vs ", comparison_of_interest)))    
  dev.off() 
  
}

create_mean_prob_heatmap(sample_wise_results.transcriptomics,
                         sample_wise_results.proteomics,
                         comparison_of_interest = "POSTOPE_TPVsREC_TP",
                         classes = c("REC_TP", "POSTOPE_TP"),
                         sample_type = "test")
create_mean_prob_heatmap(sample_wise_results.transcriptomics,
                         sample_wise_results.proteomics,
                         comparison_of_interest = "PREOPEVsPOSTOPE_TP",
                         classes = c("POSTOPE_TP", "PREOPE"),
                         sample_type = "test")
create_mean_prob_heatmap(sample_wise_results.transcriptomics,
                         sample_wise_results.proteomics,
                         comparison_of_interest = "PREOPEVsREC_TP",
                         classes = c("REC_TP", "PREOPE"),
                         sample_type = "test")

create_mean_prob_heatmap(sample_wise_results.transcriptomics,
                         sample_wise_results.proteomics,
                         comparison_of_interest = "POSTOPE_TPVsREC_TP",
                         classes = c("REC_TP", "POSTOPE_TP"),
                         sample_type = "train")
create_mean_prob_heatmap(sample_wise_results.transcriptomics,
                         sample_wise_results.proteomics,
                         comparison_of_interest = "PREOPEVsPOSTOPE_TP",
                         classes = c("POSTOPE_TP", "PREOPE"),
                         sample_type = "train")
create_mean_prob_heatmap(sample_wise_results.transcriptomics,
                         sample_wise_results.proteomics,
                         comparison_of_interest = "PREOPEVsREC_TP",
                         classes = c("REC_TP", "PREOPE"),
                         sample_type = "train")

# 
# sample_wise_results.train <- sample_wise_results %>%
#   filter(Type == "train")
# sample_wise_results.test <- sample_wise_results %>%
#   filter(Type == "test")
# 
# current_comparison <- "CFRDVsNGT"
# classes <- c("NGT", "CFRD")
# for(i in c(1:30)){
#   #i <- 1
#   results_subset <- sample_wise_results.test %>%
#     filter(Iter == i, comparison == current_comparison)
#   
#   pr <- ROCR::prediction(results_subset$PredProb, 
#                          results_subset$TrueLabel, label.ordering = classes)
#   #compute ROC curve, and AUC  
#   prf <- ROCR::performance(pr, measure = "tpr", x.measure = "fpr")
#   if(i == 1){
#     plot(prf, col = i)  
#   } else{
#     plot(prf, col = i, add = TRUE)
#   }
#   
#   print(ROCR::performance(pr, measure = "auc")@y.values[[1]])
# }
