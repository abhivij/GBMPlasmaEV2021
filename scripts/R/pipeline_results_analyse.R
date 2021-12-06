setwd("~/UNSW/VafaeeLab/GBMPlasmaEV")
library(tidyverse)
library(viridis)
library(ComplexHeatmap)
source("scripts/R/utils.R")
source("scripts/R/dataset_pipeline_arguments.R")

# data_info <- read.table('data_info.csv', sep = ',', header = TRUE)
# fsm_info <- read.table('fsm_info.csv', sep = ',', header = TRUE)
# model_results <- read.table('model_results.csv', sep = ',', header = TRUE)

# setwd("~/UNSW/VafaeeLab/GBMPlasmaEV/")

dataset_vec <- c("GBMPlasmaEV_transcriptomic_PREOPEVsMET",
                 "GBMPlasmaEV_transcriptomic_PREOPEVsHC",
                 "GBMPlasmaEV_transcriptomic_METVsHC",
                 "GBMPlasmaEV_proteomic_norm_quantile_PREOPEVsMET",
                 "GBMPlasmaEV_proteomic_norm_quantile_PREOPEVsHC",
                 "GBMPlasmaEV_proteomic_norm_quantile_METVsHC",
                 "GBMPlasmaEV_proteomic_impute75fil_norm_quantile_PREOPEVsMET",
                 "GBMPlasmaEV_proteomic_impute75fil_norm_quantile_PREOPEVsHC",
                 "GBMPlasmaEV_proteomic_impute75fil_norm_quantile_METVsHC",
                 "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsMET",
                 "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsHC",
                 "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_METVsHC"                 
                 )
dataset_vec <- dataset_vec[c(1:3, 10:12)]

fsm_vec <- c("all", 
  "t-test", "t-test_BH",
  "t-test_pval_0.025", "t-test_pval_0.01", "t-test_pval_0.005",
  "wilcoxontest", "wilcoxontest_BH",
  "wilcoxontest_pval_0.025", "wilcoxontest_pval_0.001", "wilcoxontest_pval_0.005",
  "ranger_impu_cor", 
  "mrmr10", "mrmr20",
  "mrmr30", "mrmr50", 
  "mrmr75", "mrmr100",
  "RF_RFE", "ga_rf")

# fsm_vec[c(1:2, 4, 6:12)]


data_info <- data_info %>%
  filter(DataSetId %in% dataset_vec)
fsm_info <- fsm_info %>%
  filter(DataSetId %in% dataset_vec) 
# %>%
#   filter(FSM %in% fsm_vec[c(1:2, 4, 6:12)])
model_results <- model_results %>%
  filter(DataSetId %in% dataset_vec)
# %>%
#   filter(FSM %in% fsm_vec[c(1:2, 4, 6:12)])

model_results <- model_results %>%
  mutate(FSM = factor(FSM)) %>%
  mutate(Model = factor(Model)) %>%
  mutate(DataSetId = factor(DataSetId, levels = dataset_vec))


dparg_vec = c(139:141, 145:147, 151, 153, 155, 157)
results_dir = "fem_pipeline_results"
dir_path = "plots/FEMPipeline/heatmap/transcriptomic_simple_norm_all/"
dataset_replace_string = "proteomic_impute50fil_"

plot_heatmap <- function(dparg_vec,
                         results_dir = "fem_pipeline_results",
                         dir_path = "",
                         dataset_replace_string = ""){
  data_info <- read.table(paste(results_dir, "data_info.csv", sep = "/"), 
                          sep = ',', header = TRUE)
  fsm_info <- read.table(paste(results_dir, "fsm_info.csv", sep = "/"),
                         sep = ',', header = TRUE)
  model_results <- read.table(paste(results_dir, "model_results.csv", sep = "/"),
                              sep = ',', header = TRUE)
  
  if(!dir.exists(dir_path)){
    dir.create(dir_path)
  }
  
  dataset_id_vec <- c()
  for(id in dparg_vec){
    dataset_id <- paste(dataset_pipeline_arguments[[id]]$dataset_id,
                        dataset_pipeline_arguments[[id]]$classification_criteria,
                        sep = "_")
    print(dataset_id)
    dataset_id_vec <- append(dataset_id_vec, dataset_id)
  }
  dataset_id_vec 
  
  replaced_dataset_id_vec <- gsub("GBMPlasmaEV_", "", dataset_id_vec, fixed = TRUE)
  replaced_dataset_id_vec <- gsub(dataset_replace_string, "", 
                                  replaced_dataset_id_vec, fixed = TRUE)
  
  model_results <- model_results %>%
    filter(DataSetId %in% dataset_id_vec) %>%
    mutate(DataSetId = gsub("GBMPlasmaEV_", "", DataSetId, fixed = TRUE)) %>%
    mutate(DataSetId = gsub(dataset_replace_string, "", DataSetId, fixed = TRUE))
  
  # model <- "L2 Regularized logistic regression"
  for(model in unique(model_results$Model)){
    print(model)
    
    model_results_sub <- model_results %>%
      filter(Model == model)
    
    data_to_plot <- model_results_sub %>%
      select(DataSetId, FSM, Mean_AUC) %>%
      # pivot_wider(names_from = DataSetId, values_from = Mean_AUC, values_fn = length) %>%          
      pivot_wider(names_from = DataSetId, values_from = Mean_AUC) %>%
      column_to_rownames("FSM")
    data_to_plot <- data.matrix(data_to_plot)
    
    file_name_curr <- paste0("Mean_AUC",
                               gsub(" ", "_", model),
                               ".png")
    file_name_curr <- paste0(dir_path, file_name_curr)

    png(file_name_curr, units = "cm", width = 20, height = 15, res = 1200)
    ht <- Heatmap(data_to_plot, name = "Mean AUC",
                  col = magma(5),
                  rect_gp = gpar(col = "white", lwd = 1),
                  cluster_columns = FALSE,
                  column_order = replaced_dataset_id_vec[replaced_dataset_id_vec %in% colnames(data_to_plot)],
                  column_names_rot = 60,
                  column_names_gp = gpar(fontsize = 10),
                  row_title = "Feature Selection Methods",
                  row_names_side = "left",
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(sprintf("%.3f", data_to_plot[i, j]), x, y, gp = gpar(fontsize = 7, col = "slateblue3"))
                  })
    draw(ht, column_title = paste("Mean AUC results from", model))
    dev.off()
    
  }
  
  dir_path <- paste0(dir_path, "common_barplot/")
  if(!dir.exists(dir_path)){
    dir.create(dir_path)
  }
  classification_models <- unique(model_results$Model)
  for (cm in classification_models) {
    individual_model <- model_results %>%
      filter(Model == cm)
    model_barplot <- ggplot(individual_model, aes(x=DataSetId, fill=FSM, y=Mean_AUC)) +
      geom_bar(stat="identity", position="dodge") +
      geom_errorbar( aes(x=DataSetId, ymin=X95.CI_AUC_lower, ymax=X95.CI_AUC_upper), position="dodge") +
      scale_fill_viridis_d() +
      theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
            axis.text.y = element_text(size=rel(1.2), face="italic", hjust=0.95),
            strip.text = element_text(size=rel(1.2), face="bold"),
            legend.title = element_text(size=rel(1.2)),
            legend.text = element_text(size=rel(1.1))
      ) +
      facet_wrap(facets = vars(Model))  
    
    file_name <- paste0(dir_path, str_replace(cm, " ", ""), "_barplot.png")
    ggsave(file_name, model_barplot, width=20, height=10, dpi=500)
  }  
}


plot_heatmap(dparg_vec = c(41:43, 47:49, 53:55),
             results_dir = "fem_pipeline_results",
             dir_path = "plots/FEMPipeline/heatmap/proteomic_norm_compare/",
             dataset_replace_string = "proteomic_impute50fil_")

plot_heatmap(dparg_vec = c(41:43, 125:127, 131, 133, 135, 137),
             results_dir = "fem_pipeline_results",
             dir_path = "plots/FEMPipeline/heatmap/proteomic_quantile/",
             dataset_replace_string = "proteomic_impute50fil_")

plot_heatmap(dparg_vec = c(139:141),
             results_dir = "fem_pipeline_results",
             dir_path = "plots/FEMPipeline/heatmap/transcriptomic_simple_norm/",
             dataset_replace_string = "transcriptomic_")

plot_heatmap(dparg_vec = c(139:141, 145:147, 151, 153, 155, 157),
             results_dir = "fem_pipeline_results",
             dir_path = "plots/FEMPipeline/heatmap/transcriptomic_simple_norm_all/",
             dataset_replace_string = "transcriptomic_simple_norm_")


classification_models <- unique(model_results$Model)
for (cm in classification_models) {
  individual_model <- model_results %>%
    filter(Model == cm)
  model_barplot <- ggplot(individual_model, aes(x=DataSetId, fill=FSM, y=Mean_AUC)) +
    geom_bar(stat="identity", position="dodge") +
    geom_errorbar( aes(x=DataSetId, ymin=X95.CI_AUC_lower, ymax=X95.CI_AUC_upper), position="dodge") +
    scale_fill_viridis_d() +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
          axis.text.y = element_text(size=rel(1.2), face="italic", hjust=0.95),
          strip.text = element_text(size=rel(1.2), face="bold"),
          legend.title = element_text(size=rel(1.2)),
          legend.text = element_text(size=rel(1.1))
    ) +
    facet_wrap(facets = vars(Model))  
  plotname <- paste(str_replace(cm, " ", ""), "barplot.png", sep = "_")
  dir_name <- "plots/FEMPipeline/barplot_AUC"
  if(!dir.exists(dir_name)){
    dir.create(dir_name, recursive = TRUE)
  }
  plotname <- paste(dir_name, plotname, sep = "/")
  ggsave(plotname, model_barplot, width=20, height=10, dpi=500)
}


plot_features_count_barplot <- function(dir_path = "") {
  tr_fsm_info <- fsm_info %>%
    filter(DataSetId %in% dataset_vec[1:3])
  pr_fsm_info <- fsm_info %>%
    filter(DataSetId %in% dataset_vec[4:12])
  
  if (dir_path != "" & !dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE)
  }

  features_barplot_tr <- ggplot(tr_fsm_info, aes(x=DataSetId, fill=FSM, y=Mean_Number.of.features)) +
    geom_bar(stat="identity", position="dodge") +
    geom_errorbar( aes(x=DataSetId, ymin=X95.CI_Number.of.features_lower, ymax=X95.CI_Number.of.features_upper), position="dodge") +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
          axis.text.y = element_text(size=rel(1.2), face="italic", hjust=0.95),
          axis.title.x = element_text(size=rel(1.5)),
          axis.title.y = element_text(size=rel(1.5), angle=90),
          legend.title = element_text(size=rel(1.5)),
          legend.text = element_text(size=rel(1.2))) +
    scale_y_log10() +
    scale_fill_viridis_d() +
    ylab("Mean number of features used")
  ggsave(append_path(dir_path, "features_count_barplot_tr.png"),
         features_barplot_tr, width=12, height=12, dpi=500)
  
  features_barplot_pr <- ggplot(pr_fsm_info, aes(x=DataSetId, fill=FSM, y=Mean_Number.of.features)) +
    geom_bar(stat="identity", position="dodge") +
    geom_errorbar( aes(x=DataSetId, ymin=X95.CI_Number.of.features_lower, ymax=X95.CI_Number.of.features_upper), position="dodge") +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
          axis.text.y = element_text(size=rel(1.2), face="italic", hjust=0.95),
          axis.title.x = element_text(size=rel(1.5)),
          axis.title.y = element_text(size=rel(1.5), angle=90),
          legend.title = element_text(size=rel(1.5)),
          legend.text = element_text(size=rel(1.2))) +
    scale_y_log10() +
    scale_fill_viridis_d() +
    ylab("Mean number of features used")
  ggsave(append_path(dir_path, "features_count_barplot_pr.png"),
         features_barplot_pr, width=12, height=12, dpi=500)
}

# dir_path = "plots/FEMPipeline/feature_count"
plot_features_count_barplot(dir_path = "plots/FEMPipeline/feature_count")


plot_JI_heatmap <- function(filename = "all_ji.csv", dir = "JI", 
                            heatmapfiledir = "",
                            heatmapfilename = "ji.svg",
                            dparg_vec,
                            fsm_allowed,
                            dataset_replace_string = "") {
  
  dataset_allowed <- c()
  for(id in dparg_vec){
    dataset_id <- paste(dataset_pipeline_arguments[[id]]$dataset_id,
                        dataset_pipeline_arguments[[id]]$classification_criteria,
                        sep = "_")
    dataset_allowed <- append(dataset_allowed, dataset_id)
  }
  print(dataset_allowed)
  
  all_ji_df <- read.table(append_path(dir, filename), sep = ',', header = TRUE)
  all_ji_df <- all_ji_df %>%
    filter(FSM1 == FSM2) %>%
    select(-c(FSM2)) %>%
    filter(FSM1 %in% fsm_allowed) %>%
    filter(DataSetId %in% dataset_allowed) %>%
    mutate(FSM1 = factor(FSM1, levels = fsm_allowed)) %>%
    mutate(DataSetId = gsub("GBMPlasmaEV_", "", DataSetId)) %>%
    mutate(DataSetId = gsub(dataset_replace_string, "", DataSetId)) %>%
    pivot_wider(names_from = DataSetId, values_from = JI) %>%
    column_to_rownames(var = "FSM1") %>%
    drop_na()
  
  all_ji_mat <- data.matrix(all_ji_df)
  
  row_ha <- rowAnnotation(methods = anno_boxplot(all_ji_mat))
  col_ha <- HeatmapAnnotation(datasets = anno_boxplot(all_ji_mat))
  
  if(!dir.exists(heatmapfiledir)){
    dir.create(heatmapfiledir, recursive = TRUE)
  }
  
  svg(append_path(heatmapfiledir, heatmapfilename), 
      width = 8, height = 8)
  ht <- Heatmap(all_ji_mat, name = "Jaccard Index",
                col = viridis(10),
                rect_gp = gpar(col = "white", lwd = 1),
                row_names_side = "left", show_row_dend = FALSE, show_column_dend = FALSE,
                column_names_rot = 45,
                row_names_max_width = max_text_width(rownames(all_ji_mat),
                                                     gp = gpar(fontsize = 12)),
                column_title = "Average pairwise Jaccard Index",
                column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                top_annotation = col_ha, right_annotation = row_ha,
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.3f", all_ji_mat[i, j]), x, y, gp = gpar(fontsize = 7))
                })
  draw(ht)
  dev.off()
}

setwd("~/UNSW/VafaeeLab/GBMPlasmaEV/")

plot_JI_heatmap(filename = "all_ji.csv",
                dir = "fem_pipeline_results/JI/",
                heatmapfilename = "prot_JI.svg",
                heatmapfiledir = "plots/FEMPipeline/JI",
                fsm_allowed = fsm_vec,
                dparg_vec = c(41:43, 126:127, 131, 133, 137),
                dataset_replace_string = "proteomic_impute50fil_")

##########################


setwd("~/UNSW/VafaeeLab/GBMPlasmaEV/")
library(tidyverse)
library(viridis)
library(ComplexHeatmap)
source("scripts/R/utils.R")



setwd("~/UNSW/VafaeeLab/GBMPlasmaEV/")


# dparg_vec <- c(31:35)
# criteria <- 28

dparg_vec <- c(159:162)
criteria <- 29
results_dir = "fem_pipeline_results_subset"

#make sure that only subsets from same classification criteria runs in one call,
#   else the dir path will be incorrect
plot_heatmap_and_var_plot <- function(dparg_vec, criteria, 
                                      results_dir = "fem_pipeline_results"){
  data_info <- read.table(paste(results_dir, "data_info.csv", sep = "/"), 
                          sep = ',', header = TRUE)
  fsm_info <- read.table(paste(results_dir, "fsm_info.csv", sep = "/"),
                         sep = ',', header = TRUE)
  model_results <- read.table(paste(results_dir, "model_results.csv", sep = "/"),
                              sep = ',', header = TRUE)
  
  
  dataset_id_vec <- c()
  for(id in dparg_vec){
    dataset_id <- paste(dataset_pipeline_arguments[[id]]$dataset_id,
                        dataset_pipeline_arguments[[id]]$classification_criteria,
                        sep = "_")
    print(dataset_id)
    dataset_id_vec <- append(dataset_id_vec, dataset_id)
  }
  dataset_id_vec 
  
  subset_file_name_components <- strsplit(dataset_id_vec[1],split = "_")[[1]]
  
  model_results <- model_results %>%
    filter(DataSetId %in% dataset_id_vec) %>%
    mutate(DataSetId = gsub("GBMPlasmaEV_", "", DataSetId, fixed = TRUE)) %>%
    mutate(DataSetId = gsub(subset_file_name_components[2], "", DataSetId, fixed = TRUE)) %>%
    mutate(DataSetId = gsub(subset_file_name_components[3], "", DataSetId, fixed = TRUE)) %>%
    mutate(DataSetId = gsub("__", "", DataSetId, fixed = TRUE))
  
  model_results <- model_results %>%
    mutate(DataSetId = gsub("_PREOPEVsMET", "", DataSetId, fixed = TRUE))
  
  data_to_plot <- model_results %>%
    select(DataSetId, Model, Mean_AUC) %>%
    pivot_wider(names_from = DataSetId, values_from = Mean_AUC) %>%
    column_to_rownames("Model")
  data_to_plot <- data.matrix(data_to_plot)
  
  dir_path <- "plots/FEMPipeline/common_heatmap/tr_preope_met/"
  
  dir_path <- paste0("plots/FEMPipeline/common_heatmap/",
                     paste(subset_file_name_components[2],
                           subset_file_name_components[3],
                           sep = "_"),
                     "/")
  if(!dir.exists(dir_path)){
    dir.create(dir_path)
  }
  file_name <- paste0(dir_path, criteria, "_AUC.png")
  
  heatmap_title <- paste("Mean AUC results from selected common features")
  
  png(file_name, units = "cm", width = 20, height = 15, res = 1200)
  ht <- Heatmap(data_to_plot, name = "Mean AUC",
                col = magma(5),
                rect_gp = gpar(col = "white", lwd = 1),
                column_names_rot = 60,
                column_names_gp = gpar(fontsize = 10),
                row_title = "Feature Selection Methods",
                row_names_side = "left",
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.4f", data_to_plot[i, j]), x, y, gp = gpar(fontsize = 10, col = "slateblue3"))
                })
  draw(ht, column_title = heatmap_title)
  dev.off()
  
  
  dir_path <- paste0("plots/FEMPipeline/common_barplot/",
                     paste(subset_file_name_components[2],
                           subset_file_name_components[3],
                           sep = "_"),
                     "/")
  if(!dir.exists(dir_path)){
    dir.create(dir_path)
  }
  classification_models <- unique(model_results$Model)
  for (cm in classification_models) {
    individual_model <- model_results %>%
      filter(Model == cm)
    model_barplot <- ggplot(individual_model, aes(x=DataSetId, fill=FSM, y=Mean_AUC)) +
      geom_bar(stat="identity", position="dodge") +
      geom_errorbar( aes(x=DataSetId, ymin=X95.CI_AUC_lower, ymax=X95.CI_AUC_upper), position="dodge") +
      scale_fill_viridis_d() +
      theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
            axis.text.y = element_text(size=rel(1.2), face="italic", hjust=0.95),
            strip.text = element_text(size=rel(1.2), face="bold"),
            legend.title = element_text(size=rel(1.2)),
            legend.text = element_text(size=rel(1.1))
      ) +
      facet_wrap(facets = vars(Model))  

    file_name <- paste0(dir_path, criteria, "_", str_replace(cm, " ", ""), "_barplot.png")
    ggsave(file_name, model_barplot, width=20, height=10, dpi=500)
  }  
}

plot_heatmap_and_var_plot(c(31:35), 28)

plot_heatmap_and_var_plot(c(73:74), 29, "fem_pipeline_results_subset")

plot_heatmap_and_var_plot(c(31:35), 28, "fem_pipeline_results_Nov18/fem_pipeline_results")
plot_heatmap_and_var_plot(c(36:40), 29, "fem_pipeline_results_Nov18/fem_pipeline_results")

plot_heatmap_and_var_plot(c(75:79), 28, "fem_pipeline_results_Nov18/fem_pipeline_results")
plot_heatmap_and_var_plot(c(80:84), 29, "fem_pipeline_results_Nov18/fem_pipeline_results")

plot_heatmap_and_var_plot(c(85:89), 28, "fem_pipeline_results_Nov18/fem_pipeline_results")
plot_heatmap_and_var_plot(c(90:94), 29, "fem_pipeline_results_Nov18/fem_pipeline_results")

plot_heatmap_and_var_plot(c(95:99), 28, "fem_pipeline_results_Nov18/fem_pipeline_results")
plot_heatmap_and_var_plot(c(100:104), 29, "fem_pipeline_results_Nov18/fem_pipeline_results")

plot_heatmap_and_var_plot(c(105:109), 28, "fem_pipeline_results_Nov18/fem_pipeline_results")
plot_heatmap_and_var_plot(c(110:114), 29, "fem_pipeline_results_Nov18/fem_pipeline_results")

plot_heatmap_and_var_plot(c(115:119), 28, "fem_pipeline_results_Nov18/fem_pipeline_results")
plot_heatmap_and_var_plot(c(120:124), 29, "fem_pipeline_results_Nov18/fem_pipeline_results")


fsm <- "ranger"
min_iter_feature_presence <- 28
file_name <- "fem_pipeline_results/GBMPlasmaEV_transcriptomic_PREOPEVsMET_ranger_impu_cor_28_PREOPEVsMET_feature_imp.csv" 


plot_feature_imp <- function(fsm, min_iter_feature_presence, file_name){
  feature_imp_data <- read.table(file_name,
                                 sep = ",", header = TRUE)
  feature_imp_data <- feature_imp_data %>%
    select(-c(FSM))
  
  data_to_plot <- feature_imp_data
  
  ggplot(data_to_plot, aes(x = feature, 
                           y = MeanDecreaseGini)) +
    geom_boxplot() +
    ggtitle(paste0("Common features selected by ", fsm, " across ", min_iter_feature_presence, " iterations")) +
    theme(axis.text.x = element_text(size=rel(1.2), angle = 90, hjust = 1),
          axis.text.y = element_text(size=rel(1.2)),
          axis.title.x = element_text(size=rel(1.5)),
          axis.title.y = element_text(size=rel(1.5))) 
  
  dir_path <- paste0("plots/FEMPipeline/feature_imp/")
  if(!dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE)
  }
  
  plot_file_name <- paste0(dir_path,
                           fsm, "_", min_iter_feature_presence,
                           ".png")
  ggsave(plot_file_name)
  
}


plot_feature_imp("ranger", 28, 
                 "fem_pipeline_results/GBMPlasmaEV_transcriptomic_PREOPEVsMET_ranger_impu_cor_28_PREOPEVsMET_feature_imp.csv")

plot_feature_imp("ranger", 29, 
                 "fem_pipeline_results/GBMPlasmaEV_transcriptomic_PREOPEVsMET_ranger_impu_cor_29_PREOPEVsMET_feature_imp.csv")
