setwd("~/UNSW/VafaeeLab/GBMPlasmaEV")
library(tidyverse)
library(viridis)
library(ComplexHeatmap)
source("scripts/R/utils.R")
source("scripts/R/dataset_pipeline_arguments.R")
source("scripts/R/dataset_pipeline_arguments_transcriptomic.R")
source("scripts/R/dataset_pipeline_arguments_proteomic.R")



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
                         dataset_pipeline_arguments = dataset_pipeline_arguments,
                         results_dir = "fem_pipeline_results",
                         dir_path = "",
                         dataset_replace_string = ""){
  data_info <- read.table(paste(results_dir, "data_info.csv", sep = "/"), 
                          sep = ',', header = TRUE)
  fsm_info <- read.table(paste(results_dir, "fsm_info.csv", sep = "/"),
                         sep = ',', header = TRUE)
  model_results <- read.table(paste(results_dir, "model_results_test.csv", sep = "/"),
                              sep = ',', header = TRUE)
  
  if(!dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE)
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
  
  ds_heatmap_dir_path <- paste0(dir_path, "/ds/") 
  if(!dir.exists(ds_heatmap_dir_path)){
    dir.create(ds_heatmap_dir_path)
  }
  for(ds in replaced_dataset_id_vec){
    print(ds)
    
    model_results_sub <- model_results %>%
      filter(Model != "Simple logistic regression") %>%
      mutate(Model = gsub("Regularized logistic regression", "regu log_reg", Model)) %>%
      filter(DataSetId == ds)
    
    data_to_plot <- model_results_sub %>%
      select(Model, FSM, Mean_AUC) %>%
      # pivot_wider(names_from = DataSetId, values_from = Mean_AUC, values_fn = length) %>%          
      pivot_wider(names_from = Model, values_from = Mean_AUC) %>%
      column_to_rownames("FSM")
    data_to_plot <- data.matrix(data_to_plot)
    
    file_name_curr <- paste0("Mean_AUC",
                             gsub(" ", "_", ds),
                             ".png")
    
    file_name_curr <- paste0(ds_heatmap_dir_path, file_name_curr)
    
    png(file_name_curr, units = "cm", width = 20, height = 15, res = 1200)
    ht <- Heatmap(data_to_plot, name = "Mean AUC",
                  col = magma(5),
                  rect_gp = gpar(col = "white", lwd = 1),
                  column_names_rot = 60,
                  column_names_gp = gpar(fontsize = 10),
                  row_title = "Feature Selection Methods",
                  row_names_side = "left",
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(sprintf("%.3f", data_to_plot[i, j]), x, y, gp = gpar(fontsize = 7, col = "slateblue3"))
                  })
    draw(ht, column_title = paste("Mean AUC results for", ds))
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


plot_heatmap(dparg_vec = c(139:141, 145:147, 151, 153, 155, 157,
                           223, 225, 227, 229, 231),   #transcriptomic new set of comparisons
             results_dir = "fem_pipeline_results",
             dir_path = "plots/FEMPipeline/heatmap/transcriptomic_simple_norm_all_2022Feb23/",
             dataset_replace_string = "transcriptomic_simple_norm_")

plot_heatmap(dparg_vec = c(41:43, 125:127, 131, 133, 135, 137,
                           233, 235, 237, 239, 241),  #proteomic new set of comparisons)
             results_dir = "fem_pipeline_results",
             dir_path = "plots/FEMPipeline/heatmap/proteomic_quantile_all_2022Feb23/",
             dataset_replace_string = "proteomic_impute50fil_")


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


plot_features_count_barplot <- function(dparg_vec, dataset_pipeline_arguments,
                                        results_dir_path, output_dir_path) {
  
  data_info <- read.table(paste(results_dir_path, "data_info.csv", sep = "/"), 
                          sep = ',', header = TRUE)
  fsm_info <- read.table(paste(results_dir_path, "fsm_info.csv", sep = "/"),
                         sep = ',', header = TRUE)
  model_results <- read.table(paste(results_dir_path, "model_results_test.csv", sep = "/"),
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
  
  tr_fsm_info <- fsm_info %>%
    filter(DataSetId %in% dataset_id_vec)
  # pr_fsm_info <- fsm_info %>%
  #   filter(DataSetId %in% dataset_vec[4:12])
  
  if (output_dir_path != "" & !dir.exists(output_dir_path)){
    dir.create(output_dir_path, recursive = TRUE)
  }
  
  features_barplot_tr <- ggplot(tr_fsm_info, aes(x=DataSetId, fill=FSM, y=Mean_Number.of.features)) +
    geom_bar(stat="identity", position="dodge") +
    geom_errorbar( aes(x=DataSetId, ymin=X95.CI_Number.of.features_lower, ymax=X95.CI_Number.of.features_upper), position="dodge") +
    theme(axis.text.x = element_text(size=rel(1.2), angle=45, hjust=1, vjust=1),
          axis.text.y = element_text(size=rel(1.2), face="italic", hjust=0.95),
          axis.title.x = element_text(size=rel(1.5)),
          axis.title.y = element_text(size=rel(1.5), angle=90),
          legend.title = element_text(size=rel(1.5)),
          legend.text = element_text(size=rel(1.2))) +
    scale_fill_viridis_d() +
    scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, by = 10)) +
    ylab("Mean number of features used")
  ggsave(append_path(output_dir_path, "features_count_barplot_tr.png"),
         features_barplot_tr, width=12, height=12, dpi=500)
  
  # features_barplot_pr <- ggplot(pr_fsm_info, aes(x=DataSetId, fill=FSM, y=Mean_Number.of.features)) +
  #   geom_bar(stat="identity", position="dodge") +
  #   geom_errorbar( aes(x=DataSetId, ymin=X95.CI_Number.of.features_lower, ymax=X95.CI_Number.of.features_upper), position="dodge") +
  #   theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
  #         axis.text.y = element_text(size=rel(1.2), face="italic", hjust=0.95),
  #         axis.title.x = element_text(size=rel(1.5)),
  #         axis.title.y = element_text(size=rel(1.5), angle=90),
  #         legend.title = element_text(size=rel(1.5)),
  #         legend.text = element_text(size=rel(1.2))) +
  #   scale_y_log10() +
  #   scale_fill_viridis_d() +
  #   ylab("Mean number of features used")
  # ggsave(append_path(output_dir_path, "features_count_barplot_pr.png"),
  #        features_barplot_pr, width=12, height=12, dpi=500)
}

# output_dir_path = "plots/FEMPipeline/feature_count"
plot_features_count_barplot(output_dir_path = "plots/FEMPipeline/feature_count")


plot_JI_heatmap <- function(filename = "all_ji.csv", dir = "JI",
                            dataset_pipeline_arguments = dataset_pipeline_arguments,
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
    column_to_rownames(var = "FSM1")
  
  all_ji_df[is.na(all_ji_df)] <- 0
  
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
                dparg_vec = c(41:43, 125:127, 131, 133, 135, 137,
                              233, 235, 237, 239, 241),
                dataset_replace_string = "proteomic_impute50fil_")

plot_JI_heatmap(filename = "all_ji.csv",
                dir = "fem_pipeline_results/JI/",
                heatmapfilename = "tra_JI.svg",
                heatmapfiledir = "plots/FEMPipeline/JI",
                fsm_allowed = fsm_vec,
                dparg_vec = c(139:141, 145:147, 151, 153, 155, 157,
                              223, 225, 227, 229, 231),
                dataset_replace_string = "transcriptomic_simple_norm_")

##########################


setwd("~/UNSW/VafaeeLab/GBMPlasmaEV/")
library(tidyverse)
library(viridis)
library(ComplexHeatmap)
source("scripts/R/utils.R")



# dparg_vec <- c(31:35)
# criteria <- 28

# dparg_vec <- c(163:165)
# criteria <- 29
# results_dir = "fem_pipeline_results_subset"
# dataset_replace_string = "transcriptomic_simple_norm_PREOPEVsMET_"
# heatmap_file_name = "tr_PREOPEVsMET.png"


plot_common_feature_heatmap <- function(dparg_vec, 
                                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                                      results_dir = "fem_pipeline_results",
                                      dataset_replace_string = "",
                                      heatmap_file_name,
                                      dir_path = "plots/FEMPipeline_new_quant/common_heatmap/"){
  data_info <- read.table(paste(results_dir, "data_info.csv", sep = "/"), 
                          sep = ',', header = TRUE)
  fsm_info <- read.table(paste(results_dir, "fsm_info.csv", sep = "/"),
                         sep = ',', header = TRUE)
  model_results <- read.table(paste(results_dir, "model_results_test.csv", sep = "/"),
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
  
  model_results <- model_results %>%
    filter(DataSetId %in% dataset_id_vec) %>%
    mutate(DataSetId = gsub(dataset_replace_string, "", DataSetId, fixed = TRUE)) %>%
    mutate(DataSetId = gsub("GBMPlasmaEV_", "", DataSetId, fixed = TRUE))
  
  data_to_plot <- model_results %>%
    select(DataSetId, Model, Mean_AUC) %>%
    pivot_wider(names_from = DataSetId, values_from = Mean_AUC) %>%
    column_to_rownames("Model")
  data_to_plot <- data.matrix(data_to_plot)
  
  if(!dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE)
  }
  file_path <- paste0(dir_path, heatmap_file_name)
  
  heatmap_title <- paste("Mean AUC results from selected common features")
  
  png(file_path, units = "cm", width = 20, height = 15, res = 1200)
  ht <- Heatmap(data_to_plot, name = "Mean AUC",
                col = magma(5),
                rect_gp = gpar(col = "white", lwd = 1),
                column_names_rot = 30,
                column_names_gp = gpar(fontsize = 10),
                row_title = "Classification Models",
                row_names_side = "left",
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.4f", data_to_plot[i, j]), x, y, gp = gpar(fontsize = 10, col = "slateblue3"))
                })
  draw(ht, column_title = heatmap_title)
  dev.off()
}


plot_common_feature_heatmap(c(163:165),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "transcriptomic_simple_norm_PREOPEVsMET_",
                            heatmap_file_name = "tr_PREOPEVsMET.png"
)
plot_common_feature_heatmap(c(166:169),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "transcriptomic_simple_norm_PREOPEVsHC_",
                            heatmap_file_name = "tr_PREOPEVsHC.png"
)
plot_common_feature_heatmap(c(170:173),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "transcriptomic_simple_norm_METVsHC_",
                            heatmap_file_name = "tr_METVsHC.png"
)
plot_common_feature_heatmap(c(174:175),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "transcriptomic_simple_norm_PREOPEVsPOSTOPE_T_",
                            heatmap_file_name = "tr_PREOPEVsPOSTOPE_T.png"
)
plot_common_feature_heatmap(c(176:177),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "transcriptomic_simple_norm_PREOPEVsPOSTOPE_P_",
                            heatmap_file_name = "tr_PREOPEVsPOSTOPE_P.png"
)
plot_common_feature_heatmap(c(178:179),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "transcriptomic_simple_norm_POSTOPE_TVsPOSTOPE_P_",
                            heatmap_file_name = "tr_POSTOPE_TVsPOSTOPE_P.png"
)
plot_common_feature_heatmap(c(180:182),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "transcriptomic_simple_norm_POSTOPE_TVsREC_T_",
                            heatmap_file_name = "tr_POSTOPE_TVsREC_T.png"
)
plot_common_feature_heatmap(c(183:184),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "transcriptomic_simple_norm_POSTOPE_PVsREC_P_",
                            heatmap_file_name = "tr_POSTOPE_PVsREC_P.png"
)
plot_common_feature_heatmap(c(185:186),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "transcriptomic_simple_norm_POSTOPE_TVsPREREC_",
                            heatmap_file_name = "tr_POSTOPE_TVsPREREC.png"
)
plot_common_feature_heatmap(c(187:189),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "transcriptomic_simple_norm_PREOPEVsREC_TP_",
                            heatmap_file_name = "tr_PREOPEVsREC_TP.png"
)


plot_common_feature_heatmap(c(190:192),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "proteomic_impute50fil_quantile_PREOPEVsMET_",
                            heatmap_file_name = "pr_PREOPEVsMET.png"
)
plot_common_feature_heatmap(c(193:198),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "proteomic_impute50fil_quantile_PREOPEVsHC_",
                            heatmap_file_name = "pr_PREOPEVsHC.png"
)
plot_common_feature_heatmap(c(199:201),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "proteomic_impute50fil_quantile_METVsHC_",
                            heatmap_file_name = "pr_METVsHC.png"
)

plot_common_feature_heatmap(c(202:205),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "proteomic_impute50fil_quantile_PREOPEVsPOSTOPE_P_",
                            heatmap_file_name = "pr_PREOPEVsPOSTOPE_P.png"
)
plot_common_feature_heatmap(c(206:207),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "proteomic_impute50fil_quantile_POSTOPE_TVsPOSTOPE_P_",
                            heatmap_file_name = "pr_POSTOPE_TVsPOSTOPE_P.png"
)
plot_common_feature_heatmap(c(208:210),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "proteomic_impute50fil_quantile_POSTOPE_TVsREC_T_",
                            heatmap_file_name = "pr_POSTOPE_TVsREC_T.png"
)
plot_common_feature_heatmap(c(211:214),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "proteomic_impute50fil_quantile_POSTOPE_PVsREC_P_",
                            heatmap_file_name = "pr_POSTOPE_PVsREC_P.png"
)

plot_common_feature_heatmap(c(215:218),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "proteomic_impute50fil_quantile_PREOPEVsREC_TP_",
                            heatmap_file_name = "pr_PREOPEVsREC_TP.png"
)

plot_common_feature_heatmap(c(219:220),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "proteomic_impute50fil_quantile_PREOPEVsPOSTOPE_T_",
                            heatmap_file_name = "pr_PREOPEVsPOSTOPE_T.png"
)

plot_common_feature_heatmap(c(221:222),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "proteomic_impute50fil_quantile_POSTOPE_TVsPREREC_",
                            heatmap_file_name = "pr_POSTOPE_TVsPREREC.png"
)



plot_common_feature_heatmap(c(243:245),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "transcriptomic_simple_norm_PREOPEVsPOSTOPE_TP_",
                            heatmap_file_name = "tr_PREOPEVsPOSTOPE_TP.png"
)
plot_common_feature_heatmap(c(246:247),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "transcriptomic_simple_norm_POSTOPE_TPVsREC_TP_",
                            heatmap_file_name = "tr_POSTOPE_TPVsREC_TP.png"
)


plot_common_feature_heatmap(c(248:251),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "proteomic_impute50fil_quantile_PREOPEVsPOSTOPE_TP_",
                            heatmap_file_name = "pr_PREOPEVsPOSTOPE_TP.png"
)

plot_common_feature_heatmap(c(252:254),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "proteomic_impute50fil_quantile_POSTOPE_TPVsREC_TP_",
                            heatmap_file_name = "pr_POSTOPE_TPVsREC_TP.png"
)
# plot_heatmap_and_var_plot(c(31:35), 28)
# 
# plot_heatmap_and_var_plot(c(73:74), 29, "fem_pipeline_results_subset")
# 
# plot_heatmap_and_var_plot(c(31:35), 28, "fem_pipeline_results_Nov18/fem_pipeline_results")
# plot_heatmap_and_var_plot(c(36:40), 29, "fem_pipeline_results_Nov18/fem_pipeline_results")
# 
# plot_heatmap_and_var_plot(c(75:79), 28, "fem_pipeline_results_Nov18/fem_pipeline_results")
# plot_heatmap_and_var_plot(c(80:84), 29, "fem_pipeline_results_Nov18/fem_pipeline_results")
# 
# plot_heatmap_and_var_plot(c(85:89), 28, "fem_pipeline_results_Nov18/fem_pipeline_results")
# plot_heatmap_and_var_plot(c(90:94), 29, "fem_pipeline_results_Nov18/fem_pipeline_results")
# 
# plot_heatmap_and_var_plot(c(95:99), 28, "fem_pipeline_results_Nov18/fem_pipeline_results")
# plot_heatmap_and_var_plot(c(100:104), 29, "fem_pipeline_results_Nov18/fem_pipeline_results")
# 
# plot_heatmap_and_var_plot(c(105:109), 28, "fem_pipeline_results_Nov18/fem_pipeline_results")
# plot_heatmap_and_var_plot(c(110:114), 29, "fem_pipeline_results_Nov18/fem_pipeline_results")
# 
# plot_heatmap_and_var_plot(c(115:119), 28, "fem_pipeline_results_Nov18/fem_pipeline_results")
# plot_heatmap_and_var_plot(c(120:124), 29, "fem_pipeline_results_Nov18/fem_pipeline_results")


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




####################### initial cohort plots on newly quantified results
#### plot_heatmap(dparg_vec = c(41:43, 47:49, 53:55),
plot_heatmap(
  dparg_vec = c(1, 5, 9),
  dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
  results_dir = "fem_pipeline_results_tr",
  dir_path = "plots/FEMPipeline_new_quant/heatmap/transcriptomic/",
  dataset_replace_string = "_initial"
)

plot_JI_heatmap(filename = "all_ji.csv",
                dir = "fem_pipeline_results_tr/JI/",
                dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                heatmapfilename = "tr_JI.svg",
                heatmapfiledir = "plots/FEMPipeline_new_quant/",
                fsm_allowed = fsm_vector,
                dparg_vec = c(1, 5, 9),
                dataset_replace_string = "_initial")

plot_features_count_barplot(dparg_vec = c(1, 5, 9), 
                            dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                            results_dir_path = "fem_pipeline_results_tr", 
                            output_dir_path = "plots/FEMPipeline_new_quant/features_count/transcriptomic/")


plot_common_feature_heatmap(c(13, 14),
                            dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                            results_dir = "fem_pipeline_results_tr_subset",
                            dataset_replace_string = "GBM_tr_initial_",
                            heatmap_file_name = "tr_PREOPEVsPOSTOPE_TP.png")

plot_common_feature_heatmap(c(15, 16),
                            dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                            results_dir = "fem_pipeline_results_tr_subset",
                            dataset_replace_string = "GBM_tr_initial_",
                            heatmap_file_name = "tr_POSTOPE_TPVsREC_TP.png")

plot_common_feature_heatmap(c(17, 18),
                            dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                            results_dir = "fem_pipeline_results_tr_subset",
                            dataset_replace_string = "GBM_tr_initial_",
                            heatmap_file_name = "tr_PREOPEVsREC_TP.png")


plot_common_feature_heatmap(dparg_vec = c(17, 18),
                            dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                            results_dir = "fem_pipeline_results_tr_subset",
                            dataset_replace_string = "GBM_tr_initial_",
                            heatmap_file_name = "tr_PREOPEVsREC_TP.png", 
                            dir_path = )


  
  

##########

#proteomic - with no norm

plot_heatmap(
  dparg_vec = c(1, 5, 9),
  dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
  results_dir = "fem_pipeline_results_pr",
  dir_path = "plots/FEMPipeline_prot_no_norm/",
  dataset_replace_string = "GBM_initial_proteomic_impute50fil_"
)




##########

#proteomic - with quantile norm with train param

plot_heatmap(
  dparg_vec = c(13, 17, 21),
  dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
  results_dir = "fem_pipeline_results_pr",
  dir_path = "plots/FEMPipeline_prot_quantile_norm_with_train_param/",
  dataset_replace_string = "GBM_initial_proteomic_impute50fil_"
)




#### proteomic subsets

plot_common_feature_heatmap(c(25, 26),
                            dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                            results_dir = "fem_pipeline_results_proteomic_subset",
                            dataset_replace_string = "GBM_initial_proteomic_impute50fil_no_norm_PREOPEVsPOSTOPE_TP_",
                            dir_path = "plots/FEMPipeline_prot_no_norm/common_heatmap/",
                            heatmap_file_name = "PREOPEVsPOSTOPE_TP.png"
)

plot_common_feature_heatmap(c(27, 28),
                            dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                            results_dir = "fem_pipeline_results_proteomic_subset",
                            dataset_replace_string = "GBM_initial_proteomic_impute50fil_no_norm_POSTOPE_TPVsREC_TP_",
                            dir_path = "plots/FEMPipeline_prot_no_norm/common_heatmap/",
                            heatmap_file_name = "POSTOPE_TPVsREC_TP.png"
)

plot_common_feature_heatmap(c(29, 30),
                            dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                            results_dir = "fem_pipeline_results_proteomic_subset",
                            dataset_replace_string = "GBM_initial_proteomic_impute50fil_no_norm_PREOPEVsREC_TP_",
                            dir_path = "plots/FEMPipeline_prot_no_norm/common_heatmap/",
                            heatmap_file_name = "PREOPEVsREC_TP.png"
)

plot_common_feature_heatmap(c(31, 32),
                            dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                            results_dir = "fem_pipeline_results_proteomic_subset",
                            dataset_replace_string = "GBM_initial_proteomic_impute50fil_quantile_train_param_PREOPEVsPOSTOPE_TP_",
                            dir_path = "plots/FEMPipeline_prot_quantile_norm_with_train_param/common_heatmap/",
                            heatmap_file_name = "PREOPEVsPOSTOPE_TP.png"
)

plot_common_feature_heatmap(c(33, 34),
                            dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                            results_dir = "fem_pipeline_results_proteomic_subset",
                            dataset_replace_string = "GBM_initial_proteomic_impute50fil_quantile_train_param_POSTOPE_TPVsREC_TP_",
                            dir_path = "plots/FEMPipeline_prot_quantile_norm_with_train_param/common_heatmap/",
                            heatmap_file_name = "POSTOPE_TPVsREC_TP.png"
)

plot_common_feature_heatmap(c(35, 36),
                            dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                            results_dir = "fem_pipeline_results_proteomic_subset",
                            dataset_replace_string = "GBM_initial_proteomic_impute50fil_quantile_train_param_PREOPEVsREC_TP_",
                            dir_path = "plots/FEMPipeline_prot_quantile_norm_with_train_param/common_heatmap/",
                            heatmap_file_name = "PREOPEVsREC_TP.png"
)


####################


##########

#proteomic common - with no norm

plot_heatmap(
  dparg_vec = c(37, 41, 45),
  dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
  results_dir = "fem_pipeline_results_pr_common",
  dir_path = "plots/FEMPipeline_prot_common_no_norm/",
  dataset_replace_string = "GBM_initial_proteomic_impute50fil_"
)

plot_common_feature_heatmap(c(61),
                            dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                            results_dir = "fem_pipeline_results_proteomic_common_none_subset",
                            dataset_replace_string = "GBM_initial_proteomic_impute50fil_common_no_norm_POSTOPE_TPVsREC_TP_",
                            dir_path = "plots/FEMPipeline_prot_common_no_norm/common_heatmap/",
                            heatmap_file_name = "POSTOPE_TPVsREC_TP.png"
)
plot_common_feature_heatmap(c(62, 63),
                            dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                            results_dir = "fem_pipeline_results_proteomic_common_none_subset",
                            dataset_replace_string = "GBM_initial_proteomic_impute50fil_common_no_norm_PREOPEVsPOSTOPE_TP_",
                            dir_path = "plots/FEMPipeline_prot_common_no_norm/common_heatmap/",
                            heatmap_file_name = "PREOPEVsPOSTOPE_TP.png"
)
plot_common_feature_heatmap(c(64),
                            dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                            results_dir = "fem_pipeline_results_proteomic_common_none_subset",
                            dataset_replace_string = "GBM_initial_proteomic_impute50fil_common_no_norm_PREOPEVsREC_TP_",
                            dir_path = "plots/FEMPipeline_prot_common_no_norm/common_heatmap/",
                            heatmap_file_name = "PREOPEVsREC_TP.png"
)


##########

#proteomic common - with quantile norm with train param

plot_heatmap(
  dparg_vec = c(49, 53, 57),
  dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
  results_dir = "fem_pipeline_results_pr_common",
  dir_path = "plots/FEMPipeline_prot_common_quantile_norm_with_train_param/",
  dataset_replace_string = "GBM_initial_proteomic_impute50fil_"
)

plot_common_feature_heatmap(c(65, 66),
                            dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                            results_dir = "fem_pipeline_results_proteomic_common_quantile_train_param_subset",
                            dataset_replace_string = "GBM_initial_proteomic_impute50fil_common_quantile_train_param_POSTOPE_TPVsREC_TP_",
                            dir_path = "plots/FEMPipeline_prot_common_quantile_norm_with_train_param/common_heatmap/",
                            heatmap_file_name = "POSTOPE_TPVsREC_TP.png"
)
plot_common_feature_heatmap(c(67),
                            dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                            results_dir = "fem_pipeline_results_proteomic_common_quantile_train_param_subset",
                            dataset_replace_string = "GBM_initial_proteomic_impute50fil_common_quantile_train_param_PREOPEVsPOSTOPE_TP_",
                            dir_path = "plots/FEMPipeline_prot_common_quantile_norm_with_train_param/common_heatmap/",
                            heatmap_file_name = "PREOPEVsPOSTOPE_TP.png"
)
plot_common_feature_heatmap(c(68),
                            dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                            results_dir = "fem_pipeline_results_proteomic_common_quantile_train_param_subset",
                            dataset_replace_string = "GBM_initial_proteomic_impute50fil_common_quantile_train_param_PREOPEVsREC_TP_",
                            dir_path = "plots/FEMPipeline_prot_common_quantile_norm_with_train_param/common_heatmap/",
                            heatmap_file_name = "PREOPEVsREC_TP.png"
)
