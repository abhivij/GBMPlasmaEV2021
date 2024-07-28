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

sample_type <- "test"
comparison_of_interest <- "POSTOPE_TPVsREC_TP"
classes <- c("REC_TP", "POSTOPE_TP")
result_file_dir <- "Data/prediction_result/integration/stacked/"
sample_type <- "test"

#classes assumed to be given as (neg_class, pos_class)
create_mean_prob_heatmap <- function(result_file_dir,
                                     comparison_of_interest,
                                     classes,
                                     sample_type){
  
  results.tra <- read.csv(paste0(result_file_dir, comparison_of_interest, "_tra.csv")) %>%
    filter(Type == sample_type) %>%
    group_by(sample, TrueLabel) %>%
    summarize(mean_prob = mean(Pred_prob))
  results.prot <- read.csv(paste0(result_file_dir, comparison_of_interest, "_prot.csv")) %>%
    filter(Type == sample_type) %>%
    group_by(sample, TrueLabel) %>%
    summarize(mean_prob = mean(Pred_prob))
  results.both <- read.csv(paste0(result_file_dir, comparison_of_interest, "_both.csv")) %>%
    filter(Type == sample_type) %>%
    group_by(sample, TrueLabel) %>%
    summarize(mean_prob = mean(Pred_prob))

  results <- rbind(results.tra %>% mutate(ensemble_type = "stacked with\ntranscriptomics\nmodels"),
                   results.prot %>% mutate(ensemble_type = "stacked with\nproteomics\nmodels"),
                   results.both %>% mutate(ensemble_type = "stacked with\nmulti-omic\nmodels")) %>%
    ungroup()
  
  data_to_plot <- results %>%
    dplyr::select(c(sample, ensemble_type, mean_prob)) %>%
    pivot_wider(names_from = ensemble_type, values_from = mean_prob) %>%
    mutate("Sample_char_part" = str_extract(sample, "[A-Z]+"),
           "Sample_num_part" = as.numeric(str_extract(sample, "[0-9]+"))) %>%
    arrange(Sample_char_part, Sample_num_part) %>%
    dplyr::select(-c(Sample_char_part, Sample_num_part)) %>%
    column_to_rownames("sample")
  data_to_plot <- data.matrix(data_to_plot)
  
  labels <- results %>%
    dplyr::select(c(sample, TrueLabel)) %>%
    distinct() %>%
    ungroup()
  
  meta_data.row <- data.frame(sample = rownames(data_to_plot)) %>% 
    inner_join(labels) %>%
    mutate(TrueLabel = factor(TrueLabel, levels = rev(classes)))

  row_col <- list()
  row_col[["TrueLabel"]] <- c("#440154", "#fde725")
  names(row_col[["TrueLabel"]]) <- classes
  # column_col <- list("Omics type" = c("proteomics" = "skyblue1",
  #                                     "transcriptomics" = "indianred1"))
  # column_col[["Model"]] <- brewer.pal(n = 7, name = "Paired")
  # names(column_col[["Model"]]) <- model_names
  # 
  ht <- Heatmap(data_to_plot, name = "Mean Prediction probability",
          col = viridis(5),
          rect_gp = gpar(col = "white", lwd = 1),
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          row_title = "Samples",
          row_names_side = "left",
          row_split = meta_data.row$TrueLabel,
          column_title = NULL,
          show_column_names = TRUE,
          row_names_gp = gpar(fontsize = 5),
          column_names_rot = 45,
          column_names_centered = TRUE,
          width = unit(9, "cm")) + 
    HeatmapAnnotation("TrueLabel" = meta_data.row$TrueLabel, 
                                 which = "row", 
                      col = row_col)
  
  png(paste0("Data/prediction_result/integration/stacked/", comparison_of_interest, 
             sample_type,
             ".png"), units = "cm", width = 20, height = 15, res = 1200)  
  draw(ht, column_title = paste(sample_type, sub("Vs", " Vs ", comparison_of_interest)))    
  dev.off() 
}

create_mean_prob_heatmap(result_file_dir = "Data/prediction_result/integration/stacked/",
                         comparison_of_interest = "POSTOPE_TPVsREC_TP", 
                         classes = c("REC_TP", "POSTOPE_TP"),
                         sample_type = "train")
create_mean_prob_heatmap(result_file_dir = "Data/prediction_result/integration/stacked/",
                         comparison_of_interest = "PREOPEVsPOSTOPE_TP", 
                         classes = c("POSTOPE_TP", "PREOPE"),
                         sample_type = "train")
create_mean_prob_heatmap(result_file_dir = "Data/prediction_result/integration/stacked/",
                         comparison_of_interest = "PREOPEVsREC_TP", 
                         classes = c("REC_TP", "PREOPE"),
                         sample_type = "train")

create_mean_prob_heatmap(result_file_dir = "Data/prediction_result/integration/stacked/",
                         comparison_of_interest = "POSTOPE_TPVsREC_TP", 
                         classes = c("REC_TP", "POSTOPE_TP"),
                         sample_type = "test")
create_mean_prob_heatmap(result_file_dir = "Data/prediction_result/integration/stacked/",
                         comparison_of_interest = "PREOPEVsPOSTOPE_TP", 
                         classes = c("POSTOPE_TP", "PREOPE"),
                         sample_type = "test")
create_mean_prob_heatmap(result_file_dir = "Data/prediction_result/integration/stacked/",
                         comparison_of_interest = "PREOPEVsREC_TP", 
                         classes = c("REC_TP", "PREOPE"),
                         sample_type = "test")



#heatmap with results from best biomarkers from DE+ML method for prot and tra
#based on manually created file best_features_meanAUC.csv
data_to_plot <- read.csv("DE_results_2024/features/best_features_meanAUC.csv") %>%
  dplyr::select(c(1:3)) %>%
  mutate(Comparison = sub("Vs", " Vs ", Comparison)) %>%
  pivot_wider(names_from = OmicsType, values_from = MeanAUC) %>%
  column_to_rownames("Comparison")
data_to_plot <- as.matrix(data_to_plot)

num_features_to_plot <- read.csv("DE_results_2024/features/best_features_meanAUC.csv") %>%
  dplyr::select(c(1, 2, 5)) %>%
  mutate(Comparison = sub("Vs", " Vs ", Comparison)) %>%
  pivot_wider(names_from = OmicsType, values_from = NumFeatures) %>%
  column_to_rownames("Comparison")
data_to_plot <- as.matrix(data_to_plot)

png("DE_results_2024/features/MeanAUC_with_feature_count.png", units = "cm", width = 30, height = 25, res = 1200)
ht <- Heatmap(data_to_plot, name = "Mean AUC",
              col = plasma(5),
              rect_gp = gpar(col = "white", lwd = 1),
              column_names_rot = 0,
              column_names_centered = TRUE,
              column_names_gp = gpar(fontsize = 12),
              row_names_gp = gpar(fontsize = 12),
              row_title_gp = gpar(fontsize = 15),
              row_title = "Comparison",
              row_names_side = "left",
              column_title_gp = gpar(fontsize = 15), 
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(paste(sprintf("%.3f", data_to_plot[i, j]),
                                paste0("( ", num_features_to_plot[i, j], " biomarkers )"),
                                sep = "\n"), 
                          x, y, gp = gpar(fontsize = 10, col = "gray30", fontface = "bold"))
              })
draw(ht, column_title = "Mean AUC results for best biomarkers")
dev.off()


pdf("DE_results_2024/features/MeanAUC_with_feature_count.pdf")
ht <- Heatmap(data_to_plot, name = "Mean AUC",
              col = plasma(5),
              rect_gp = gpar(col = "white", lwd = 1),
              column_names_rot = 0,
              column_names_centered = TRUE,
              column_names_gp = gpar(fontsize = 12),
              row_names_gp = gpar(fontsize = 12),
              row_title_gp = gpar(fontsize = 15),
              row_title = "Comparison",
              row_names_side = "left",
              column_title_gp = gpar(fontsize = 15), 
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(paste(sprintf("%.3f", data_to_plot[i, j]),
                                paste0("( ", num_features_to_plot[i, j], " biomarkers )"),
                                sep = "\n"), 
                          x, y, gp = gpar(fontsize = 10, col = "gray30", fontface = "bold"))
              })
draw(ht, column_title = "Mean AUC results for best biomarkers")
dev.off()


png("DE_results_2024/features/MeanAUC.png", units = "cm", width = 30, height = 25, res = 1200)
ht <- Heatmap(data_to_plot, name = "Mean AUC",
              col = plasma(5),
              rect_gp = gpar(col = "white", lwd = 1),
              column_names_rot = 0,
              column_names_centered = TRUE,
              column_names_gp = gpar(fontsize = 12),
              row_names_gp = gpar(fontsize = 12),
              row_title_gp = gpar(fontsize = 15),
              row_title = "Comparison",
              row_names_side = "left",
              column_title_gp = gpar(fontsize = 15), 
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.3f", data_to_plot[i, j]),
                          x, y, gp = gpar(fontsize = 10, col = "gray30", fontface = "bold"))
              })
draw(ht, column_title = "Mean AUC results for best biomarkers")
dev.off()


pdf("DE_results_2024/features/MeanAUC.pdf")
ht <- Heatmap(data_to_plot, name = "Mean AUC",
              col = plasma(5),
              rect_gp = gpar(col = "white", lwd = 1),
              column_names_rot = 0,
              column_names_centered = TRUE,
              column_names_gp = gpar(fontsize = 12),
              row_names_gp = gpar(fontsize = 12),
              row_title_gp = gpar(fontsize = 15),
              row_title = "Comparison",
              row_names_side = "left",
              column_title_gp = gpar(fontsize = 15), 
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.3f", data_to_plot[i, j]),
                          x, y, gp = gpar(fontsize = 10, col = "gray30", fontface = "bold"))
              })
draw(ht, column_title = "Mean AUC results for best biomarkers")
dev.off()

############################################################################
# barplot with error bar - performance of best biomarkers
data_to_plot <- read.csv("DE_results_2024/features/best_features_meanAUC.csv") %>%
  dplyr::select(c(1:3)) %>%
  mutate(Comparison = sub("Vs", " Vs ", Comparison))

ggplot(data_to_plot, aes(x = Comparison, fill = Comparison, y = MeanAUC)) +
  geom_bar(stat="identity", position="dodge") +
  xlab("") +
  ylab("Mean AUC") +
  coord_cartesian(ylim = c(0.9, 1.0)) +
  guides(fill=guide_legend(title = "Comparison")) +
  # scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=rel(1.2), hjust=0.95),
        axis.title.x = element_text(size=rel(1.2)),
        axis.title.y = element_text(size=rel(1.2)),
        strip.text = element_text(size=rel(1.2)),
        legend.title = element_text(size=rel(1.2)),
        legend.text = element_text(size=rel(1.1)),
        panel.background = element_rect(colour = "grey50", fill = "white"),
        strip.background = element_rect(colour = "grey50", fill = "white")
  ) +
  facet_wrap(facets = vars(OmicsType))
ggsave("DE_results_2024/features/MeanAUC_barplot_after0.9.png")
ggsave("DE_results_2024/features/MeanAUC_barplot_after0.9.pdf")


ggplot(data_to_plot, aes(x = Comparison, fill = Comparison, y = MeanAUC)) +
  geom_bar(stat="identity", position="dodge") +
  xlab("") +
  ylab("Mean AUC") +
  guides(fill=guide_legend(title = "Comparison")) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=rel(1.2), hjust=0.95),
        axis.title.x = element_text(size=rel(1.2)),
        axis.title.y = element_text(size=rel(1.2)),
        strip.text = element_text(size=rel(1.2)),
        legend.title = element_text(size=rel(1.2)),
        legend.text = element_text(size=rel(1.1)),
        panel.background = element_rect(colour = "grey50", fill = "white"),
        strip.background = element_rect(colour = "grey50", fill = "white")
  ) +
  facet_wrap(facets = vars(OmicsType))
ggsave("DE_results_2024/features/MeanAUC_barplot.png")
ggsave("DE_results_2024/features/MeanAUC_barplot.pdf")

