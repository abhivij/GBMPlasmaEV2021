library(tidyverse)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

source("scripts/R/results_analysis/plot_biomarker_linegraph.R")

# comparison <- "METVsHC"
# omics_type = "transcriptomic"
# phenotype_column = "HC_MET_PREOPE_REC_TP"
# conditions <- c("HC", "MET", "PREOPE", "REC_TP")
# best_features_file_path <- "Data/selected_features/best_features_with_add_col.csv"
# plot_dir_path <- "plots/FEMPipeline/biomarker_linegraph/HC_MET_PREOPE_REC_TP"

plot_biomarker_linegraph(conditions = c("HC", "MET", "PREOPE", "REC_TP"),
                         comparison = "METVsHC",
                         omics_type = "transcriptomic",
                         phenotype_column = "HC_MET_PREOPE_REC_TP",
                         best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                         plot_dir_path = "plots/FEMPipeline/biomarker_linegraph/HC_MET_PREOPE_REC_TP")
plot_biomarker_linegraph(conditions = c("HC", "MET", "PREOPE", "REC_TP"),
                         comparison = "METVsHC",
                         omics_type = "proteomic",
                         phenotype_column = "HC_MET_PREOPE_REC_TP",
                         best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                         plot_dir_path = "plots/FEMPipeline/biomarker_linegraph/HC_MET_PREOPE_REC_TP")
