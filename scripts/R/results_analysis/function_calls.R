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



plot_biomarker_linegraph(conditions = c("PREOPE", "POSTOPE_TP", "PREREC", "REC_TP"),
                         comparison = "PREOPEVsPOSTOPE_TP",
                         omics_type = "transcriptomic",
                         phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP",
                         best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                         plot_dir_path = "plots/FEMPipeline/biomarker_linegraph/PREOPE_POSTOPE_TP_PREREC_REC_TP")
plot_biomarker_linegraph(conditions = c("PREOPE", "POSTOPE_TP", "PREREC", "REC_TP"),
                         comparison = "POSTOPE_TPVsREC_TP",
                         omics_type = "transcriptomic",
                         phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP",
                         best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                         plot_dir_path = "plots/FEMPipeline/biomarker_linegraph/PREOPE_POSTOPE_TP_PREREC_REC_TP")
plot_biomarker_linegraph(conditions = c("PREOPE", "POSTOPE_P", "POSTOPE_T", "PREREC", "REC_TP"),
                         comparison = "PREOPEVsPOSTOPE_TP",
                         omics_type = "transcriptomic",
                         phenotype_column = "PREOPE_POSTOPE_T_POSTOPE_P_PREREC_REC_TP",
                         best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                         plot_dir_path = "plots/FEMPipeline/biomarker_linegraph/PREOPE_POSTOPE_P_POSTOPE_T_PREREC_REC_TP")


create_protein_biomarker_mapping()

plot_biomarker_linegraph(conditions = c("PREOPE", "POSTOPE_TP", "PREREC", "REC_TP"),
                         comparison = "PREOPEVsPOSTOPE_TP",
                         omics_type = "proteomic",
                         phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP",
                         best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                         plot_dir_path = "plots/FEMPipeline/biomarker_linegraph/PREOPE_POSTOPE_TP_PREREC_REC_TP")
plot_biomarker_linegraph(conditions = c("PREOPE", "POSTOPE_TP", "PREREC", "REC_TP"),
                         comparison = "POSTOPE_TPVsREC_TP",
                         omics_type = "proteomic",
                         phenotype_column = "PREOPE_POSTOPE_TP_PREREC_REC_TP",
                         best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                         plot_dir_path = "plots/FEMPipeline/biomarker_linegraph/PREOPE_POSTOPE_TP_PREREC_REC_TP")
plot_biomarker_linegraph(conditions = c("PREOPE", "POSTOPE_P", "POSTOPE_T", "PREREC", "REC_TP"),
                         comparison = "PREOPEVsPOSTOPE_TP",
                         omics_type = "proteomic",
                         phenotype_column = "PREOPE_POSTOPE_T_POSTOPE_P_PREREC_REC_TP",
                         best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                         plot_dir_path = "plots/FEMPipeline/biomarker_linegraph/PREOPE_POSTOPE_P_POSTOPE_T_PREREC_REC_TP")

