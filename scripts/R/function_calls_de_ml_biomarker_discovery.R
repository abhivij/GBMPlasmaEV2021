#source('scripts/R/de_ml_biomarker_discovery.R')

RFE_from_ranked_list(data_file_path <- "Data/RNA_all/newquant_Nov2023_umi_counts_PREOPE_MET_HC_filter90.csv",
                     phenotype_file_path <- "Data/transcriptomic_phenotype_PREOPE_MET_HC_withaddicolumn.txt",
                     results_file_path = "DE_results_2024/RFE_results_tra_PREOPEVsHC.csv",
                     omics_type = "transcriptomics",
                     comparison = "PREOPEVsHC",
                     conditions = c("HC", "PREOPE"),
                     ranked_feature_file_path = "DE_results_2024/tra_result_PREOPEVsHC_agg.csv",
                     meanAUC_decrease_cutoff = 0.01,
                     pos_logFC_min_proportion = 0.25)
