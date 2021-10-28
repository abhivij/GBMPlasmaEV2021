dataset_pipeline_arguments <- list(
  
  #1
  #transcriptomic PREOPEVsMET
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all", 
                       "t-test", "t-test_BH",
                       "t-test_pval_0.025", "t-test_pval_0.01", "t-test_pval_0.005",
                       "wilcoxontest", "wilcoxontest_BH",
                       "wilcoxontest_pval_0.025", "wilcoxontest_pval_0.001", "wilcoxontest_pval_0.005",
                       "ranger_impu_cor", 
                       "mrmr10", "mrmr20",
                       "mrmr30", "mrmr50", 
                       "mrmr75", "mrmr100")),

  #2
  #transcriptomic PREOPEVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all", 
                       "t-test", "t-test_BH",
                       "t-test_pval_0.025", "t-test_pval_0.01", "t-test_pval_0.005",
                       "wilcoxontest", "wilcoxontest_BH",
                       "wilcoxontest_pval_0.025", "wilcoxontest_pval_0.001", "wilcoxontest_pval_0.005",
                       "ranger_impu_cor", 
                       "mrmr10", "mrmr20",
                       "mrmr30", "mrmr50", 
                       "mrmr75", "mrmr100")),
  
  #3
  #transcriptomic METVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all", 
                       "t-test", "t-test_BH",
                       "t-test_pval_0.025", "t-test_pval_0.01", "t-test_pval_0.005",
                       "wilcoxontest", "wilcoxontest_BH",
                       "wilcoxontest_pval_0.025", "wilcoxontest_pval_0.001", "wilcoxontest_pval_0.005",
                       "ranger_impu_cor", 
                       "mrmr10", "mrmr20",
                       "mrmr30", "mrmr50", 
                       "mrmr75", "mrmr100")),  

  
  #4
  #transcriptomic PREOPEVsMET
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf")),
  
  #5
  #transcriptomic PREOPEVsHC
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf")),
  
  #6
  #transcriptomic METVsHC
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf")),
  
  #7
  #proteomic PREOPEVsMET
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_quantile",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all", 
                       "t-test", "t-test_BH",
                       "t-test_pval_0.025", "t-test_pval_0.01", "t-test_pval_0.005",
                       "wilcoxontest", "wilcoxontest_BH",
                       "wilcoxontest_pval_0.025", "wilcoxontest_pval_0.001", "wilcoxontest_pval_0.005",
                       "ranger_impu_cor", 
                       "mrmr10", "mrmr20",
                       "mrmr30", "mrmr50", 
                       "mrmr75", "mrmr100"),
       perform_filter = FALSE,
       norm = "quantile"),
  
  #8
  #proteomic PREOPEVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_quantile",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all", 
                       "t-test", "t-test_BH",
                       "t-test_pval_0.025", "t-test_pval_0.01", "t-test_pval_0.005",
                       "wilcoxontest", "wilcoxontest_BH",
                       "wilcoxontest_pval_0.025", "wilcoxontest_pval_0.001", "wilcoxontest_pval_0.005",
                       "ranger_impu_cor", 
                       "mrmr10", "mrmr20",
                       "mrmr30", "mrmr50", 
                       "mrmr75", "mrmr100"),
       perform_filter = FALSE,
       norm = "quantile"),
  
  #9
  #proteomic METVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_quantile",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all", 
                       "t-test", "t-test_BH",
                       "t-test_pval_0.025", "t-test_pval_0.01", "t-test_pval_0.005",
                       "wilcoxontest", "wilcoxontest_BH",
                       "wilcoxontest_pval_0.025", "wilcoxontest_pval_0.001", "wilcoxontest_pval_0.005",
                       "ranger_impu_cor", 
                       "mrmr10", "mrmr20",
                       "mrmr30", "mrmr50", 
                       "mrmr75", "mrmr100"),
       perform_filter = FALSE,
       norm = "quantile"),  
  
  
  #10
  #proteomic PREOPEVsMET
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_quantile",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "quantile"),
  
  #11
  #proteomic PREOPEVsHC
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_quantile",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "quantile"),
  
  #12
  #proteomic METVsHC
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_quantile",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "quantile"),
  
  #13
  #proteomic PREOPEVsMET
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_norm_quantile",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all", 
                       "t-test", "t-test_BH",
                       "t-test_pval_0.025", "t-test_pval_0.01", "t-test_pval_0.005",
                       "wilcoxontest", "wilcoxontest_BH",
                       "wilcoxontest_pval_0.025", "wilcoxontest_pval_0.001", "wilcoxontest_pval_0.005",
                       "ranger_impu_cor", 
                       "mrmr10", "mrmr20",
                       "mrmr30", "mrmr50", 
                       "mrmr75", "mrmr100"),
       perform_filter = FALSE,
       norm = "norm_quantile"),
  
  #14
  #proteomic PREOPEVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_norm_quantile",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all", 
                       "t-test", "t-test_BH",
                       "t-test_pval_0.025", "t-test_pval_0.01", "t-test_pval_0.005",
                       "wilcoxontest", "wilcoxontest_BH",
                       "wilcoxontest_pval_0.025", "wilcoxontest_pval_0.001", "wilcoxontest_pval_0.005",
                       "ranger_impu_cor", 
                       "mrmr10", "mrmr20",
                       "mrmr30", "mrmr50", 
                       "mrmr75", "mrmr100"),
       perform_filter = FALSE,
       norm = "norm_quantile"),
  
  #15
  #proteomic METVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_norm_quantile",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all", 
                       "t-test", "t-test_BH",
                       "t-test_pval_0.025", "t-test_pval_0.01", "t-test_pval_0.005",
                       "wilcoxontest", "wilcoxontest_BH",
                       "wilcoxontest_pval_0.025", "wilcoxontest_pval_0.001", "wilcoxontest_pval_0.005",
                       "ranger_impu_cor", 
                       "mrmr10", "mrmr20",
                       "mrmr30", "mrmr50", 
                       "mrmr75", "mrmr100"),
       perform_filter = FALSE,
       norm = "norm_quantile"),  
  
  
  #16
  #proteomic PREOPEVsMET
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_norm_quantile",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "norm_quantile"),
  
  #17
  #proteomic PREOPEVsHC
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_norm_quantile",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "norm_quantile"),
  
  #18
  #proteomic METVsHC
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_norm_quantile",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "norm_quantile"),      
  
  #19
  # imputed data after 75% NA filter
  #proteomic PREOPEVsMET
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute75fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute75fil_norm_quantile",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all", 
                       "t-test", "t-test_BH",
                       "t-test_pval_0.025", "t-test_pval_0.01", "t-test_pval_0.005",
                       "wilcoxontest", "wilcoxontest_BH",
                       "wilcoxontest_pval_0.025", "wilcoxontest_pval_0.001", "wilcoxontest_pval_0.005",
                       "ranger_impu_cor", 
                       "mrmr10", "mrmr20",
                       "mrmr30", "mrmr50", 
                       "mrmr75", "mrmr100"),
       perform_filter = FALSE,
       norm = "norm_quantile"),
  
  #20
  #proteomic PREOPEVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute75fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute75fil_norm_quantile",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all", 
                       "t-test", "t-test_BH",
                       "t-test_pval_0.025", "t-test_pval_0.01", "t-test_pval_0.005",
                       "wilcoxontest", "wilcoxontest_BH",
                       "wilcoxontest_pval_0.025", "wilcoxontest_pval_0.001", "wilcoxontest_pval_0.005",
                       "ranger_impu_cor", 
                       "mrmr10", "mrmr20",
                       "mrmr30", "mrmr50", 
                       "mrmr75", "mrmr100"),
       perform_filter = FALSE,
       norm = "norm_quantile"),
  
  #21
  #proteomic METVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute75fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute75fil_norm_quantile",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all", 
                       "t-test", "t-test_BH",
                       "t-test_pval_0.025", "t-test_pval_0.01", "t-test_pval_0.005",
                       "wilcoxontest", "wilcoxontest_BH",
                       "wilcoxontest_pval_0.025", "wilcoxontest_pval_0.001", "wilcoxontest_pval_0.005",
                       "ranger_impu_cor", 
                       "mrmr10", "mrmr20",
                       "mrmr30", "mrmr50", 
                       "mrmr75", "mrmr100"),
       perform_filter = FALSE,
       norm = "norm_quantile"),  
  
  
  #22
  #proteomic PREOPEVsMET
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute75fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute75fil_norm_quantile",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "norm_quantile"),
  
  #23
  #proteomic PREOPEVsHC
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute75fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute75fil_norm_quantile",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "norm_quantile"),
  
  #24
  #proteomic METVsHC
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute75fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute75fil_norm_quantile",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "norm_quantile"),
  
  #25
  #proteomic PREOPEVsMET
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all", 
                       "t-test", "t-test_BH",
                       "t-test_pval_0.025", "t-test_pval_0.01", "t-test_pval_0.005",
                       "wilcoxontest", "wilcoxontest_BH",
                       "wilcoxontest_pval_0.025", "wilcoxontest_pval_0.001", "wilcoxontest_pval_0.005",
                       "ranger_impu_cor", 
                       "mrmr10", "mrmr20",
                       "mrmr30", "mrmr50", 
                       "mrmr75", "mrmr100"),
       perform_filter = FALSE,
       norm = "norm_quantile"),
  
  #26
  #proteomic PREOPEVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all", 
                       "t-test", "t-test_BH",
                       "t-test_pval_0.025", "t-test_pval_0.01", "t-test_pval_0.005",
                       "wilcoxontest", "wilcoxontest_BH",
                       "wilcoxontest_pval_0.025", "wilcoxontest_pval_0.001", "wilcoxontest_pval_0.005",
                       "ranger_impu_cor", 
                       "mrmr10", "mrmr20",
                       "mrmr30", "mrmr50", 
                       "mrmr75", "mrmr100"),
       perform_filter = FALSE,
       norm = "norm_quantile"),
  
  #27
  #proteomic METVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all", 
                       "t-test", "t-test_BH",
                       "t-test_pval_0.025", "t-test_pval_0.01", "t-test_pval_0.005",
                       "wilcoxontest", "wilcoxontest_BH",
                       "wilcoxontest_pval_0.025", "wilcoxontest_pval_0.001", "wilcoxontest_pval_0.005",
                       "ranger_impu_cor", 
                       "mrmr10", "mrmr20",
                       "mrmr30", "mrmr50", 
                       "mrmr75", "mrmr100"),
       perform_filter = FALSE,
       norm = "norm_quantile"),  
  
  
  #28
  #proteomic PREOPEVsMET
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "norm_quantile"),
  
  #29
  #proteomic PREOPEVsHC
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "norm_quantile"),
  
  #30
  #proteomic METVsHC
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "norm_quantile"),
  
  #31
  #transcriptomic PREOPEVsMET
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_mrmr30_28_umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_mrmr30_28",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #32
  #transcriptomic PREOPEVsMET
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_mrmr100_28_umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_mrmr100_28",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #33
  #transcriptomic PREOPEVsMET
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_wilcoxontest_28_umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_wilcoxontest_28",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #34
  #transcriptomic PREOPEVsMET
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_ranger_impu_cor_28_umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_ranger_impu_cor_28",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),
  
  #35
  #transcriptomic PREOPEVsMET
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_mrmr30_29_umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_mrmr30_29",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #36
  #transcriptomic PREOPEVsMET
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_mrmr100_29_umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_mrmr100_29",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #37
  #transcriptomic PREOPEVsMET
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_wilcoxontest_29_umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_wilcoxontest_29",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #38
  #transcriptomic PREOPEVsMET
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_ranger_impu_cor_29_umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_ranger_impu_cor_29",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE)    

)  