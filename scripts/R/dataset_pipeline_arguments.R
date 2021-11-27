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
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_t-test_28_umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_t-test_28",
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
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_ranger_impu_cor_28_umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_ranger_impu_cor_28",
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
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_mrmr30_29_umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_mrmr30_29",
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
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_mrmr100_29_umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_mrmr100_29",
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
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_wilcoxontest_29_umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_wilcoxontest_29",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #39
  #transcriptomic PREOPEVsMET
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_t-test_29_umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_t-test_29",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE), 
    
  #40
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
       classifier_feature_imp = TRUE),

  #41
  #proteomic PREOPEVsMET
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_quantile",
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
  
  #42
  #proteomic PREOPEVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_quantile",
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
  
  #43
  #proteomic METVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_quantile",
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
  
  
  #44
  #proteomic PREOPEVsMET
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_quantile",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "quantile"),
  
  #45
  #proteomic PREOPEVsHC
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_quantile",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "quantile"),
  
  #46
  #proteomic METVsHC
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_quantile",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "quantile"),
  
  #47
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
  
  #48
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
  
  #49
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
  
  
  #50
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
  
  #51
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
  
  #52
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
  
  #53
  #proteomic PREOPEVsMET
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_vsn",
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
       norm = "vsn"),
  
  #54
  #proteomic PREOPEVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_vsn",
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
       norm = "vsn"),
  
  #55
  #proteomic METVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_vsn",
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
       norm = "vsn"),  
  
  
  #56
  #proteomic PREOPEVsMET
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_vsn",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "vsn"),
  
  #57
  #proteomic PREOPEVsHC
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_vsn",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "vsn"),
  
  #58
  #proteomic METVsHC
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_vsn",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "vsn"),

  #59
  #transcriptomic PREOPEVsPOSTOPE_T
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "PREOPEVsPOSTOPE_T",
       classes = c("POSTOPE_T", "PREOPE"),
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
  
  #60
  #transcriptomic PREOPEVsPOSTOPE_P
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "PREOPEVsPOSTOPE_P",
       classes = c("POSTOPE_P", "PREOPE"),
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
  
  #61
  #transcriptomic POSTOPE_TVsPOSTOPE_P
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "POSTOPE_TVsPOSTOPE_P",
       classes = c("POSTOPE_P", "POSTOPE_T"),
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
  
  
  #62
  #transcriptomic PREOPEVsPOSTOPE_T
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "PREOPEVsPOSTOPE_T",
       classes = c("POSTOPE_T", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf")),
  
  #63
  #transcriptomic PREOPEVsPOSTOPE_P
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "PREOPEVsPOSTOPE_P",
       classes = c("POSTOPE_P", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf")),
  
  #64
  #transcriptomic POSTOPE_TVsPOSTOPE_P
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "POSTOPE_TVsPOSTOPE_P",
       classes = c("POSTOPE_P", "POSTOPE_T"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf")),
  
  
  #65
  #transcriptomic POSTOPE_TVsREC_T
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "POSTOPE_TVsREC_T",
       classes = c("REC_T", "POSTOPE_T"),
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
  
  #66
  #transcriptomic POSTOPE_TVsREC_T
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "POSTOPE_TVsREC_T",
       classes = c("REC_T", "POSTOPE_T"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf")),
  
  #67
  #transcriptomic POSTOPE_PVsREC_P
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "POSTOPE_PVsREC_P",
       classes = c("REC_P", "POSTOPE_P"),
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
  
  #68
  #transcriptomic POSTOPE_PVsREC_P
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "POSTOPE_PVsREC_P",
       classes = c("REC_P", "POSTOPE_P"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf")),
  
  #69
  #transcriptomic POSTOPE_TVsPREREC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "POSTOPE_TVsPREREC",
       classes = c("PREREC", "POSTOPE_T"),
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
  
  #70
  #transcriptomic POSTOPE_TVsPREREC
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "POSTOPE_TVsPREREC",
       classes = c("PREREC", "POSTOPE_T"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf")),
  
  #71
  #transcriptomic PREOPEVsREC_TP
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "PREOPEVsREC_TP",
       classes = c("REC_TP", "PREOPE"),
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
  
  #72
  #transcriptomic PREOPEVsREC_TP
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "PREOPEVsREC_TP",
       classes = c("REC_TP", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf")),  
  
  
  #73
  #transcriptomic PREOPEVsMET
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_ranger_impu_cor_29_umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_ranger_impu_cor_29",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results_subset",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  

  #74
  #transcriptomic PREOPEVsMET
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_ranger_impu_cor_29_umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsMET_ranger_impu_cor_29_random2000",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results_subset",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE,
       random_seed = 2000),
  
  #75
  #transcriptomic PREOPEVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsHC_mrmr30_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsHC_mrmr30_28",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #76
  #transcriptomic PREOPEVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsHC_mrmr100_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsHC_mrmr100_28",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #77
  #transcriptomic PREOPEVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsHC_wilcoxontest_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsHC_wilcoxontest_28",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),
  
  #78
  #transcriptomic PREOPEVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsHC_t-test_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsHC_t-test_28",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #79
  #transcriptomic PREOPEVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsHC_ranger_impu_cor_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsHC_ranger_impu_cor_28",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),
  
  #80
  #transcriptomic PREOPEVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsHC_mrmr30_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsHC_mrmr30_29",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #81
  #transcriptomic PREOPEVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsHC_mrmr100_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsHC_mrmr100_29",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #82
  #transcriptomic PREOPEVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsHC_wilcoxontest_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsHC_wilcoxontest_29",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #83
  #transcriptomic PREOPEVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsHC_t-test_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsHC_t-test_29",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE), 
  
  #84
  #transcriptomic PREOPEVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_PREOPEVsHC_ranger_impu_cor_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_PREOPEVsHC_ranger_impu_cor_29",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),
  
  #85
  #transcriptomic METVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_METVsHC_mrmr30_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_METVsHC_mrmr30_28",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #86
  #transcriptomic METVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_METVsHC_mrmr100_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_METVsHC_mrmr100_28",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #87
  #transcriptomic METVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_METVsHC_wilcoxontest_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_METVsHC_wilcoxontest_28",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),
  
  #88
  #transcriptomic METVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_METVsHC_t-test_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_METVsHC_t-test_28",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #89
  #transcriptomic METVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_METVsHC_ranger_impu_cor_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_METVsHC_ranger_impu_cor_28",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),
  
  #90
  #transcriptomic METVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_METVsHC_mrmr30_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_METVsHC_mrmr30_29",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #91
  #transcriptomic METVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_METVsHC_mrmr100_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_METVsHC_mrmr100_29",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #92
  #transcriptomic METVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_METVsHC_wilcoxontest_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_METVsHC_wilcoxontest_29",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #93
  #transcriptomic METVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_METVsHC_t-test_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_METVsHC_t-test_29",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE), 
  
  #94
  #transcriptomic METVsHC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "GBMPlasmaEV_transcriptomic_METVsHC_ranger_impu_cor_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic_METVsHC_ranger_impu_cor_29",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),
  
  #95
  #proteomic PREOPEVsMET
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsMET_mrmr30_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_PREOPEVsMET_mrmr30_28",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #96
  #proteomic PREOPEVsMET
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsMET_mrmr100_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_PREOPEVsMET_mrmr100_28",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #97
  #proteomic PREOPEVsMET
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsMET_wilcoxontest_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_PREOPEVsMET_wilcoxontest_28",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),
  
  #98
  #proteomic PREOPEVsMET
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsMET_t-test_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_PREOPEVsMET_t-test_28",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #99
  #proteomic PREOPEVsMET
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsMET_ranger_impu_cor_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_PREOPEVsMET_ranger_impu_cor_28",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),
  
  #100
  #proteomic PREOPEVsMET
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsMET_mrmr30_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_PREOPEVsMET_mrmr30_29",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #101
  #proteomic PREOPEVsMET
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsMET_mrmr100_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_PREOPEVsMET_mrmr100_29",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #102
  #proteomic PREOPEVsMET
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsMET_wilcoxontest_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_PREOPEVsMET_wilcoxontest_29",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #103
  #proteomic PREOPEVsMET
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsMET_t-test_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_PREOPEVsMET_t-test_29",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE), 
  
  #104
  #proteomic PREOPEVsMET
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsMET_ranger_impu_cor_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_PREOPEVsMET_ranger_impu_cor_29",
       classification_criteria = "PREOPEVsMET",
       classes = c("MET", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),
  
  #105
  #proteomic PREOPEVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsHC_mrmr30_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_PREOPEVsHC_mrmr30_28",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #106
  #proteomic PREOPEVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsHC_mrmr100_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_PREOPEVsHC_mrmr100_28",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #107
  #proteomic PREOPEVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsHC_wilcoxontest_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_PREOPEVsHC_wilcoxontest_28",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),
  
  #108
  #proteomic PREOPEVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsHC_t-test_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_PREOPEVsHC_t-test_28",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #109
  #proteomic PREOPEVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsHC_ranger_impu_cor_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_PREOPEVsHC_ranger_impu_cor_28",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),
  
  #110
  #proteomic PREOPEVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsHC_mrmr30_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_PREOPEVsHC_mrmr30_29",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #111
  #proteomic PREOPEVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsHC_mrmr100_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_PREOPEVsHC_mrmr100_29",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #112
  #proteomic PREOPEVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsHC_wilcoxontest_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_PREOPEVsHC_wilcoxontest_29",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #113
  #proteomic PREOPEVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsHC_t-test_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_PREOPEVsHC_t-test_29",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE), 
  
  #114
  #proteomic PREOPEVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_PREOPEVsHC_ranger_impu_cor_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_PREOPEVsHC_ranger_impu_cor_29",
       classification_criteria = "PREOPEVsHC",
       classes = c("HC", "PREOPE"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),
  
  #115
  #proteomic METVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_METVsHC_mrmr30_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_METVsHC_mrmr30_28",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #116
  #proteomic METVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_METVsHC_mrmr100_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_METVsHC_mrmr100_28",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #117
  #proteomic METVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_METVsHC_wilcoxontest_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_METVsHC_wilcoxontest_28",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),
  
  #118
  #proteomic METVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_METVsHC_t-test_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_METVsHC_t-test_28",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #119
  #proteomic METVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_METVsHC_ranger_impu_cor_28.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_METVsHC_ranger_impu_cor_28",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),
  
  #120
  #proteomic METVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_METVsHC_mrmr30_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_METVsHC_mrmr30_29",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #121
  #proteomic METVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_METVsHC_mrmr100_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_METVsHC_mrmr100_29",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #122
  #proteomic METVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_METVsHC_wilcoxontest_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_METVsHC_wilcoxontest_29",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),  
  
  #123
  #proteomic METVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_METVsHC_t-test_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_METVsHC_t-test_29",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE), 
  
  #124
  #proteomic METVsHC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein",
       read_count_file_name = "GBMPlasmaEV_proteomic_impute50fil_norm_quantile_METVsHC_ranger_impu_cor_29.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_METVsHC_ranger_impu_cor_29",
       classification_criteria = "METVsHC",
       classes = c("HC", "MET"),
       cores = 4,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("all"),
       classifier_feature_imp = TRUE),
  

  #125
  #proteomic PREOPEVsPOSTOPE_T
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_quantile",
       classification_criteria = "PREOPEVsPOSTOPE_T",
       classes = c("POSTOPE_T", "PREOPE"),
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
  
  #126
  #proteomic PREOPEVsPOSTOPE_P
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_quantile",
       classification_criteria = "PREOPEVsPOSTOPE_P",
       classes = c("POSTOPE_P", "PREOPE"),
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
  
  #127
  #proteomic POSTOPE_TVsPOSTOPE_P
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_quantile",
       classification_criteria = "POSTOPE_TVsPOSTOPE_P",
       classes = c("POSTOPE_P", "POSTOPE_T"),
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
  
  
  #128
  #proteomic PREOPEVsPOSTOPE_T
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_quantile",
       classification_criteria = "PREOPEVsPOSTOPE_T",
       classes = c("POSTOPE_T", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "quantile"),
  
  #129
  #proteomic PREOPEVsPOSTOPE_P
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_quantile",
       classification_criteria = "PREOPEVsPOSTOPE_P",
       classes = c("POSTOPE_P", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "quantile"),
  
  #130
  #proteomic POSTOPE_TVsPOSTOPE_P
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_quantile",
       classification_criteria = "POSTOPE_TVsPOSTOPE_P",
       classes = c("POSTOPE_P", "POSTOPE_T"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "quantile"),
  
  #131
  #proteomic POSTOPE_TVsREC_T
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_quantile",
       classification_criteria = "POSTOPE_TVsREC_T",
       classes = c("REC_T", "POSTOPE_T"),
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
  
  #132
  #proteomic POSTOPE_TVsREC_T
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_quantile",
       classification_criteria = "POSTOPE_TVsREC_T",
       classes = c("REC_T", "POSTOPE_T"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "quantile"),
  
  #133
  #proteomic POSTOPE_PVsREC_P
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_quantile",
       classification_criteria = "POSTOPE_PVsREC_P",
       classes = c("REC_P", "POSTOPE_P"),
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
  
  #134
  #proteomic POSTOPE_PVsREC_P
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_quantile",
       classification_criteria = "POSTOPE_PVsREC_P",
       classes = c("REC_P", "POSTOPE_P"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "quantile"),      
  
  #135
  #proteomic POSTOPE_TVsPREREC
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_quantile",
       classification_criteria = "POSTOPE_TVsPREREC",
       classes = c("PREREC", "POSTOPE_T"),
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
  
  #136
  #proteomic POSTOPE_TVsPREREC
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_quantile",
       classification_criteria = "POSTOPE_TVsPREREC",
       classes = c("PREREC", "POSTOPE_T"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "quantile"),
  
  #137
  #proteomic PREOPEVsREC_TP
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q7_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_quantile",
       classification_criteria = "PREOPEVsREC_TP",
       classes = c("REC_TP", "PREOPE"),
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
  
  #138
  #proteomic PREOPEVsREC_TP
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q7_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_proteomic_impute50fil_quantile",
       classification_criteria = "PREOPEVsREC_TP",
       classes = c("REC_TP", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       perform_filter = FALSE,
       norm = "quantile"),

  #139
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
                       "mrmr75", "mrmr100"),
       norm = "norm_log_cpm_simple"),
  
  #140
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
                       "mrmr75", "mrmr100"),
       norm = "norm_log_cpm_simple"),
  
  #141
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
                       "mrmr75", "mrmr100"),
       norm = "norm_log_cpm_simple"),  
  
  
  #142
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
       fems_to_run = c("RF_RFE", "ga_rf"),
       norm = "norm_log_cpm_simple"),
  
  #143
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
       fems_to_run = c("RF_RFE", "ga_rf"),
       norm = "norm_log_cpm_simple"),
  
  #144
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
       fems_to_run = c("RF_RFE", "ga_rf"),
       norm = "norm_log_cpm_simple"),
  
  #145
  #transcriptomic PREOPEVsPOSTOPE_T
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "PREOPEVsPOSTOPE_T",
       classes = c("POSTOPE_T", "PREOPE"),
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
       norm = "norm_log_cpm_simple"),
  
  #146
  #transcriptomic PREOPEVsPOSTOPE_P
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "PREOPEVsPOSTOPE_P",
       classes = c("POSTOPE_P", "PREOPE"),
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
       norm = "norm_log_cpm_simple"),
  
  #147
  #transcriptomic POSTOPE_TVsPOSTOPE_P
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "POSTOPE_TVsPOSTOPE_P",
       classes = c("POSTOPE_P", "POSTOPE_T"),
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
       norm = "norm_log_cpm_simple"),  
  
  
  #148
  #transcriptomic PREOPEVsPOSTOPE_T
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "PREOPEVsPOSTOPE_T",
       classes = c("POSTOPE_T", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       norm = "norm_log_cpm_simple"),
  
  #149
  #transcriptomic PREOPEVsPOSTOPE_P
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "PREOPEVsPOSTOPE_P",
       classes = c("POSTOPE_P", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       norm = "norm_log_cpm_simple"),
  
  #150
  #transcriptomic POSTOPE_TVsPOSTOPE_P
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "POSTOPE_TVsPOSTOPE_P",
       classes = c("POSTOPE_P", "POSTOPE_T"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       norm = "norm_log_cpm_simple"),
  
  
  #151
  #transcriptomic POSTOPE_TVsREC_T
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "POSTOPE_TVsREC_T",
       classes = c("REC_T", "POSTOPE_T"),
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
       norm = "norm_log_cpm_simple"),
  
  #152
  #transcriptomic POSTOPE_TVsREC_T
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "POSTOPE_TVsREC_T",
       classes = c("REC_T", "POSTOPE_T"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       norm = "norm_log_cpm_simple"),
  
  #153
  #transcriptomic POSTOPE_PVsREC_P
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "POSTOPE_PVsREC_P",
       classes = c("REC_P", "POSTOPE_P"),
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
       norm = "norm_log_cpm_simple"),
  
  #154
  #transcriptomic POSTOPE_PVsREC_P
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "POSTOPE_PVsREC_P",
       classes = c("REC_P", "POSTOPE_P"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       norm = "norm_log_cpm_simple"),
  
  #155
  #transcriptomic POSTOPE_TVsPREREC
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "POSTOPE_TVsPREREC",
       classes = c("PREREC", "POSTOPE_T"),
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
       norm = "norm_log_cpm_simple"),
  
  #156
  #transcriptomic POSTOPE_TVsPREREC
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "POSTOPE_TVsPREREC",
       classes = c("PREREC", "POSTOPE_T"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       norm = "norm_log_cpm_simple"),
  
  #157
  #transcriptomic PREOPEVsREC_TP
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "PREOPEVsREC_TP",
       classes = c("REC_TP", "PREOPE"),
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
       norm = "norm_log_cpm_simple"),
  
  #158
  #transcriptomic PREOPEVsREC_TP
  # to run "RF_RFE", "ga_rf"
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "GBMPlasmaEV_transcriptomic",
       classification_criteria = "PREOPEVsREC_TP",
       classes = c("REC_TP", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       fems_to_run = c("RF_RFE", "ga_rf"),
       norm = "norm_log_cpm_simple")     

)  
