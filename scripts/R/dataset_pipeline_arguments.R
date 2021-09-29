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
       fems_to_run = c("all", "t-test", "t-test_BH",
                       "wilcoxontest", "wilcoxontest_BH",
                       "ranger_impu_cor", 
                       "mrmr30", "mrmr50", "mrmr75", "mrmr100")),

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
       fems_to_run = c("all", "t-test", "t-test_BH",
                       "wilcoxontest", "wilcoxontest_BH",
                       "ranger_impu_cor", 
                       "mrmr30", "mrmr50", "mrmr75", "mrmr100")),
  
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
       fems_to_run = c("all", "t-test", "t-test_BH",
                       "wilcoxontest", "wilcoxontest_BH",
                       "ranger_impu_cor", 
                       "mrmr30", "mrmr50", "mrmr75", "mrmr100")),  

  
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
       fems_to_run = c("all", "t-test", "t-test_BH",
                       "wilcoxontest", "wilcoxontest_BH",
                       "ranger_impu_cor", 
                       "mrmr30", "mrmr50", "mrmr75", "mrmr100"),
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
       fems_to_run = c("all", "t-test", "t-test_BH",
                       "wilcoxontest", "wilcoxontest_BH",
                       "ranger_impu_cor", 
                       "mrmr30", "mrmr50", "mrmr75", "mrmr100"),
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
       fems_to_run = c("all", "t-test", "t-test_BH",
                       "wilcoxontest", "wilcoxontest_BH",
                       "ranger_impu_cor", 
                       "mrmr30", "mrmr50", "mrmr75", "mrmr100"),
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
       fems_to_run = c("all", "t-test", "t-test_BH",
                       "wilcoxontest", "wilcoxontest_BH",
                       "ranger_impu_cor", 
                       "mrmr30", "mrmr50", "mrmr75", "mrmr100"),
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
       fems_to_run = c("all", "t-test", "t-test_BH",
                       "wilcoxontest", "wilcoxontest_BH",
                       "ranger_impu_cor", 
                       "mrmr30", "mrmr50", "mrmr75", "mrmr100"),
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
       fems_to_run = c("all", "t-test", "t-test_BH",
                       "wilcoxontest", "wilcoxontest_BH",
                       "ranger_impu_cor", 
                       "mrmr30", "mrmr50", "mrmr75", "mrmr100"),
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
       norm = "norm_quantile")      
    

)  