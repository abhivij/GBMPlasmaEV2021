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
       fems_to_run = c("RF_RFE", "ga_rf"))           

)  