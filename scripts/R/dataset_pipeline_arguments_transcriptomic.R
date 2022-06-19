dataset_pipeline_arguments_transcriptomic <- list(
  
  #1
  #transcriptomic PREOPEVsPOSTOPE_TP
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts_initial_cohort.csv",
       sep = ",",
       dataset_id = "GBM_tr_initial",
       classification_criteria = "PREOPEVsPOSTOPE_TP",
       classes = c("POSTOPE_TP", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_tr",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("all", 
                       "t-test", "wilcoxontest",
                       "ranger_impu_cor",
                       "ranger_pos_impu_cor",
                       "mrmr10", "mrmr20",
                       "mrmr30", "mrmr50", 
                       "mrmr75", "mrmr100")),
  
  #2
  #transcriptomic PREOPEVsPOSTOPE_TP
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts_initial_cohort.csv",
       sep = ",",
       dataset_id = "GBM_tr_initial",
       classification_criteria = "PREOPEVsPOSTOPE_TP",
       classes = c("POSTOPE_TP", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_tr",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("mrmr_perc50")),
  
  #3
  #transcriptomic PREOPEVsPOSTOPE_TP
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts_initial_cohort.csv",
       sep = ",",
       dataset_id = "GBM_tr_initial",
       classification_criteria = "PREOPEVsPOSTOPE_TP",
       classes = c("POSTOPE_TP", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_tr",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("RF_RFE")),
  
  #4
  #transcriptomic PREOPEVsPOSTOPE_TP
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts_initial_cohort.csv",
       sep = ",",
       dataset_id = "GBM_tr_initial",
       classification_criteria = "PREOPEVsPOSTOPE_TP",
       classes = c("POSTOPE_TP", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_tr",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("ga_rf")),  
  
  
  
  #5
  #transcriptomic POSTOPE_TPVsREC_TP
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts_initial_cohort.csv",
       sep = ",",
       dataset_id = "GBM_tr_initial",
       classification_criteria = "POSTOPE_TPVsREC_TP",
       classes = c("REC_TP", "POSTOPE_TP"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_tr",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("all", 
                       "t-test", "wilcoxontest",
                       "ranger_impu_cor",
                       "ranger_pos_impu_cor",
                       "mrmr10", "mrmr20",
                       "mrmr30", "mrmr50", 
                       "mrmr75", "mrmr100")),
  
  #6
  #transcriptomic POSTOPE_TPVsREC_TP
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts_initial_cohort.csv",
       sep = ",",
       dataset_id = "GBM_tr_initial",
       classification_criteria = "POSTOPE_TPVsREC_TP",
       classes = c("REC_TP", "POSTOPE_TP"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_tr",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("mrmr_perc50")),
  
  #7
  #transcriptomic POSTOPE_TPVsREC_TP
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts_initial_cohort.csv",
       sep = ",",
       dataset_id = "GBM_tr_initial",
       classification_criteria = "POSTOPE_TPVsREC_TP",
       classes = c("REC_TP", "POSTOPE_TP"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_tr",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("RF_RFE")),
  
  #8
  #transcriptomic POSTOPE_TPVsREC_TP
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts_initial_cohort.csv",
       sep = ",",
       dataset_id = "GBM_tr_initial",
       classification_criteria = "POSTOPE_TPVsREC_TP",
       classes = c("REC_TP", "POSTOPE_TP"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_tr",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("ga_rf")),    
  
  
  #9
  #transcriptomic PREOPEVsREC_TP
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts_initial_cohort.csv",
       sep = ",",
       dataset_id = "GBM_tr_initial",
       classification_criteria = "PREOPEVsREC_TP",
       classes = c("REC_TP", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_tr",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("all", 
                       "t-test", "wilcoxontest",
                       "ranger_impu_cor",
                       "ranger_pos_impu_cor",
                       "mrmr10", "mrmr20",
                       "mrmr30", "mrmr50", 
                       "mrmr75", "mrmr100")),
  
  #10
  #transcriptomic PREOPEVsREC_TP
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts_initial_cohort.csv",
       sep = ",",
       dataset_id = "GBM_tr_initial",
       classification_criteria = "PREOPEVsREC_TP",
       classes = c("REC_TP", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_tr",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("mrmr_perc50")),
  
  #11
  #transcriptomic PREOPEVsREC_TP
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts_initial_cohort.csv",
       sep = ",",
       dataset_id = "GBM_tr_initial",
       classification_criteria = "PREOPEVsREC_TP",
       classes = c("REC_TP", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_tr",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("RF_RFE")),
  
  #12
  #transcriptomic PREOPEVsREC_TP
  list(phenotype_file_name = "Data/transcriptomic_phenotype.txt",
       read_count_dir_path = "Data/RNA",
       read_count_file_name = "umi_counts_initial_cohort.csv",
       sep = ",",
       dataset_id = "GBM_tr_initial",
       classification_criteria = "PREOPEVsREC_TP",
       classes = c("REC_TP", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_tr",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("ga_rf"))
  
)
