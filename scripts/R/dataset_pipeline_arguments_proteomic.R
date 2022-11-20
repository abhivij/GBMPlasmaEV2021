dataset_pipeline_arguments_proteomic <- list(
  
  #with no norm
  #1
  #proteomic PREOPEVsPOSTOPE_TP
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBM_initial_proteomic_impute50fil_no_norm",
       classification_criteria = "PREOPEVsPOSTOPE_TP",
       classes = c("POSTOPE_TP", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_pr",
       norm = "none",
       fems_to_run = c("all", 
                       "t-test", "wilcoxontest",
                       "ranger_pos_impu_cor",
                       "mrmr30", "mrmr50", 
                       "mrmr75", "mrmr100")),
  
  #2
  #proteomic PREOPEVsPOSTOPE_TP
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBM_initial_proteomic_impute50fil_no_norm",
       classification_criteria = "PREOPEVsPOSTOPE_TP",
       classes = c("POSTOPE_TP", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_pr",
       norm = "none",
       fems_to_run = c("mrmr_perc50")),
  
  #3
  #proteomic PREOPEVsPOSTOPE_TP
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBM_initial_proteomic_impute50fil_no_norm",
       classification_criteria = "PREOPEVsPOSTOPE_TP",
       classes = c("POSTOPE_TP", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_pr",
       norm = "none",
       fems_to_run = c("RF_RFE")),
  
  #4
  #proteomic PREOPEVsPOSTOPE_TP
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBM_initial_proteomic_impute50fil_no_norm",
       classification_criteria = "PREOPEVsPOSTOPE_TP",
       classes = c("POSTOPE_TP", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_pr",
       norm = "none",
       fems_to_run = c("ga_rf")),  
  
  
  
  #5
  #proteomic POSTOPE_TPVsREC_TP
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBM_initial_proteomic_impute50fil_no_norm",
       classification_criteria = "POSTOPE_TPVsREC_TP",
       classes = c("REC_TP", "POSTOPE_TP"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_pr",
       norm = "none",
       fems_to_run = c("all", 
                       "t-test", "wilcoxontest",
                       "ranger_pos_impu_cor",
                       "mrmr30", "mrmr50", 
                       "mrmr75", "mrmr100")),
  
  #6
  #proteomic POSTOPE_TPVsREC_TP
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBM_initial_proteomic_impute50fil_no_norm",
       classification_criteria = "POSTOPE_TPVsREC_TP",
       classes = c("REC_TP", "POSTOPE_TP"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_pr",
       norm = "none",
       fems_to_run = c("mrmr_perc50")),
  
  #7
  #proteomic POSTOPE_TPVsREC_TP
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBM_initial_proteomic_impute50fil_no_norm",
       classification_criteria = "POSTOPE_TPVsREC_TP",
       classes = c("REC_TP", "POSTOPE_TP"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_pr",
       norm = "none",
       fems_to_run = c("RF_RFE")),
  
  #8
  #proteomic POSTOPE_TPVsREC_TP
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBM_initial_proteomic_impute50fil_no_norm",
       classification_criteria = "POSTOPE_TPVsREC_TP",
       classes = c("REC_TP", "POSTOPE_TP"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_pr",
       norm = "none",
       fems_to_run = c("ga_rf")),    
  
  
  #9
  #proteomic PREOPEVsREC_TP
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBM_initial_proteomic_impute50fil_no_norm",
       classification_criteria = "PREOPEVsREC_TP",
       classes = c("REC_TP", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_pr",
       norm = "none",
       fems_to_run = c("all", 
                       "t-test", "wilcoxontest",
                       "ranger_pos_impu_cor",
                       "mrmr30", "mrmr50", 
                       "mrmr75", "mrmr100")),
  
  #10
  #proteomic PREOPEVsREC_TP
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBM_initial_proteomic_impute50fil_no_norm",
       classification_criteria = "PREOPEVsREC_TP",
       classes = c("REC_TP", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_pr",
       norm = "none",
       fems_to_run = c("mrmr_perc50")),
  
  #11
  #proteomic PREOPEVsREC_TP
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBM_initial_proteomic_impute50fil_no_norm",
       classification_criteria = "PREOPEVsREC_TP",
       classes = c("REC_TP", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_pr",
       norm = "none",
       fems_to_run = c("RF_RFE")),
  
  #12
  #proteomic PREOPEVsREC_TP
  list(phenotype_file_name = "Data/proteomic_phenotype.txt",
       read_count_dir_path = "Data/Protein/formatted_data",
       read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
       sep = ",",
       dataset_id = "GBM_initial_proteomic_impute50fil_no_norm",
       classification_criteria = "PREOPEVsREC_TP",
       classes = c("REC_TP", "PREOPE"),
       cores = 16,
       results_dir_path = "fem_pipeline_results_pr",
       norm = "none",
       fems_to_run = c("ga_rf")),
  
#with quantile train param norm

#1
#proteomic PREOPEVsPOSTOPE_TP
list(phenotype_file_name = "Data/proteomic_phenotype.txt",
     read_count_dir_path = "Data/Protein/formatted_data",
     read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
     sep = ",",
     dataset_id = "GBM_initial_proteomic_impute50fil_quantile_train_param",
     classification_criteria = "PREOPEVsPOSTOPE_TP",
     classes = c("POSTOPE_TP", "PREOPE"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_pr",
     norm = "quantile_train_param",
     fems_to_run = c("all", 
                     "t-test", "wilcoxontest",
                     "ranger_pos_impu_cor",
                     "mrmr30", "mrmr50", 
                     "mrmr75", "mrmr100")),

#2
#proteomic PREOPEVsPOSTOPE_TP
list(phenotype_file_name = "Data/proteomic_phenotype.txt",
     read_count_dir_path = "Data/Protein/formatted_data",
     read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
     sep = ",",
     dataset_id = "GBM_initial_proteomic_impute50fil_quantile_train_param",
     classification_criteria = "PREOPEVsPOSTOPE_TP",
     classes = c("POSTOPE_TP", "PREOPE"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_pr",
     norm = "quantile_train_param",
     fems_to_run = c("mrmr_perc50")),

#3
#proteomic PREOPEVsPOSTOPE_TP
list(phenotype_file_name = "Data/proteomic_phenotype.txt",
     read_count_dir_path = "Data/Protein/formatted_data",
     read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
     sep = ",",
     dataset_id = "GBM_initial_proteomic_impute50fil_quantile_train_param",
     classification_criteria = "PREOPEVsPOSTOPE_TP",
     classes = c("POSTOPE_TP", "PREOPE"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_pr",
     norm = "quantile_train_param",
     fems_to_run = c("RF_RFE")),

#4
#proteomic PREOPEVsPOSTOPE_TP
list(phenotype_file_name = "Data/proteomic_phenotype.txt",
     read_count_dir_path = "Data/Protein/formatted_data",
     read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
     sep = ",",
     dataset_id = "GBM_initial_proteomic_impute50fil_quantile_train_param",
     classification_criteria = "PREOPEVsPOSTOPE_TP",
     classes = c("POSTOPE_TP", "PREOPE"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_pr",
     norm = "quantile_train_param",
     fems_to_run = c("ga_rf")),  



#5
#proteomic POSTOPE_TPVsREC_TP
list(phenotype_file_name = "Data/proteomic_phenotype.txt",
     read_count_dir_path = "Data/Protein/formatted_data",
     read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
     sep = ",",
     dataset_id = "GBM_initial_proteomic_impute50fil_quantile_train_param",
     classification_criteria = "POSTOPE_TPVsREC_TP",
     classes = c("REC_TP", "POSTOPE_TP"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_pr",
     norm = "quantile_train_param",
     fems_to_run = c("all", 
                     "t-test", "wilcoxontest",
                     "ranger_pos_impu_cor",
                     "mrmr30", "mrmr50", 
                     "mrmr75", "mrmr100")),

#6
#proteomic POSTOPE_TPVsREC_TP
list(phenotype_file_name = "Data/proteomic_phenotype.txt",
     read_count_dir_path = "Data/Protein/formatted_data",
     read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
     sep = ",",
     dataset_id = "GBM_initial_proteomic_impute50fil_quantile_train_param",
     classification_criteria = "POSTOPE_TPVsREC_TP",
     classes = c("REC_TP", "POSTOPE_TP"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_pr",
     norm = "quantile_train_param",
     fems_to_run = c("mrmr_perc50")),

#7
#proteomic POSTOPE_TPVsREC_TP
list(phenotype_file_name = "Data/proteomic_phenotype.txt",
     read_count_dir_path = "Data/Protein/formatted_data",
     read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
     sep = ",",
     dataset_id = "GBM_initial_proteomic_impute50fil_quantile_train_param",
     classification_criteria = "POSTOPE_TPVsREC_TP",
     classes = c("REC_TP", "POSTOPE_TP"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_pr",
     norm = "quantile_train_param",
     fems_to_run = c("RF_RFE")),

#8
#proteomic POSTOPE_TPVsREC_TP
list(phenotype_file_name = "Data/proteomic_phenotype.txt",
     read_count_dir_path = "Data/Protein/formatted_data",
     read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
     sep = ",",
     dataset_id = "GBM_initial_proteomic_impute50fil_quantile_train_param",
     classification_criteria = "POSTOPE_TPVsREC_TP",
     classes = c("REC_TP", "POSTOPE_TP"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_pr",
     norm = "quantile_train_param",
     fems_to_run = c("ga_rf")),    


#9
#proteomic PREOPEVsREC_TP
list(phenotype_file_name = "Data/proteomic_phenotype.txt",
     read_count_dir_path = "Data/Protein/formatted_data",
     read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
     sep = ",",
     dataset_id = "GBM_initial_proteomic_impute50fil_quantile_train_param",
     classification_criteria = "PREOPEVsREC_TP",
     classes = c("REC_TP", "PREOPE"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_pr",
     norm = "quantile_train_param",
     fems_to_run = c("all", 
                     "t-test", "wilcoxontest",
                     "ranger_pos_impu_cor",
                     "mrmr30", "mrmr50", 
                     "mrmr75", "mrmr100")),

#10
#proteomic PREOPEVsREC_TP
list(phenotype_file_name = "Data/proteomic_phenotype.txt",
     read_count_dir_path = "Data/Protein/formatted_data",
     read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
     sep = ",",
     dataset_id = "GBM_initial_proteomic_impute50fil_quantile_train_param",
     classification_criteria = "PREOPEVsREC_TP",
     classes = c("REC_TP", "PREOPE"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_pr",
     norm = "quantile_train_param",
     fems_to_run = c("mrmr_perc50")),

#11
#proteomic PREOPEVsREC_TP
list(phenotype_file_name = "Data/proteomic_phenotype.txt",
     read_count_dir_path = "Data/Protein/formatted_data",
     read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
     sep = ",",
     dataset_id = "GBM_initial_proteomic_impute50fil_quantile_train_param",
     classification_criteria = "PREOPEVsREC_TP",
     classes = c("REC_TP", "PREOPE"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_pr",
     norm = "quantile_train_param",
     fems_to_run = c("RF_RFE")),

#12
#proteomic PREOPEVsREC_TP
list(phenotype_file_name = "Data/proteomic_phenotype.txt",
     read_count_dir_path = "Data/Protein/formatted_data",
     read_count_file_name = "Q1-6_nonorm_formatted_impute50fil.csv",
     sep = ",",
     dataset_id = "GBM_initial_proteomic_impute50fil_quantile_train_param",
     classification_criteria = "PREOPEVsREC_TP",
     classes = c("REC_TP", "PREOPE"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_pr",
     norm = "quantile_train_param",
     fems_to_run = c("ga_rf"))

)
