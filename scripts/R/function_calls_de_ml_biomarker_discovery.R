source('scripts/R/de_ml_biomarker_discovery.R')

get_RFE_selected_features(ranked_feature_file_path = "DE_results_2024/tra_result_PREOPEVsMET_agg.csv",
                          RFE_results_001_file_path = "DE_results_2024/RFE_results_tra_001_PREOPEVsMET.csv",
                          RFE_results_003_file_path = "DE_results_2024/RFE_results_tra_003_PREOPEVsMET.csv",
                          RFE_results_005_file_path = "DE_results_2024/RFE_results_tra_005_PREOPEVsMET.csv",
                          output_file_path = "DE_results_2024/features/features_RFE_tra_PREOPEVsMET.csv")
get_RFE_selected_features(ranked_feature_file_path = "DE_results_2024/tra_result_PREOPEVsHC_agg.csv",
                          RFE_results_001_file_path = "DE_results_2024/RFE_results_tra_001_PREOPEVsHC.csv",
                          RFE_results_003_file_path = "DE_results_2024/RFE_results_tra_003_PREOPEVsHC.csv",
                          RFE_results_005_file_path = "DE_results_2024/RFE_results_tra_005_PREOPEVsHC.csv",
                          output_file_path = "DE_results_2024/features/features_RFE_tra_PREOPEVsHC.csv")
get_RFE_selected_features(ranked_feature_file_path = "DE_results_2024/tra_result_METVsHC_agg.csv",
                          RFE_results_001_file_path = "DE_results_2024/RFE_results_tra_001_METVsHC.csv",
                          RFE_results_003_file_path = "DE_results_2024/RFE_results_tra_003_METVsHC.csv",
                          RFE_results_005_file_path = "DE_results_2024/RFE_results_tra_005_METVsHC.csv",
                          output_file_path = "DE_results_2024/features/features_RFE_tra_METVsHC.csv")

get_RFE_selected_features(ranked_feature_file_path = "DE_results_2024/prot_result_PREOPEVsMET_agg.csv",
                          RFE_results_001_file_path = "DE_results_2024/RFE_results_prot_filtered_001_PREOPEVsMET.csv",
                          RFE_results_003_file_path = "DE_results_2024/RFE_results_prot_filtered_003_PREOPEVsMET.csv",
                          RFE_results_005_file_path = "DE_results_2024/RFE_results_prot_filtered_005_PREOPEVsMET.csv",
                          output_file_path = "DE_results_2024/features/features_RFE_prot_PREOPEVsMET.csv", 
                          filter_low_scored_features = TRUE)
get_RFE_selected_features(ranked_feature_file_path = "DE_results_2024/prot_result_PREOPEVsHC_agg.csv",
                          RFE_results_001_file_path = "DE_results_2024/RFE_results_prot_filtered_001_PREOPEVsHC.csv",
                          RFE_results_003_file_path = "DE_results_2024/RFE_results_prot_filtered_003_PREOPEVsHC.csv",
                          RFE_results_005_file_path = "DE_results_2024/RFE_results_prot_filtered_005_PREOPEVsHC.csv",
                          output_file_path = "DE_results_2024/features/features_RFE_prot_PREOPEVsHC.csv", 
                          filter_low_scored_features = TRUE)
get_RFE_selected_features(ranked_feature_file_path = "DE_results_2024/prot_result_METVsHC_agg.csv",
                          RFE_results_001_file_path = "DE_results_2024/RFE_results_prot_filtered_001_METVsHC.csv",
                          RFE_results_003_file_path = "DE_results_2024/RFE_results_prot_filtered_003_METVsHC.csv",
                          RFE_results_005_file_path = "DE_results_2024/RFE_results_prot_filtered_005_METVsHC.csv",
                          output_file_path = "DE_results_2024/features/features_RFE_prot_METVsHC.csv", 
                          filter_low_scored_features = TRUE)
