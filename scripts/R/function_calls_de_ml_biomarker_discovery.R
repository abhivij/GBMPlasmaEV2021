source('scripts/R/de_ml_biomarker_discovery.R')



get_RFE_selected_features(ranked_feature_file_path = "DE_results_2024/tra_result_PREOPEVsHC_agg.csv",
                          RFE_results_file_path = "DE_results_2024/RFE_results_tra_PREOPEVsHC.csv",
                          stop_iter = 202,
                          output_file_path = "DE_results_2024/RFE_features_tra_PREOPEVsHC.csv")
get_RFE_selected_features(ranked_feature_file_path = "DE_results_2024/tra_result_PREOPEVsMET_agg.csv",
                          RFE_results_file_path = "DE_results_2024/RFE_results_tra_PREOPEVsMET.csv",
                          stop_iter = 207,
                          output_file_path = "DE_results_2024/RFE_features_tra_PREOPEVsMET.csv")
get_RFE_selected_features(ranked_feature_file_path = "DE_results_2024/tra_result_METVsHC_agg.csv",
                          RFE_results_file_path = "DE_results_2024/RFE_results_tra_METVsHC.csv",
                          stop_iter = 207,
                          output_file_path = "DE_results_2024/RFE_features_tra_METVsHC.csv")
