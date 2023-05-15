source("scripts/R/integration/base_model_pipeline.R")

base_model_pipeline(comparison = "POSTOPE_TPVsREC_TP",
                    conditions = c("REC_TP", "POSTOPE_TP"),
                    data_file_path_prot <- "Data/Protein/combined_data.combat.POSTOPE_TPVsREC_TP.csv",
                    data_file_path_tra <- "Data/RNA/combined_data.combat.POSTOPE_TPVsREC_TP.csv")

base_model_pipeline(comparison = "PREOPEVsPOSTOPE_TP",
                    conditions = c("POSTOPE_TP", "PREOPE"),
                    data_file_path_prot <- "Data/Protein/combined_data.combat.PREOPEVsPOSTOPE_TP.csv",
                    data_file_path_tra <- "Data/RNA/combined_data.combat.PREOPEVsPOSTOPE_TP.csv")

base_model_pipeline(comparison = "PREOPEVsREC_TP",
                    conditions = c("REC_TP", "PREOPE"),
                    data_file_path_prot <- "Data/Protein/combined_data.combat.PREOPEVsREC_TP.csv",
                    data_file_path_tra <- "Data/RNA/combined_data.combat.PREOPEVsREC_TP.csv")

models <- read.csv("Data/prediction_result/integration/POSTOPE_TPVsREC_TP.csv") %>%
  dplyr::select(model) %>%
  distinct()
models <- models$model
omics_types <- c("prot", "tra")
sample_types <- c("train", "test") 

for(omics_type in omics_types){
  for(model in models){
    for(sample_type in sample_types){
      compute_metrics(comparison = "POSTOPE_TPVsREC_TP",
                      conditions = c("REC_TP", "POSTOPE_TP"),
                      o_type = omics_type,
                      m = model,
                      sample_type = sample_type,
                      result_file_path = "Data/prediction_result/integration/POSTOPE_TPVsREC_TP.csv",
                      metric_output_file_path = "Data/prediction_result/integration/metrics_base_models.csv")
      compute_metrics(comparison = "PREOPEVsPOSTOPE_TP",
                      conditions = c("POSTOPE_TP", "PREOPE"),
                      o_type = omics_type,
                      m = model,
                      sample_type = sample_type,
                      result_file_path = "Data/prediction_result/integration/PREOPEVsPOSTOPE_TP.csv",
                      metric_output_file_path = "Data/prediction_result/integration/metrics_base_models.csv")
      compute_metrics(comparison = "PREOPEVsREC_TP",
                      conditions = c("REC_TP", "PREOPE"),
                      o_type = omics_type,
                      m = model,
                      sample_type = sample_type,
                      result_file_path = "Data/prediction_result/integration/PREOPEVsREC_TP.csv",
                      metric_output_file_path = "Data/prediction_result/integration/metrics_base_models.csv")
    }
  }  
}

metrics <- read.csv("Data/prediction_result/integration/metrics_base_models.csv") %>%
  arrange(Comparison, Omics_type, Type, model)
write.csv(metrics, "Data/prediction_result/integration/metrics_base_models_arranged.csv", row.names = FALSE)



#level 2 stack

for(o in c("prot", "tra", "both")){
  level2_stack(comparison = "POSTOPE_TPVsREC_TP",
               conditions = c("REC_TP", "POSTOPE_TP"),
               data_file_path = "Data/prediction_result/integration/POSTOPE_TPVsREC_TP.csv",
               o_type = o)
  level2_stack(comparison = "PREOPEVsPOSTOPE_TP",
               conditions = c("POSTOPE_TP", "PREOPE"),
               data_file_path = "Data/prediction_result/integration/PREOPEVsPOSTOPE_TP.csv",
               o_type = o)
  level2_stack(comparison = "PREOPEVsREC_TP",
               conditions = c("REC_TP", "PREOPE"),
               data_file_path = "Data/prediction_result/integration/PREOPEVsREC_TP.csv",
               o_type = o)
}

for(o in c("tra", "prot", "both")){
  for(sample_type in sample_types){
    print(paste(o, sample_type))
    compute_metrics(comparison = "POSTOPE_TPVsREC_TP",
                    conditions = c("REC_TP", "POSTOPE_TP"),
                    o_type = paste(o, "stacked", sep = "_"),
                    m = "L2 Regularized logistic regression",
                    sample_type = sample_type,
                    result_file_path = paste0("Data/prediction_result/integration/stacked/POSTOPE_TPVsREC_TP",
                                              "_", o,
                                              ".csv"),
                    metric_output_file_path = "Data/prediction_result/integration/stacked/metrics.csv")
    compute_metrics(comparison = "PREOPEVsPOSTOPE_TP",
                    conditions = c("POSTOPE_TP", "PREOPE"),
                    o_type = paste(o, "stacked", sep = "_"),
                    m = "L2 Regularized logistic regression",
                    sample_type = sample_type,
                    result_file_path = paste0("Data/prediction_result/integration/stacked/PREOPEVsPOSTOPE_TP",
                                              "_", o,
                                              ".csv"),
                    metric_output_file_path = "Data/prediction_result/integration/stacked/metrics.csv")
    compute_metrics(comparison = "PREOPEVsREC_TP",
                    conditions = c("REC_TP", "PREOPE"),
                    o_type = paste(o, "stacked", sep = "_"),
                    m = "L2 Regularized logistic regression",
                    sample_type = sample_type,
                    result_file_path = paste0("Data/prediction_result/integration/stacked/PREOPEVsREC_TP",
                                              "_", o,
                                              ".csv"),
                    metric_output_file_path = "Data/prediction_result/integration/stacked/metrics.csv")
    print(paste(o, sample_type, "done------------------"))
  }
}
