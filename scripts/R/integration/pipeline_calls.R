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
