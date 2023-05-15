library(tidyverse)
library(sva)
library(readxl)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

phenotype_file_path_prot <- "Data/proteomic_phenotype_combined.txt"
phenotype_file_path_tra <- "Data/transcriptomic_phenotype_combined.txt"

label_prot <- read.table(phenotype_file_path_prot, header=TRUE, sep="\t") %>%
  mutate(condition = case_when(!is.na(POSTOPE_TPVsREC_TP) ~ POSTOPE_TPVsREC_TP,
                               !is.na(PREOPEVsPOSTOPE_TP) ~ PREOPEVsPOSTOPE_TP,
                               TRUE ~ PREOPEVsREC_TP)) %>%
  mutate(condition_cohort = paste0(data_cohort, "_", condition))
label_tra <- read.table(phenotype_file_path_tra, header=TRUE, sep="\t") %>%
  mutate(condition = case_when(!is.na(POSTOPE_TPVsREC_TP) ~ POSTOPE_TPVsREC_TP,
                               !is.na(PREOPEVsPOSTOPE_TP) ~ PREOPEVsPOSTOPE_TP,
                               TRUE ~ PREOPEVsREC_TP)) %>%
  mutate(condition_cohort = paste0(data_cohort, "_", condition))


sum(!is.na(label_prot$condition))
sum(!is.na(label_tra$condition))

summary(factor(label_prot$condition_cohort))
summary(factor(label_tra$condition_cohort))



#########
#prot
data_file_path <- "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv"
validation_data_file_path <- "Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil.csv"

data <- read.csv(data_file_path, row.names = 1)
validation_data <- read.csv(validation_data_file_path, row.names = 1)
colnames(validation_data)[colnames(validation_data) == "SB12_01"] = "SB12"
#use SB22.02
colnames(validation_data)[colnames(validation_data) == "SB22.02"] = "SBtobeused22"
colnames(validation_data)[colnames(validation_data) == "SB22"] = "SB22_dont_include"
colnames(validation_data)[colnames(validation_data) == "SBtobeused22"] = "SB22"
validation_data <- validation_data %>% dplyr::select(-c("SB7")) 

length(intersect(rownames(data), rownames(validation_data)))
#4117

#tra
data_file_path <- "Data/RNA/umi_counts_initial_cohort.csv"
validation_data_file_path <- "Data/RNA/umi_counts_validation_cohort.csv"      

data <- read.csv(data_file_path, row.names = 1)
validation_data <- read.csv(validation_data_file_path, row.names = 1)

length(intersect(rownames(data), rownames(validation_data)))
#272
