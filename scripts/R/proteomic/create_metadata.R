library(tidyverse)
library(readxl)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

source("scripts/R/utils.R")

#get group mappings from Q1-6 proteomics output
protein_data <- read.csv(file = "Data/Protein/norm_output/norm_annotatedQ1-6_NA_equalizeMedians.csv")

group_mapping <- protein_data %>%
  select(SUBJECT_ORIGINAL, GROUP_ORIGINAL) %>%
  arrange(SUBJECT_ORIGINAL)
rm(protein_data)

metadata <- group_mapping %>%
  rename(GROUP_Q1to6 = GROUP_ORIGINAL)

#get group mappings from Q7 proteomics output
protein_data <- read.csv(file = "Data/Protein/norm_output/norm_annotatedQ7_NA_equalizeMedians.csv")

group_mapping <- protein_data %>%
  select(SUBJECT_ORIGINAL, GROUP_ORIGINAL) %>%
  arrange(SUBJECT_ORIGINAL)
rm(protein_data)

metadata <- metadata %>%
  inner_join(group_mapping) %>%
  rename(GROUP_Q7 = GROUP_ORIGINAL)

metadata <- metadata %>%
  mutate(GROUP_Q1to6 = gsub("-", "_", GROUP_Q1to6)) %>%
  mutate(GROUP_Q7 = gsub("-", "_", GROUP_Q7))

#order metadata as HC1, HC2, HC3, ...
#   and NOT HC1, HC10, HC11 ... HC19, HC2, HC20, ...
# to get this ordering split string part and num part, and order separately

metadata <- metadata %>%
  separate(SUBJECT_ORIGINAL, c("modified_sample", NA), remove = FALSE, sep = "-", fill = "right") %>%
  mutate(str_part = gsub("[0-9]", "", modified_sample)) %>%
  mutate(num_part = strtoi(gsub("[^0-9]", "", modified_sample), base = 10)) 
metadata <- metadata %>%
  arrange(str_part, num_part) %>%
  select(-c(str_part, num_part, modified_sample))

write.csv(metadata, "Data/proteomic_sample_metadata.csv", row.names = FALSE)


#create phenotype file to use in FEMPipeline
phenotype_info <- metadata %>%
  rename("Sample" = "SUBJECT_ORIGINAL") %>%
  mutate(Biomarker = "Protein", .after = "Sample") %>%  
  mutate(Technology = "SWATH-MS", .after = "Biomarker") 

phenotype_info <- insert_comparison_columns(phenotype_info,
                                            comparison_list = list(c("PREOPE", "MET"), 
                                                                   c("PREOPE", "HC"), 
                                                                   c("MET", "HC"),
                                                                   
                                                                   c("PREOPE", "POSTOPE_T"), 
                                                                   c("PREOPE", "POSTOPE_P"), 
                                                                   c("POSTOPE_T", "POSTOPE_P"),
                                                                   
                                                                   c("POSTOPE_T", "REC_T"),
                                                                   c("POSTOPE_P", "REC_P"),
                                                                   c("POSTOPE_T", "PREREC")
                                                                   ), 
                                            class_column_name = "GROUP_Q1to6")
phenotype_info <- insert_comparison_columns(phenotype_info,
                                            comparison_list = list(c("PREOPE", "REC_TP")), 
                                            class_column_name = "GROUP_Q7")
phenotype_info <- insert_comparison_columns(phenotype_info,
                                            comparison_list = list(c("POSTOPE_T", "HC"),
                                                                   c("POSTOPE_P", "HC")
                                            ), 
                                            class_column_name = "GROUP_Q1to6")
phenotype_info <- phenotype_info %>%   
  mutate(POSTOPE_TPVsHC = case_when(
    GROUP_Q1to6 == "POSTOPE_T" | GROUP_Q1to6 == "POSTOPE_P" ~ "POSTOPE_TP",
    GROUP_Q1to6 == "HC" ~ "HC",
    TRUE ~ NA_character_)) %>%
  mutate(PREOPEVsPOSTOPE_TP = case_when(
    GROUP_Q1to6 == "PREOPE" ~ "PREOPE",
    GROUP_Q1to6 == "POSTOPE_T" | GROUP_Q1to6 == "POSTOPE_P" ~ "POSTOPE_TP",
    TRUE ~ NA_character_)) %>%
  mutate(POSTOPE_TPVsREC_TP = case_when(
    #not using GROUP_Q7 because POSTOPE_T, POSTOPE_P are not properly marked in it
    GROUP_Q1to6 == "POSTOPE_T" | GROUP_Q1to6 == "POSTOPE_P" ~ "POSTOPE_TP",
    GROUP_Q1to6 == "REC_T" | GROUP_Q1to6 == "REC_P" ~ "REC_TP",
    TRUE ~ NA_character_))  

phenotype_info <- phenotype_info %>%
  mutate(HC_MET_PREOPE_REC_TP = case_when(
    GROUP_Q1to6 == "HC" ~ "HC",
    GROUP_Q1to6 == "MET" ~ "MET",
    GROUP_Q1to6 == "PREOPE" ~ "PREOPE",
    GROUP_Q1to6 == "REC_T" | GROUP_Q1to6 == "REC_P" ~ "REC_TP",
    TRUE ~ NA_character_)) %>%
  mutate(PREOPE_POSTOPE_TP_PREREC_REC_TP = case_when(
    GROUP_Q1to6 == "PREOPE" ~ "PREOPE",
    GROUP_Q1to6 == "POSTOPE_T" | GROUP_Q1to6 == "POSTOPE_P" ~ "POSTOPE_TP",
    GROUP_Q1to6 == "PREREC" ~ "PREREC",
    GROUP_Q1to6 == "REC_T" | GROUP_Q1to6 == "REC_P" ~ "REC_TP",
    TRUE ~ NA_character_)) %>%
  mutate(PREOPE_POSTOPE_T_POSTOPE_P_PREREC_REC_TP = case_when(
    GROUP_Q1to6 == "PREOPE" ~ "PREOPE",
    GROUP_Q1to6 == "POSTOPE_T" ~ "POSTOPE_T",
    GROUP_Q1to6 == "POSTOPE_P" ~ "POSTOPE_P",
    GROUP_Q1to6 == "PREREC" ~ "PREREC",
    GROUP_Q1to6 == "REC_T" | GROUP_Q1to6 == "REC_P" ~ "REC_TP",
    TRUE ~ NA_character_))

write.table(phenotype_info, 
            file = "Data/proteomic_phenotype.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)


# test if read using FEM pipeline works

# phenotype_file <- "Data/proteomic_phenotype.txt"
# classification_criteria <- "PREOPEVsMET"
# filter <- expression(TRUE)
# phenotype <- read.table(phenotype_file, header=TRUE, sep="\t")
# all_equal(phenotype, phenotype_info)
# 
# data <- read.table("Data/Protein/formatted_data/Q1-6_nonorm_formatted.csv", header=TRUE, sep=",", 
#                    row.names=1, skip=0,
#                    nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
# 
# extracted_samples <- phenotype %>% subset(eval(filter))
# extracted_samples <- extracted_samples[!is.na(extracted_samples[classification_criteria]), ]
# extracted_samples$Sample <- factor(extracted_samples$Sample)
# 
# filtered_samples_read_count <- data %>% dplyr::select(extracted_samples$Sample)
# 
# filtered_samples_output_labels <- extracted_samples[, c('Sample', classification_criteria)]
# colnames(filtered_samples_output_labels) <- c("Sample", "Label")




#############
#creating phenotype for validation cohort and combined cohort

phenotype <- read.table("Data/proteomic_phenotype.txt", header=TRUE, sep="\t")
validation_metadata <- read.csv("Data/RNA_validation/metadata_glionet.csv")

validation_phenotype <- validation_metadata %>%
  rename("Sample" = "sample_id") %>%
  relocate(Sample, .before = patient_id) %>%
  dplyr::select(-c(sample_instance)) %>%
  mutate("POSTOPE_TPVsREC_TP" = case_when(category_old_name == "POSTOPE_TP" ~ "POSTOPE_TP",
                                          category_old_name == "REC_TP" ~ "REC_TP",
                                          TRUE ~ NA_character_),
         "PREOPEVsPOSTOPE_TP" = case_when(category_old_name == "PREOPE" ~ "PREOPE",
                                          category_old_name == "POSTOPE_TP" ~ "POSTOPE_TP",
                                          TRUE ~ NA_character_),
         "PREOPEVsREC_TP" = case_when(category_old_name == "PREOPE" ~ "PREOPE",
                                      category_old_name == "REC_TP" ~ "REC_TP",
                                      TRUE ~ NA_character_)
  ) %>%
  filter(Sample != "SB7")

write.table(validation_phenotype, 
            file = "Data/proteomic_phenotype_validation.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)

combined_phenotype <- rbind(phenotype %>% 
                              dplyr::select(Sample, POSTOPE_TPVsREC_TP, 
                                            PREOPEVsPOSTOPE_TP, PREOPEVsREC_TP) %>%
                              mutate(data_cohort = "initial"),
                            validation_phenotype %>% 
                              dplyr::select(Sample, POSTOPE_TPVsREC_TP, 
                                            PREOPEVsPOSTOPE_TP, PREOPEVsREC_TP) %>%
                              mutate(data_cohort = "validation"))
write.table(combined_phenotype, 
            file = "Data/proteomic_phenotype_combined.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)

################################################################
# # test if read using FEM pipeline works validation cohort
# 
# phenotype_file <- "Data/proteomic_phenotype_validation.txt"
# classification_criteria <- "POSTOPE_TPVsREC_TP"
# filter <- expression(TRUE)
# phenotype <- read.table(phenotype_file, header=TRUE, sep="\t")
# 
# data <- read.table("Data/Protein/formatted_data/newcohort_common_correctedsamples.csv", header=TRUE, sep=",",
#                    row.names=1, skip=0,
#                    nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
# 
# extracted_samples <- phenotype %>% subset(eval(filter))
# extracted_samples <- extracted_samples[!is.na(extracted_samples[classification_criteria]), ]
# extracted_samples$Sample <- factor(extracted_samples$Sample)
# 
# filtered_samples_read_count <- data %>% dplyr::select(extracted_samples$Sample)
# 
# filtered_samples_output_labels <- extracted_samples[, c('Sample', classification_criteria)]
# colnames(filtered_samples_output_labels) <- c("Sample", "Label")
# 
################################################################
# # test if read using FEM pipeline works validation cohort post combat
# 
# phenotype_file <- "Data/proteomic_phenotype_validation.txt"
# classification_criteria <- "POSTOPE_TPVsREC_TP"
# filter <- expression(TRUE)
# phenotype <- read.table(phenotype_file, header=TRUE, sep="\t")
# 
# data <- read.table("Data/Protein/validation_data.combat.POSTOPE_TPVsREC_TP.csv", header=TRUE, sep=",",
#                    row.names=1, skip=0,
#                    nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
# 
# extracted_samples <- phenotype %>% subset(eval(filter))
# extracted_samples <- extracted_samples[!is.na(extracted_samples[classification_criteria]), ]
# extracted_samples$Sample <- factor(extracted_samples$Sample)
# 
# filtered_samples_read_count <- data %>% dplyr::select(extracted_samples$Sample)
# 
# filtered_samples_output_labels <- extracted_samples[, c('Sample', classification_criteria)]
# colnames(filtered_samples_output_labels) <- c("Sample", "Label")
# 
# 
################################################################
# # test if read using FEM pipeline works combined cohort
# 
# phenotype_file <- "Data/proteomic_phenotype_combined.txt"
# classification_criteria <- "POSTOPE_TPVsREC_TP"
# filter <- expression(TRUE)
# phenotype <- read.table(phenotype_file, header=TRUE, sep="\t")
# 
# data <- read.table("Data/Protein/combined_data.combat.POSTOPE_TPVsREC_TP.csv", header=TRUE, sep=",",
#                    row.names=1, skip=0,
#                    nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
# 
# extracted_samples <- phenotype %>% subset(eval(filter))
# extracted_samples <- extracted_samples[!is.na(extracted_samples[classification_criteria]), ]
# extracted_samples$Sample <- factor(extracted_samples$Sample)
# 
# filtered_samples_read_count <- data %>% dplyr::select(extracted_samples$Sample)
# 
# filtered_samples_output_labels <- extracted_samples[, c('Sample', classification_criteria)]
# colnames(filtered_samples_output_labels) <- c("Sample", "Label")


###############################

#creating combined phenotype file for PREOPE Vs MET Vs HC

phenotype <- read.table("Data/proteomic_phenotype.txt", header=TRUE, sep="\t")
validation_metadata <- read.csv("Data/RNA_validation/metadata_glionet.csv")

validation_phenotype <- validation_metadata %>%
  rename("Sample" = "sample_id") %>%
  relocate(Sample, .before = patient_id) %>%
  dplyr::select(-c(sample_instance))

validation_phenotype <- insert_comparison_columns(validation_phenotype,
                                                  comparison_list = list(c("PREOPE", "MET"), 
                                                                         c("PREOPE", "HC"), 
                                                                         c("MET", "HC")), 
                                                  class_column_name = "category_old_name")


combined_phenotype <- rbind(phenotype %>% 
                              dplyr::select(Sample, PREOPEVsMET, PREOPEVsHC, METVsHC) %>%
                              mutate(data_cohort = "initial"),
                            validation_phenotype %>% 
                              dplyr::select(Sample, PREOPEVsMET, PREOPEVsHC, METVsHC) %>%
                              mutate(data_cohort = "validation")) %>%
  filter(!(is.na(PREOPEVsMET) & is.na(PREOPEVsHC) & is.na(METVsHC)))

PREOPE_MET_HC_metadata <- read.csv("Data/PREOPE_MET_HC/meta_data_updated.csv")

PREOPE_MET_HC_metadata <- PREOPE_MET_HC_metadata %>%
  full_join(combined_phenotype) %>%
  mutate(Subgroup = gsub(" ", "_", Subgroup, fixed = TRUE))

PREOPE_MET_HC_phenotype <- insert_comparison_columns(PREOPE_MET_HC_metadata,
                                                     comparison_list = list(c("Melanoma_met", "Other")), 
                                                     class_column_name = "Subgroup")

write.table(PREOPE_MET_HC_phenotype, 
            file = "Data/proteomic_phenotype_PREOPE_MET_HC.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)
