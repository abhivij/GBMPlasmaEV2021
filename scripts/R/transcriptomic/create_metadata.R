library(tidyverse)
library(readxl)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

source("scripts/R/utils.R")

umi_counts <- read.csv("Data/RNA/umi_counts.csv", row.names = 1)

#get group mappings from Q1-6 proteomics output
protein_data <- read.csv(file = "Data/Protein/norm_output/norm_annotatedQ1-6_NA_equalizeMedians.csv")

group_mapping <- protein_data %>%
  select(SUBJECT_ORIGINAL, GROUP_ORIGINAL) %>%
  arrange(SUBJECT_ORIGINAL)
rm(protein_data)

metadata <- data.frame(SUBJECT_ORIGINAL = colnames(umi_counts))

#order metadata as HB1, HB2, HB3, ...
#   and NOT HB1, HB10, HB11 ... HB19, HB2, HB20, ...
# to get this ordering split string part and num part, and order separately
metadata <- metadata %>%
  mutate(str_part = gsub("[0-9]", "", SUBJECT_ORIGINAL)) %>%
  mutate(num_part = strtoi(gsub("[^0-9]", "", SUBJECT_ORIGINAL))) %>%
  arrange(str_part, num_part) %>%
  select(SUBJECT_ORIGINAL)


#in proteomics data, unlike transcriptomics, HB numbered 01, 02, 03, ...
# so changing that to match transcriptomics sample names
group_mapping <- group_mapping %>%
  mutate(SUBJECT_ORIGINAL = gsub("HB0", "HB", SUBJECT_ORIGINAL))

metadata <- metadata %>%
  inner_join(group_mapping) %>%
  rename(GROUP_Q1to6 = GROUP_ORIGINAL)

group_mapping[!group_mapping$SUBJECT_ORIGINAL %in% metadata$SUBJECT_ORIGINAL,]
write.csv(group_mapping[!group_mapping$SUBJECT_ORIGINAL %in% metadata$SUBJECT_ORIGINAL,],
          "Data/missing_in_transcriptomics.csv", row.names = FALSE)

#get group mappings from Q7 proteomics output
protein_data <- read.csv(file = "Data/Protein/norm_output/norm_annotatedQ7_NA_equalizeMedians.csv")

group_mapping <- protein_data %>%
  select(SUBJECT_ORIGINAL, GROUP_ORIGINAL) %>%
  arrange(SUBJECT_ORIGINAL)
rm(protein_data)

group_mapping <- group_mapping %>%
  mutate(SUBJECT_ORIGINAL = gsub("HB0", "HB", SUBJECT_ORIGINAL))

metadata <- metadata %>%
  inner_join(group_mapping) %>%
  rename(GROUP_Q7 = GROUP_ORIGINAL)

group_mapping[!group_mapping$SUBJECT_ORIGINAL %in% metadata$SUBJECT_ORIGINAL,]

metadata <- metadata %>%
  mutate(GROUP_Q1to6 = gsub("-", "_", GROUP_Q1to6)) %>%
  mutate(GROUP_Q7 = gsub("-", "_", GROUP_Q7))

write.csv(metadata, "Data/transcriptomic_sample_metadata.csv", row.names = FALSE)


#note : metadata from below sheet seem to be different
# categories_data <- read_excel("Data/Cohorts for comparative analyses.xlsx", 
#                               sheet = "Categories_RNAseq")
# categories_data <- categories_data %>%
#   rename("SUBJECT_ORIGINAL" = "...1")
# metadata_from_categories_data <- categories_data %>%
#   select(SUBJECT_ORIGINAL, Cohort_code, Cohort_code_2) %>%
#   mutate(SUBJECT_ORIGINAL = gsub("HB0", "HB", SUBJECT_ORIGINAL)) %>%
#   arrange(SUBJECT_ORIGINAL) %>%
#   mutate(Cohort_code = gsub(pattern = "Unclassified", replacement = "OUT", Cohort_code)) %>%
#   mutate(Cohort_code_2 = gsub(pattern = "Unclassified", replacement = "OUT", Cohort_code_2)) %>%
#   rename("GROUP_Q1to6" = "Cohort_code") %>%
#   rename("GROUP_Q7" = "Cohort_code_2")



#create phenotype file to use in FEMPipeline
phenotype_info <- metadata %>%
  rename("Sample" = "SUBJECT_ORIGINAL") %>%
  mutate(Biomarker = "sncRNA", .after = "Sample") %>%  
  mutate(Technology = "RNASeq", .after = "Biomarker") 


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
            file = "Data/transcriptomic_phenotype.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)


# test if read using FEM pipeline works

# phenotype_file <- "Data/transcriptomic_phenotype.txt"
# classification_criteria <- "PREOPEVsMET"
# filter <- expression(TRUE)
# phenotype <- read.table(phenotype_file, header=TRUE, sep="\t")
# all_equal(phenotype, phenotype_info)
# 
# data <- read.table("Data/RNA/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
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


##############################


#metadata for new cohort - glionet cohort

metadata_glionet <- read_xlsx("Data/RNA_validation/External_validation_cohort_glionet.xlsx")[, 1:7]
colnames(metadata_glionet) <- c("patient_id", "sample_id", "age", "gender", 
                                "sample_instance", "sample_category", "test_label")

#filling empty values of age and gender
#assuming the structure :age, gender present in one row and then not present in next few rows
for(i in c(1:dim(metadata_glionet)[1])){
  if(!is.na(metadata_glionet[i, "age"]) && !is.na(metadata_glionet[i, "gender"])){
    specified_age <- round(metadata_glionet[i, ]$age, 2)
    specified_gender <- metadata_glionet[i, ]$gender
    print(specified_age)
    print(specified_gender)
  }
  metadata_glionet[i, ]$age <- specified_age
  metadata_glionet[i, ]$gender <- specified_gender
}

metadata_glionet[52:55, ]$age <- round(metadata_glionet[52:55, ]$age, 1) 
metadata_glionet <- metadata_glionet %>%
  mutate(sample_id = gsub("-", "", sample_id, fixed = TRUE))


metadata_glionet <- metadata_glionet %>%
  mutate(category_old_name = case_when(is.na(sample_category) ~ "UNK",
                                       sample_category == "PRE-OP" ~ "PREOPE",
                                       sample_category == "POST-OP" ~ "POSTOPE_TP",
                                       sample_category == "RECURRENCE" ~ "REC_TP",
                                       TRUE ~ NA_character_),
         .after = sample_category)

write.csv(metadata_glionet, "Data/RNA_validation/metadata_glionet.csv", row.names = FALSE)


summary(factor(metadata_glionet$sample_category))


#creating phenotype for validation cohort and combined cohort
phenotype <- read.table("Data/transcriptomic_phenotype.txt", header=TRUE, sep="\t")
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
  )

write.table(validation_phenotype, 
            file = "Data/transcriptomic_phenotype_validation.txt", 
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
            file = "Data/transcriptomic_phenotype_combined.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)

# ################################################################
# # test if read using FEM pipeline works validation cohort
# 
# phenotype_file <- "Data/transcriptomic_phenotype_validation.txt"
# classification_criteria <- "POSTOPE_TPVsREC_TP"
# filter <- expression(TRUE)
# phenotype <- read.table(phenotype_file, header=TRUE, sep="\t")
# 
# data <- read.table("Data/RNA/umi_counts_validation_cohort_common_tr.csv", header=TRUE, sep=",",
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
# ################################################################
# # test if read using FEM pipeline works validation cohort post combat
# 
# phenotype_file <- "Data/transcriptomic_phenotype_validation.txt"
# classification_criteria <- "POSTOPE_TPVsREC_TP"
# filter <- expression(TRUE)
# phenotype <- read.table(phenotype_file, header=TRUE, sep="\t")
# 
# data <- read.table("Data/RNA/validation_data.combat.POSTOPE_TPVsREC_TP.csv", header=TRUE, sep=",",
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
# ################################################################
# # test if read using FEM pipeline works combined cohort
# 
# phenotype_file <- "Data/transcriptomic_phenotype_combined.txt"
# classification_criteria <- "POSTOPE_TPVsREC_TP"
# filter <- expression(TRUE)
# phenotype <- read.table(phenotype_file, header=TRUE, sep="\t")
# 
# data <- read.table("Data/RNA/combined_data.combat.POSTOPE_TPVsREC_TP.csv", header=TRUE, sep=",",
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
