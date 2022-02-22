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

