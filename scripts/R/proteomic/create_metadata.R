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
