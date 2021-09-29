library(tidyverse)
library(readxl)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

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


#order metadata as HC1, HC2, HC3, ...
#   and NOT HC1, HC10, HC11 ... HC19, HC2, HC20, ...
# to get this ordering split string part and num part, and order separately
# for some unknown reason, 'num_part' entries corresponding to HB08, HB09 become NA. 
#                                 So setting them manually
metadata <- metadata %>%
  separate(SUBJECT_ORIGINAL, c("modified_sample", NA), remove = FALSE, sep = "-", fill = "right") %>%
  mutate(str_part = gsub("[0-9]", "", modified_sample)) %>%
  mutate(num_part = strtoi(gsub("[^0-9]", "", modified_sample))) 
metadata[metadata$SUBJECT_ORIGINAL == "HB08", "num_part"] <- 8
metadata[metadata$SUBJECT_ORIGINAL == "HB09", "num_part"] <- 9
metadata <- metadata %>%
  arrange(str_part, num_part) %>%
  select(-c(str_part, num_part, modified_sample))

write.csv(metadata, "Data/proteomic_sample_metadata.csv", row.names = FALSE)


#create phenotype file to use in FEMPipeline
phenotype_info <- metadata %>%
  rename("Sample" = "SUBJECT_ORIGINAL") %>%
  mutate(Biomarker = "Protein", .after = "Sample") %>%  
  mutate(Technology = "SWATH-MS", .after = "Biomarker") %>%
  mutate(PREOPEVsMET = ifelse(GROUP_Q1to6 == "PREOPE", 
                              "PREOPE", ifelse(GROUP_Q1to6 == "MET", "MET", NA))) %>%
  mutate(PREOPEVsHC = ifelse(GROUP_Q1to6 == "PREOPE", 
                             "PREOPE", ifelse(GROUP_Q1to6 == "HC", "HC", NA))) %>%
  mutate(METVsHC = ifelse(GROUP_Q1to6 == "MET", 
                          "MET", ifelse(GROUP_Q1to6 == "HC", "HC", NA)))

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

