#to check if metedata from protein annotation, cohorts for comparitive analyses, replicates file are same

library(tidyverse)
library(readxl)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

replicatesQ1to6 <- read_excel(path = "Data/Protein/ReplicatesQ1-6.xlsx")

replicatesQ7 <- read_excel(path = "Data/Protein/ReplicatesQ7.xlsx")

replicatesQ1to6 <- replicatesQ1to6 %>%
  select(Replicate, Condition) %>%
  rename("sample" = "Replicate") %>%
  rename("Q1to6condition" = "Condition")

replicatesQ7 <- replicatesQ7 %>%
  select(Replicate, Condition) %>%
  rename("sample" = "Replicate") %>%
  rename("Q7condition" = "Condition")

prot_meta_data <- replicatesQ1to6 %>%
  inner_join(replicatesQ7)

prot_meta_data <- prot_meta_data %>%
  separate(sample, c("modified_sample", NA), remove = FALSE, sep = "-", fill = "right") %>%
  mutate(str_part = gsub("[0-9]", "", modified_sample)) %>%
  mutate(num_part = strtoi(gsub("[^0-9]", "", modified_sample), base = 10)) %>%
  arrange(str_part, num_part) %>%
  select(-c(str_part, num_part, modified_sample))

prot_meta_data <- data.frame(prot_meta_data)

##########################################################

rna_seq_categories <- read_excel(path = "Data/Cohorts for comparative analyses.xlsx",
                                 sheet = "Categories_RNAseq")
rna_seq_categories <- rna_seq_categories %>%
  rename("sample" = "...1")

rna_seq_categories <- rna_seq_categories %>%
  mutate(str_part = gsub("[0-9]", "", sample)) %>%
  mutate(num_part = strtoi(gsub("[^0-9]", "", sample), base = 10)) %>%
  arrange(str_part, num_part) %>%
  select(-c(str_part, num_part))

rna_seq_categories <- data.frame(rna_seq_categories)

##################get metadata from proteomic data file labels############

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

metadata <- metadata %>%
  separate(SUBJECT_ORIGINAL, c("modified_sample", NA), remove = FALSE, sep = "-", fill = "right") %>%
  mutate(str_part = gsub("[0-9]", "", modified_sample)) %>%
  mutate(num_part = strtoi(gsub("[^0-9]", "", modified_sample), base = 10)) 
metadata <- metadata %>%
  arrange(str_part, num_part) %>%
  select(-c(str_part, num_part, modified_sample))
colnames(metadata) <- c("sample", "Q1to6condition", "Q7condition")


################check equality

str(prot_meta_data)
str(metadata)
str(rna_seq_categories)

all.equal(prot_meta_data, metadata)
#TRUE

rna_metadata <- rna_seq_categories %>%
  select(sample, Cohort_code, Cohort_code_2)

colnames(rna_metadata) <- c("sample", "Q1to6condition", "Q7condition")
rna_metadata <- rna_metadata %>%
  mutate_all(str_replace_all, "Unclassified", "OUT") %>%
  mutate_all(str_replace_all, "POSTOP-P", "POSTOPE-P") %>%
  mutate(sample = gsub("HB0", "HB", sample))

prot_meta_data2 <- prot_meta_data %>%
  mutate(sample = gsub("HB0", "HB", sample)) %>%
  filter(sample %in% rna_metadata$sample)

all.equal(rna_metadata[,1:2], prot_meta_data2[,1:2])
#TRUE
all.equal(rna_metadata, prot_meta_data2)
#"Component “Q7condition”: 38 string mismatches"



q7_rna_metadata <- rna_metadata %>%
  filter(Q7condition %in% c("PREOPE", "REC-TP"))
q7_prot_metadata <- prot_meta_data2 %>%
  filter(Q7condition %in% c("PREOPE", "REC-TP"))  

all.equal(q7_prot_metadata, q7_rna_metadata)
#TRUE



##################

write.csv(rna_metadata, file = "Data/to_test_rna_metadata.csv")
write.csv(prot_meta_data2, file = "Data/to_test_prot_metadata.csv")
