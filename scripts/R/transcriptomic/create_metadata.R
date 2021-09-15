library(tidyverse)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

umi_counts <- read.csv("Data/RNA/umi_counts.csv", row.names = 1)

#get group mappings from Q1-6 proteomics output
protein_data <- read.csv(file = "Data/Protein/norm_output/norm_annotatedQ1-6_NA_equalizeMedians.csv")

group_mapping <- protein_data %>%
  select(SUBJECT_ORIGINAL, GROUP_ORIGINAL) %>%
  arrange(SUBJECT_ORIGINAL)
rm(protein_data)

metadata <- data.frame(SUBJECT_ORIGINAL = colnames(umi_counts)) %>%
  arrange(SUBJECT_ORIGINAL)

#in proteomics data, unlike transcriptomics, HB numbered 01, 02, 03, ...
# so changing that to match transcriptomics sample names
group_mapping <- group_mapping %>%
  mutate(SUBJECT_ORIGINAL = gsub("HB0", "HB", SUBJECT_ORIGINAL))

metadata <- metadata %>%
  inner_join(group_mapping) %>%
  mutate(GROUP_ORIGINAL = gsub("-", "_", GROUP_ORIGINAL)) %>%
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
  mutate(GROUP_ORIGINAL = gsub("-", "_", GROUP_ORIGINAL)) %>%
  rename(GROUP_Q7 = GROUP_ORIGINAL)

group_mapping[!group_mapping$SUBJECT_ORIGINAL %in% metadata$SUBJECT_ORIGINAL,]

write.csv(metadata, "Data/metadata.csv", row.names = FALSE)

