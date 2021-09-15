library(tidyverse)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

source("scripts/R/utils.R")

data <- read_excel("Data/RNA/158629.all_samples.summary.xlsx", sheet = "miRNA_piRNA")
mirna_data <- data[1:2505,]
pirna_data <- data[2507:2642,] %>%
  separate(miRNA, c("miRNA", NA, NA, NA), sep = "/")

dim(data)[1]
dim(mirna_data)[1] + dim(pirna_data)[1]

data <- rbind(mirna_data, pirna_data) %>%
  column_to_rownames("miRNA")

umi_counts <- data %>%
  select(ends_with("UMIs"))
colnames(umi_counts) <- gsub("-UMIs", "", colnames(umi_counts))
colnames(umi_counts) <- sapply(colnames(umi_counts), FUN = 
                                 function(x){
                                   strsplit(x, split = "_", fixed = TRUE)[[1]][1]
                                 }
)





#to get meta data start

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
  left_join(group_mapping) %>%
  mutate(GROUP_ORIGINAL = gsub("-", "_", GROUP_ORIGINAL)) %>%
  rename(GROUP_Q1to6 = GROUP_ORIGINAL)



#get group mappings from Q7 proteomics output
protein_data <- read.csv(file = "Data/Protein/norm_output/norm_annotatedQ7_NA_equalizeMedians.csv")

group_mapping <- protein_data %>%
  select(SUBJECT_ORIGINAL, GROUP_ORIGINAL) %>%
  arrange(SUBJECT_ORIGINAL)
rm(protein_data)

group_mapping <- group_mapping %>%
  mutate(SUBJECT_ORIGINAL = gsub("HB0", "HB", SUBJECT_ORIGINAL))

metadata <- metadata %>%
  left_join(group_mapping) %>%
  mutate(GROUP_ORIGINAL = gsub("-", "_", GROUP_ORIGINAL)) %>%
  rename(GROUP_Q7 = GROUP_ORIGINAL)


write.csv(metadata, "Data/metadata.csv", row.names = FALSE)

