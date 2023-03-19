library(tidyverse)
library(readxl)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

###############################################

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

write.csv(umi_counts, "Data/RNA/umi_counts.csv")



###############
#validation rna data glionet

data <- read_excel("Data/RNA_validation/218924.all_samples.summary.xlsx", 
                   sheet = "miRNA_piRNA")
mirna_data <- data[1:226,]
pirna_data <- data[228:395,] %>%
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

umi_counts1 <- umi_counts

write.csv(umi_counts, "Data/RNA_validation/umi_counts.csv")





data <- read_excel("Data/qiagen_results_test/gbmplasmev_validation_cohort_rna/218950.all_samples.summary.xlsx", 
                   sheet = "miRNA_piRNA")
mirna_data <- data[1:226,]
pirna_data <- data[228:395,] %>%
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
colnames(umi_counts) <- paste("S", colnames(umi_counts), sep = "")

samples <- paste("SB", c(1:55), sep = "")

umi_counts <- umi_counts[, samples]

umi_counts1 <- umi_counts1[, samples]

all.equal(umi_counts, umi_counts1)


##################


#initial cohort again

data <- read_excel("Data/qiagen_results_test/GBMPlasmaEV_initial_cohort_RNAData/219155.all_samples.summary.xlsx", 
                   sheet = "miRNA_piRNA")
mirna_data <- data[1:324,]
pirna_data <- data[326:461,] %>%
  separate(miRNA, c("miRNA", NA, NA, NA), sep = "/")

dim(data)[1]
dim(mirna_data)[1] + dim(pirna_data)[1]

data <- rbind(mirna_data, pirna_data) %>%
  column_to_rownames("miRNA")

umi_counts_a <- data %>%
  select(ends_with("UMIs"))
colnames(umi_counts_a) <- gsub("-UMIs", "", colnames(umi_counts_a))
colnames(umi_counts_a) <- sapply(colnames(umi_counts_a), FUN = 
                                 function(x){
                                   strsplit(x, split = "_", fixed = TRUE)[[1]][1]
                                 }
)


#initial cohort again - one lane

#initial cohort again

data <- read_excel("Data/qiagen_results_test/GBMPlasmaEV_initial_cohort_RNAData_onelane/219158.all_samples.summary.xlsx", 
                   sheet = "miRNA_piRNA")
mirna_data <- data[1:324,]
pirna_data <- data[326:461,] %>%
  separate(miRNA, c("miRNA", NA, NA, NA), sep = "/")

dim(data)[1]
dim(mirna_data)[1] + dim(pirna_data)[1]

data <- rbind(mirna_data, pirna_data) %>%
  column_to_rownames("miRNA")

umi_counts_b <- data %>%
  select(ends_with("UMIs"))
colnames(umi_counts_b) <- gsub("-UMIs", "", colnames(umi_counts_b))
colnames(umi_counts_b) <- sapply(colnames(umi_counts_b), FUN = 
                                 function(x){
                                   strsplit(x, split = "_", fixed = TRUE)[[1]][1]
                                 }
)


samples <- data.frame(samples = colnames(umi_counts_a)) %>%
  arrange(samples)

umi_counts_a <- umi_counts_a[, samples$samples]
umi_counts_b <- umi_counts_b[, samples$samples]

umi_counts <- umi_counts[, samples$samples]

all.equal(umi_counts_a, umi_counts_b)

all.equal(umi_counts_a, umi_counts)



#######################################

####### initial cohort and validation cohort new geneglobe execution #######


#initial cohort

data <- read_excel("Data/qiagen_results_test/GBMPlasmaEV_initial_cohort_RNAData/219155.all_samples.summary.xlsx", sheet = "miRNA_piRNA")
mirna_data <- data[1:324,]
pirna_data <- data[326:461,] %>%
  separate(miRNA, c("miRNA", NA, NA, NA), sep = "/")

dim(data)[1]
dim(mirna_data)[1] + dim(pirna_data)[1]

data <- rbind(mirna_data, pirna_data) %>%
  mutate(miRNA = gsub("-", "_", miRNA, fixed = TRUE)) %>%
  column_to_rownames("miRNA")

umi_counts <- data %>%
  select(ends_with("UMIs"))
colnames(umi_counts) <- gsub("-UMIs", "", colnames(umi_counts))
colnames(umi_counts) <- sapply(colnames(umi_counts), FUN = 
                                 function(x){
                                   strsplit(x, split = "_", fixed = TRUE)[[1]][1]
                                 }
)

write.csv(umi_counts, "Data/RNA/umi_counts_initial_cohort.csv")

# test read

# phenotype_file_name = "Data/transcriptomic_phenotype.txt"
# read_count_dir_path = "Data/RNA"
# read_count_file_name = "umi_counts_initial_cohort.csv"
# sep = ","
# dataset_id = "GBM_tr_initial"
# classification_criteria = "PREOPEVsPOSTOPE_TP"
# classes = c("POSTOPE_TP", "PREOPE")
# cores = 16
# results_dir_path = "fem_pipeline_results_tr"
# norm = "norm_log_cpm_simple"
# filter_expression = expression(TRUE)
# 
# read_count_file_path <- paste(read_count_dir_path, read_count_file_name, sep = "/")
# 
# data <- read.table(read_count_file_path, header=TRUE, sep=sep, row.names=1, skip=0,
#                    nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
# #data format : (transcripts x samples)
# data[is.na(data)] <- 0
# phenotype <- read.table(phenotype_file_name, header=TRUE, sep="\t")
# 
# extracted_samples <- phenotype %>% subset(eval(filter_expression))
# extracted_samples <- extracted_samples[!is.na(extracted_samples[classification_criteria]), ]
# extracted_samples$Sample <- factor(extracted_samples$Sample)
# 
# filtered_samples_read_count <- data %>% dplyr::select(extracted_samples$Sample)
# 
# #from the extracted_samples, select the 'Sample' column and classification_criteria column
# filtered_samples_output_labels <- extracted_samples[, c('Sample', classification_criteria)]
# colnames(filtered_samples_output_labels) <- c("Sample", "Label")




# validation cohort (i.e. glionet cohort)

data <- read_excel("Data/qiagen_results_test/gbmplasmev_validation_cohort_rna/218950.all_samples.summary.xlsx", 
                   sheet = "miRNA_piRNA")
mirna_data <- data[1:226,]
pirna_data <- data[228:395,] %>%
  separate(miRNA, c("miRNA", NA, NA, NA), sep = "/")

dim(data)[1]
dim(mirna_data)[1] + dim(pirna_data)[1]

data <- rbind(mirna_data, pirna_data) %>%
  mutate(miRNA = gsub("-", "_", miRNA, fixed = TRUE)) %>%
  column_to_rownames("miRNA")

umi_counts <- data %>%
  select(ends_with("UMIs"))
colnames(umi_counts) <- gsub("-UMIs", "", colnames(umi_counts))
colnames(umi_counts) <- sapply(colnames(umi_counts), FUN = 
                                 function(x){
                                   strsplit(x, split = "_", fixed = TRUE)[[1]][1]
                                 }
)

write.csv(umi_counts, "Data/RNA/umi_counts_validation_cohort.csv")



######## verifying if data output files for validation cohort obtained via 
###################### 2 different Geneglobe quantification jobs are the same

validation_data1 <- read.csv("Data/RNA_validation/umi_counts.csv", row.names = 1)
validation_data2 <- read.csv("Data/RNA/umi_counts_validation_cohort.csv", row.names = 1)

colnames(validation_data2) <- paste0("S", colnames(validation_data2))

sample_names1 <- sort(colnames(validation_data1))
sample_names2 <- sort(colnames(validation_data2))

all.equal(sample_names1, sample_names2)

validation_data1 <- validation_data1[, sample_names1]
rownames(validation_data1) <- gsub("-", "_", rownames(validation_data1), fixed = TRUE)

validation_data2 <- validation_data2[, sample_names2]

all.equal(validation_data1, validation_data2)

####################################################################################



#fetching common set of transcripts and write to file

data_file_path <- "Data/RNA/umi_counts_initial_cohort.csv"
validation_data_file_path <- "Data/RNA/umi_counts_validation_cohort.csv"      

data <- read.csv(data_file_path, row.names = 1)
validation_data <- read.csv(validation_data_file_path, row.names = 1)
colnames(validation_data) <- paste0("S", colnames(validation_data))

common <- intersect(rownames(data), rownames(validation_data))  
data.common <- data[common, ]
validation_data.common <- validation_data[common, ]

write.csv(data.common, "Data/RNA/umi_counts_initial_cohort_common_tr.csv")
write.csv(validation_data.common, "Data/RNA/umi_counts_validation_cohort_common_tr.csv")
# data.common <- read.csv("Data/RNA/umi_counts_initial_cohort_common_tr.csv",
#                         row.names = 1)
