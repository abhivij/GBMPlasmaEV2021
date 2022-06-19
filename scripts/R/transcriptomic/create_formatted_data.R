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
