library(tidyverse)
library(readxl)

library(ComplexHeatmap)
library(viridis)
library(RColorBrewer)

library(checkmate)

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



###############################################
#176 samples quantified result from RNA Seq Portal

data <- read_excel("Data/RNA_all/sRNAmatrix.xlsx")
sum(is.na(data$Name))

which(is.na(data$Name))

data[2856, ]

data <- data %>%
  filter(!is.na(Name)) %>%
  mutate(across(!contains("Name"), as.numeric))
dim(data)
# [1] 6203  177

data_grouped <- data %>% 
  group_by(Name) %>%
  summarize(n = n()) %>%
  filter(n > 1)
#0

length(unique(data$Name))
#[1] 6203


data <- data %>%
  column_to_rownames("Name")

min(data)
#0
#no -ve numbers

rsum <- rowSums(data)

data <- data %>%
  filter(rsum != 0)
dim(data)
# [1] 4804  176


meta_data <- read.csv("Data/transcriptomic_metadata_2023_176samples.csv") %>%
  filter(Condition != "OUT")

data <- data[, meta_data$Sample]
colnames(data) <- meta_data$orig_sample_name
dim(data)
# [1] 4804  165

rsum <- rowSums(data)

data <- data %>%
  filter(rsum != 0)
dim(data)
# [1] 4760  165

write.csv(data, "Data/RNA_all/umi_counts.csv")


#perform some qc to check :
#for each of the transcripts, perc of zeroes of that transcript across all samples
zero_entries <- as.data.frame(data == 0)
perc_zero <- data.frame(perc = 100 * rowSums(zero_entries) / ncol(zero_entries))

summary(perc_zero$perc)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   95.15   98.79   88.99   99.39   99.39


ggplot(perc_zero, aes(x = "", y = perc)) +
  geom_boxplot() +
  labs(x = "",
       y = "% of zeroes across all samples",
       title = "Boxplot of % of zeroes in each transcript across all samples")
ggsave("plots_RNA_all/qc/perc_zero_boxplot.png")


ggplot(perc_zero, aes(x = perc)) +
  geom_histogram(binwidth = 1) +
  labs(x = "% of zeroes across all samples",
       y = "number of transcripts",
       title = "Frequency of % of zeroes in each transcript across all samples",
       caption = paste("Total number of transcripts =", nrow(perc_zero))) +
  scale_x_continuous(breaks = seq(0, 100, by = 5)) +
  scale_y_continuous(n.breaks = 10)
ggsave("plots_RNA_all/qc/perc_zero_barplot.png")


ggplot(perc_zero, aes(x = perc)) +
  stat_ecdf(geom = "step") +
  labs(x = "% of zeroes across all samples",
       y = "Cumulative proportion of transcripts",
       title = "Cumulative distribution of % of zeroes in each transcript across all samples",
       caption = paste("Total number of transcripts =", nrow(perc_zero))) +
  scale_x_continuous(breaks = seq(0, 100, by = 5)) +
  scale_y_continuous(n.breaks = 10)
ggsave("plots_RNA_all/qc/perc_zero_cumulative.png")


#from the cumulative dist plot, about 0.15 transcripts have < 85% zeros across all samples
#verifying
sum(perc_zero$perc < 85) / nrow(perc_zero) 
# [1] 0.1569328

sum(perc_zero$perc < 90) / nrow(perc_zero) 
# [1] 0.1817227

sum(perc_zero$perc < 95) / nrow(perc_zero) 
# [1] 0.2378151

#median
sum(perc_zero$perc < 98.79) / nrow(perc_zero)
# [1] 0.555042

#3rd quartile
sum(perc_zero$perc < 99.39) / nrow(perc_zero)
# [1] 0.555042

sum(perc_zero$perc < 99.394) / nrow(perc_zero)

filt_data <- data[perc_zero$perc < 85, ]
nrow(filt_data)
#747

filt_data <- data[perc_zero$perc < 90, ]
nrow(filt_data)
#865

filt_data <- data[perc_zero$perc < 95, ]
nrow(filt_data)
#1132

#using median
filt_data <- data[perc_zero$perc < 98.79, ]
nrow(filt_data)
#2642

#using 3rd quartile
filt_data <- data[perc_zero$perc < 99.39, ]
nrow(filt_data)
#2642


#use which of the above ?
# maybe compare using expression heatmap

meta_data <- meta_data %>%
  dplyr::select(-c(Sample)) %>%
  mutate(Condition = factor(Condition),
         Cohort = factor(Cohort),
         Subgroup = factor(Subgroup))

# expr_data <- data
create_expression_heatmap <- function(expr_data, meta_data, file_name, 
                                      main_title = "", plot_dir_path = "plots_RNA_all/qc/heatmap/"){
  expr_data <- expr_data[, meta_data$orig_sample_name]
  all.equal(meta_data$orig_sample_name, colnames(expr_data))
  
  data_to_plot <- as.matrix(log2(expr_data + 2^-10))
  
  annotation_col_list <- list("Cohort" = c("cohort1" = "coral",
                                           "cohort2" = "cyan"))
  annotation_col_list[["Condition"]] <- brewer.pal(n = length(levels(meta_data$Condition)), name = "Dark2") 
  names(annotation_col_list[["Condition"]]) <- levels(meta_data$Condition)
  
  annotation_col_list[["Subgroup"]] <- brewer.pal(n = length(levels(meta_data$Subgroup)), name = "Paired") 
  names(annotation_col_list[["Subgroup"]]) <- levels(meta_data$Subgroup)
  
  ht <- Heatmap(data_to_plot, name = "Log2\nTranscriptomics\nexpression",
                col = viridis(n = 10, option = "magma"),
                show_column_names = FALSE,
                show_column_dend = FALSE,
                show_row_names = FALSE,
                show_row_dend = FALSE,
                row_title = paste0("Transcripts (", nrow(data_to_plot), ")"),
                column_title = paste0("Samples (", ncol(data_to_plot), ")"),
                bottom_annotation = HeatmapAnnotation(
                  "Condition" = meta_data$Condition,
                  "Subgroup" = meta_data$Subgroup,
                  "Cohort" = meta_data$Cohort,
                  col = annotation_col_list
                ))
  if(!dir.exists(plot_dir_path)){
    dir.create(plot_dir_path, recursive = TRUE)
  }
  file_path <- paste0(plot_dir_path, file_name)
  png(file_path, units = "cm", width = 20, height = 25, res = 1200)  
  draw(ht, column_title = main_title, column_title_gp = gpar(fontsize = 15))    
  dev.off()
}


create_expression_heatmap(data, meta_data, "1_zero_filtered.png")

filt_data <- data[perc_zero$perc < 98.79, ]
nrow(filt_data)
create_expression_heatmap(filt_data, meta_data, "2_median_filtered.png")

filt_data <- data[perc_zero$perc < 95, ]
create_expression_heatmap(filt_data, meta_data, "3_95_filtered.png")

filt_data <- data[perc_zero$perc < 90, ]
create_expression_heatmap(filt_data, meta_data, "4_90_filtered.png")

filt_data <- data[perc_zero$perc < 85, ]
create_expression_heatmap(filt_data, meta_data, "5_85_filtered.png")



###################
#filtering just PREOPE, MET, HC samples

data <- read_excel("Data/RNA_all/sRNAmatrix.xlsx")
sum(is.na(data$Name))

which(is.na(data$Name))

data[2856, ]

data <- data %>%
  filter(!is.na(Name)) %>%
  mutate(across(!contains("Name"), as.numeric)) %>%
  column_to_rownames("Name")
dim(data)
# [1] 6203  176

min(data)
#0
#no -ve numbers

meta_data <- read.csv("Data/transcriptomic_metadata_2023_176samples.csv") %>%
  filter(Condition %in% c('PREOPE', 'MET', 'HC'))

data <- data[, meta_data$Sample]
colnames(data) <- meta_data$orig_sample_name
dim(data)
#[1] 6203   67

rsum <- rowSums(data)

data <- data %>%
  filter(rsum != 0)
dim(data)
#[1] 3087   67

write.csv(data, "Data/RNA_all/newquant_Nov2023_umi_counts_PREOPE_MET_HC.csv")


#perform some qc to check :
#for each of the transcripts, perc of zeroes of that transcript across all samples
zero_entries <- as.data.frame(data == 0)
perc_zero <- data.frame(perc = 100 * rowSums(zero_entries) / ncol(zero_entries))

summary(perc_zero$perc)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   86.57   97.01   82.37   98.51   98.51


ggplot(perc_zero, aes(x = "", y = perc)) +
  geom_boxplot() +
  labs(x = "",
       y = "% of zeroes across all samples",
       title = "Boxplot of % of zeroes in each transcript across all samples")
ggsave("plots_RNA_all/PREOPE_MET_HC/qc/perc_zero_boxplot.png")


ggplot(perc_zero, aes(x = perc)) +
  geom_histogram(binwidth = 1) +
  labs(x = "% of zeroes across all samples",
       y = "number of transcripts",
       title = "Frequency of % of zeroes in each transcript across all samples",
       caption = paste("Total number of transcripts =", nrow(perc_zero))) +
  scale_x_continuous(breaks = seq(0, 100, by = 5)) +
  scale_y_continuous(n.breaks = 10)
ggsave("plots_RNA_all/PREOPE_MET_HC/qc/perc_zero_barplot.png")


ggplot(perc_zero, aes(x = perc)) +
  stat_ecdf(geom = "step") +
  labs(x = "% of zeroes across all samples",
       y = "Cumulative proportion of transcripts",
       title = "Cumulative distribution of % of zeroes in each transcript across all samples",
       caption = paste("Total number of transcripts =", nrow(perc_zero))) +
  scale_x_continuous(breaks = seq(0, 100, by = 5)) +
  scale_y_continuous(n.breaks = 10)
ggsave("plots_RNA_all/PREOPE_MET_HC/qc/perc_zero_cumulative.png")


#from the cumulative dist plot, about 0.3 transcripts have < 90% zeros across all samples
#verifying
sum(perc_zero$perc < 90) / nrow(perc_zero) 
# [1] 0.2938128

sum(perc_zero$perc < 95) / nrow(perc_zero) 
# [1] 0.3751215

#median
sum(perc_zero$perc < 97.01) / nrow(perc_zero)
# [1] 0.4331066

#3rd quartile
sum(perc_zero$perc < 98.51) / nrow(perc_zero)
# [1] 1

sum(perc_zero$perc < 98.5074) / nrow(perc_zero)
# [1] 0.5484289

sum(perc_zero$perc == 0) / nrow(perc_zero)
# [1] 0.04275996


filt_data <- data[perc_zero$perc == 0, ]
nrow(filt_data)
# [1] 132

filt_data <- data[perc_zero$perc < 90, ]
nrow(filt_data)
#907

filt_data <- data[perc_zero$perc < 95, ]
nrow(filt_data)
#1158

#using median
filt_data <- data[perc_zero$perc < 97.01, ]
nrow(filt_data)
#1337

#using value just before 3rd quartile where cumulative prop is < 1
filt_data <- data[perc_zero$perc < 98.5074, ]
nrow(filt_data)
#1693

#using 3rd quartile
filt_data <- data[perc_zero$perc < 98.51, ]
nrow(filt_data)
#3087


#use which of the above ?
# maybe compare using expression heatmap

meta_data <- meta_data %>%
  dplyr::select(-c(Sample)) %>%
  mutate(Condition = factor(Condition),
         Cohort = factor(Cohort),
         Subgroup = factor(Subgroup))




create_expression_heatmap(data, meta_data, "1_zero_filtered.png", 
                          main_title = "Transcripts non-zero in atleast 1 sample", 
                          plot_dir_path = "plots_RNA_all/PREOPE_MET_HC/qc/heatmap/")

filt_data <- data[perc_zero$perc < 98.5074, ]
nrow(filt_data)
create_expression_heatmap(filt_data, meta_data, "2_filtered.png", 
                          main_title = "Transcripts non-zero in less than 98.5074% samples", 
                          plot_dir_path = "plots_RNA_all/PREOPE_MET_HC/qc/heatmap/")


filt_data <- data[perc_zero$perc < 97.01, ]
nrow(filt_data)
create_expression_heatmap(filt_data, meta_data, "3_filtered_median.png", 
                          main_title = "Transcripts non-zero in less than\n97.01% samples (median 0%)", 
                          plot_dir_path = "plots_RNA_all/PREOPE_MET_HC/qc/heatmap/")


filt_data <- data[perc_zero$perc < 95, ]
nrow(filt_data)

non_zero <- filt_data != 0
sum(non_zero) / (nrow(non_zero) * ncol(non_zero))
# [1] 0.4358389

create_expression_heatmap(filt_data, meta_data, "4_filtered_95perc.png", 
                          main_title = "Transcripts non-zero in less than\n95% samples", 
                          plot_dir_path = "plots_RNA_all/PREOPE_MET_HC/qc/heatmap/")


filt_data <- data[perc_zero$perc < 90, ]
nrow(filt_data)

non_zero <- filt_data != 0
sum(non_zero) / (nrow(non_zero) * ncol(non_zero))
# [1] 0.5367704

create_expression_heatmap(filt_data, meta_data, "5_filtered_90perc.png", 
                          main_title = "Transcripts non-zero in less than\n90% samples", 
                          plot_dir_path = "plots_RNA_all/PREOPE_MET_HC/qc/heatmap/")


filt_data <- data[perc_zero$perc < 85, ]
nrow(filt_data)

non_zero <- filt_data != 0
sum(non_zero) / (nrow(non_zero) * ncol(non_zero))
# [1] 0.6327669

create_expression_heatmap(filt_data, meta_data, "6_filtered_85perc.png", 
                          main_title = "Transcripts non-zero in less than\n85% samples", 
                          plot_dir_path = "plots_RNA_all/PREOPE_MET_HC/qc/heatmap/")


filt_data <- data[perc_zero$perc < 80, ]
nrow(filt_data)

non_zero <- filt_data != 0
sum(non_zero) / (nrow(non_zero) * ncol(non_zero))
# [1] 0.6768023

create_expression_heatmap(filt_data, meta_data, "7_filtered_80perc.png", 
                          main_title = "Transcripts non-zero in less than\n80% samples", 
                          plot_dir_path = "plots_RNA_all/PREOPE_MET_HC/qc/heatmap/")


filt_data <- data[perc_zero$perc < 75, ]
nrow(filt_data)

non_zero <- filt_data != 0
sum(non_zero) / (nrow(non_zero) * ncol(non_zero))
# [1] 0.7180197

create_expression_heatmap(filt_data, meta_data, "8_filtered_75perc.png", 
                          main_title = "Transcripts non-zero in less than\n75% samples", 
                          plot_dir_path = "plots_RNA_all/PREOPE_MET_HC/qc/heatmap/")


filt_data <- data[perc_zero$perc ==0, ]
nrow(filt_data)
create_expression_heatmap(filt_data, meta_data, "9_filtered_0perc.png", 
                          main_title = "Transcripts non-zero in all samples", 
                          plot_dir_path = "plots_RNA_all/PREOPE_MET_HC/qc/heatmap/")


#going with <90, since thats the first to hit > 50% total non-zeros, 
#and also about 0.3 of transcripts have <90 across all samples

filt_data <- data[perc_zero$perc < 90, ]
nrow(filt_data)

write.csv(filt_data, "Data/RNA_all/newquant_Nov2023_umi_counts_PREOPE_MET_HC_filter90.csv")


######################################################

#create data subsets with DE significant features

data <- read.csv("Data/RNA_all/newquant_Nov2023_umi_counts_PREOPE_MET_HC_filter90.csv", row.names = 1)
metadata <- read.table("Data/transcriptomic_phenotype_PREOPE_MET_HC.txt", header=TRUE, sep="\t") %>%
  filter(!is.na(data_cohort)) %>%
  mutate(condition = case_when(Condition == "Pre-op" ~ "PREOPE",
                               Condition == "Metastatic" ~ "MET",
                               Condition == "Healthy Control" ~ "HC"))
data <- data[, metadata$Sample]

de_file_path <- "plots_RNA_all/PREOPE_MET_HC/volcano/volcano_1_PREOPEVsMET.csv"
comparison <- "PREOPEVsMET"
subset_file_path <- "Data/RNA_all/newquant_Nov2023_umi_counts_PREOPE_MET_HC_filter90_de_PREOPEVsMET.csv"
create_de_subset_file <- function(data, metadata, comparison, de_file_path,
                                  subset_file_path){
  de_data <- read.csv(de_file_path) %>%
    filter(significance != "Not significant") %>%
    dplyr::select(Molecule)
  
  metadata_sub <- metadata %>%
    dplyr::rename("comparison" = comparison) %>%
    dplyr::filter(!is.na(comparison))
  
  data_sub <- data[de_data$Molecule, metadata_sub$Sample]
  
  #https://github.com/r-lib/rlang/issues/1111
  checkmate::assert(all.equal(rownames(data_sub), de_data$Molecule))
  checkmate::assert(all.equal(colnames(data_sub), metadata_sub$Sample))  
  print(dim(data_sub))
  write.csv(data_sub, subset_file_path)
}

create_de_subset_file(data, metadata, 
                      comparison = "PREOPEVsMET", 
                      de_file_path = "plots_RNA_all/PREOPE_MET_HC/volcano/volcano_1_PREOPEVsMET.csv",
                      subset_file_path = "Data/RNA_all/newquant_Nov2023_umi_counts_PREOPE_MET_HC_filter90_de_PREOPEVsMET.csv")

create_de_subset_file(data, metadata, 
                      comparison = "PREOPEVsHC", 
                      de_file_path = "plots_RNA_all/PREOPE_MET_HC/volcano/volcano_2_PREOPEVsHC.csv",
                      subset_file_path = "Data/RNA_all/newquant_Nov2023_umi_counts_PREOPE_MET_HC_filter90_de_PREOPEVsHC.csv")

create_de_subset_file(data, metadata, 
                      comparison = "METVsHC", 
                      de_file_path = "plots_RNA_all/PREOPE_MET_HC/volcano/volcano_3_METVsHC.csv",
                      subset_file_path = "Data/RNA_all/newquant_Nov2023_umi_counts_PREOPE_MET_HC_filter90_de_METVsHC.csv")

