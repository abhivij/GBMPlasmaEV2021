if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MSstats")

library(MSstats)
library(tidyverse)
# library(stringr)   for str_detect()


setwd("/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV")
source("sample_column_generator.R")

data_path <- "Data/Protein/MSstatsinputannotatedQ1-6"
file_name <- list.files(path = data_path, full.names = TRUE)[1]

raw_data <- read.csv(file_name)
raw_data <- SkylinetoMSstatsFormat(raw_data)


raw_data_filtered <- raw_data %>%
  filter(BioReplicate == 'HB01')

data_process_output <- dataProcess(raw_data, logTrans = "2", normalization = "equalizeMedians")
data_process_output2 <- dataProcess(raw_data, logTrans = "2", normalization = "equalizeMedians",
                                    censoredInt = '0')
data_process_output_nonorm <- dataProcess(raw_data, logTrans = "2", normalization = FALSE)
data_process_output2_nonorm <- dataProcess(raw_data, logTrans = "2", normalization = FALSE,
                                           censoredInt = '0')
# dataProcessPlots(data_process_output, type = "ProfilePlot")
dataProcessPlots(data_process_output2, type = "QCPlot")
# dataProcessPlots(data_process_output, type = "ConditionPlot")


saveRDS(data_process_output2, file = "dataProcessOutput.rds")
data_process_output2_new <- readRDS("dataProcessOutput.rds")

levels(data_process_output$ProcessedData$GROUP_ORIGINAL)

processed_data <- data_process_output$ProcessedData 
run_level_data <- data_process_output$RunlevelData

run_level_data2 <- data_process_output_nonorm$RunlevelData

normed <- run_level_data %>% 
  select(Protein, LogIntensities, GROUP_ORIGINAL, SUBJECT_ORIGINAL)

length(levels(normed$Protein))


normed <- normed %>%
  separate(Protein, c(NA, "Protein", NA), sep = "\\|") %>% 
  pivot_wider(names_from = Protein, values_from = LogIntensities)


write.table(normed, "normQ1-6.csv", quote = FALSE, sep = ",", row.names = FALSE)


normed_abundance <- processed_data %>% 
  select(PROTEIN, FEATURE, ABUNDANCE, GROUP_ORIGINAL, SUBJECT_ORIGINAL) 
prots <- strsplit(as.character(normed_abundance$PROTEIN), split="\\|")
prots <- unlist(lapply(prots, function(x) x[3]))
prots <- gsub("_HUMAN", "", prots)
normed_abundance <- normed_abundance %>%
  mutate(PROTEIN = paste(prots, FEATURE, sep = "_")) %>%
  select(-FEATURE)
normed_abundance <- normed_abundance %>% 
  spread(PROTEIN, ABUNDANCE)

# not required since no NAs
# normed_abundance[is.na(normed_abundance)] <- 0

############PCA start###########################################

create_pca <- function(norm_data, filename, width = 30, height = 30){
  ne <- t(norm_data)
  
  # PCA on all samples and all proteins
  all_group_names <- ne["GROUP_ORIGINAL", ]
  all_sample_names <- ne["SUBJECT_ORIGINAL", ]
  colnames(ne) <- ne["SUBJECT_ORIGINAL", ]
  ne <- ne[-which(row.names(ne) %in% c("GROUP_ORIGINAL", "SUBJECT_ORIGINAL")), ]
  ne <- data.frame(ne)
  ne <- sapply(ne, function(x) as.numeric(x))  #causes rownames that is protein names to be lost
  
  p <- prcomp(t(ne), scale. = TRUE, center = TRUE)
  pca_plotdata <- data.frame(p$x) %>%
    mutate(type = all_group_names) %>%
    mutate(name = all_sample_names)
  
  pca_plotdata %>%
    ggplot(aes(x = PC1, y = PC2, colour = type, label = name)) +
    geom_point(size = 2) + 
    labs(title = "PCA for all proteins") +
    ggrepel::geom_text_repel(show.legend = FALSE)
  
  ggsave(filename, width = width, height = height, units = "cm")  
}

############PCA end###########################################

create_pca(norm_data = normed, filename = "PCA_intensity_all.jpg")

normed_q2 <- normed %>%
  filter(GROUP_ORIGINAL %in% c("PREOPE", "MET", "HC"))
create_pca(norm_data = normed_q2, filename = "PCA_intensity_Q2.jpg")

normed_qc <- normed %>%
  filter(GROUP_ORIGINAL %in% c("QC1", "QC2"))
create_pca(norm_data = normed_qc, filename = "PCA_intensity_QC.jpg", height = 10, width = 10)

normed_qc2 <- normed %>%
  filter(grepl("HB18|HC7", SUBJECT_ORIGINAL)) %>%
  inner_join(sample_column_qc, by = c("SUBJECT_ORIGINAL" = "Sample")) %>%
  select(-GROUP_ORIGINAL) %>%
  rename(GROUP_ORIGINAL = Column)
create_pca(norm_data = normed_qc2, filename = "PCA_intensity_QC2.jpg", height = 10, width = 10)

create_pca(norm_data = normed_abundance, filename = "PCA_abundance_all.jpg")

normed_abundance_q2 <- normed_abundance %>%
  filter(GROUP_ORIGINAL %in% c("PREOPE", "MET", "HC"))
create_pca(norm_data = normed_abundance_q2, filename = "PCA_abundance_Q2.jpg")

normed_abundance_qc <- normed_abundance %>%
  filter(GROUP_ORIGINAL %in% c("QC1", "QC2"))
create_pca(norm_data = normed_abundance_qc, filename = "PCA_abundance_QC.jpg", height = 10, width = 10)

normed_abundance_qc2 <- normed_abundance %>%
  filter(grepl("HB18|HC7", SUBJECT_ORIGINAL)) %>%
  inner_join(sample_column_qc, by = c("SUBJECT_ORIGINAL" = "Sample")) %>%
  select(-GROUP_ORIGINAL) %>%
  rename(GROUP_ORIGINAL = Column)
create_pca(norm_data = normed_abundance_qc2, filename = "PCA_abundance_QC2.jpg", height = 10, width = 10)
#############DE start################################

#aim is to do GBM Vs HC Vs MET

#for now doing HC Vs Preop-GBM Vs Preop-Met
# i.e. HC Vs PREOPE Vs MET
comparison <- rbind(c(-1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), # HC Vs PREOPE
                    c(-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), # HC Vs MET
                    c(0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0)) # PREOPE Vs MET


row.names(comparison) <- c("hc_preopgbm", 
                           "hc_met",
                           "preopgbm_met")
# make all the comparisons
multipleResults <- groupComparison(contrast.matrix = comparison,
                                   data=data_process_output)
############DE end##############################################



#reading sum of normalized areas

file_path <- "Data/Protein/Sumofnormalisedareas/Batch1-18_modified.csv"
sum_norm_area_data <- read.table(file_path, header = TRUE, sep = ",", 
                                 comment.char = "", na.strings = "#N/A",
                                 row.names = 1)
sum(is.na(sum_norm_area_data))
sum_norm_area_data[is.na(sum_norm_area_data)] <- 0.000001
sum_norm_area_data[is.na(sum_norm_area_data)] <- 0
sum(sum_norm_area_data == 0)
boxplot(sum_norm_area_data[1:10,])
boxplot(log2(sum_norm_area_data[1:10,]))

#read norm output data

norm_output1 <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_NA_equalizeMedians.csv")
norm_output2 <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_NA_FALSE.csv")
norm_output3 <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_column_equalizeMedians.csv")
norm_output5 <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_disease_equalizeMedians.csv")
norm_output8 <- read.csv(file = "Data/Protein/output/norm_annotatedQ7_NA_FALSE.csv")
norm_output16 <- read.csv(file = "Data/Protein/output/norm_unannotated_disease_FALSE.csv")

norm_output1 <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_NA_equalizeMedians.csv") %>%
  arrange(SUBJECT_ORIGINAL) %>%
  select(-c(GROUP_ORIGINAL))
norm_output2 <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_NA_FALSE.csv") %>%
  arrange(SUBJECT_ORIGINAL) %>%
  select(-c(GROUP_ORIGINAL))
norm_output5 <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_disease_equalizeMedians.csv") %>%
  arrange(SUBJECT_ORIGINAL) %>%
  select(-c(GROUP_ORIGINAL))
norm_output8 <- read.csv(file = "Data/Protein/output/norm_annotatedQ7_NA_FALSE.csv") %>%
  arrange(SUBJECT_ORIGINAL) %>%
  select(-c(GROUP_ORIGINAL))
norm_output16 <- read.csv(file = "Data/Protein/output/norm_unannotated_disease_FALSE.csv") %>%
  arrange(SUBJECT_ORIGINAL) %>%
  select(-c(GROUP_ORIGINAL))


modified_sum_norm_area <- log2(sum_norm_area_data)
modified_sum_norm_area <- data.frame(t(modified_sum_norm_area)) %>%
  rownames_to_column("SUBJECT_ORIGINAL") %>%
  arrange(SUBJECT_ORIGINAL)


all.equal(norm_output1, norm_output5)
all.equal(norm_output1, norm_output2)
all.equal(norm_output2, norm_output5)

all.equal(norm_output2, norm_output8)
all.equal(norm_output8, norm_output16)
all.equal(norm_output16, norm_output2)

all.equal(norm_output2, modified_sum_norm_area)
