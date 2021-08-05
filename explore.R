if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MSstats")

library(MSstats)
library(tidyverse)
# library(stringr)   for str_detect()


setwd("/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV")
source("sample_column_generator.R")

raw_data <- read.csv("MSstatsinputannotatedQ1-6/MSstats_Input_Batch1.csv")
raw_data <- SkylinetoMSstatsFormat(raw_data)


raw_data_filtered <- raw_data %>%
  filter(BioReplicate == 'HB01')

data_process_output <- dataProcess(raw_data, logTrans = "2", normalization = "equalizeMedians")

dataProcessPlots(data_process_output, type = "ProfilePlot")
# dataProcessPlots(data_process_output, type = "QCPlot")
# dataProcessPlots(data_process_output, type = "ConditionPlot")

levels(data_process_output$ProcessedData$GROUP_ORIGINAL)

processed_data <- data_process_output$ProcessedData 
run_level_data <- data_process_output$RunlevelData

normed <- run_level_data %>% 
  select(Protein, LogIntensities, GROUP_ORIGINAL, SUBJECT_ORIGINAL) 

length(levels(normed$Protein))



prots <- strsplit(as.character(normed$Protein), split="\\|")
prots <- unlist(lapply(prots, function(x) x[3]))
prots <- gsub("_HUMAN", "", prots)



normed$Protein <- prots
normed <- normed %>% 
  spread(Protein, LogIntensities)


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

normed_abundance[is.na(normed_abundance)] <- 0

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