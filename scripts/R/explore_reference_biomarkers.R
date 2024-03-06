#analyze and compare proteins and miRNAs from previous studies

library(tidyverse)
library(readxl)
library(ggplot2)
library(ggvenn)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

proteins.prev.plasmaEVs <- read_excel("Data/selected_features/features_of_interest.xlsx", sheet = "PlasmaEV_GBMVsHC")
proteins.prev.urineEVs <- read_excel("Data/selected_features/features_of_interest.xlsx", sheet = "UrineEV_GBMVsHC")

proteins.our_study <- read.csv("Data/Protein/formatted_data/all_protein_names.csv")


data <- read.csv("Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv", 
                 row.names=1)
validation_data <- read.csv("Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil.csv", 
                            row.names = 1)

common <- intersect(rownames(data), rownames(validation_data))  
#4117

proteins.our_study <- proteins.our_study %>%
  dplyr::filter(from_id %in% common)

plot_dir_path <- "plots_reference_biomarkers/venn/"
if(!dir.exists(plot_dir_path)){
  dir.create(plot_dir_path, recursive = TRUE)
}

common.prev_study <- intersect(proteins.prev.plasmaEVs$Protein, proteins.prev.urineEVs$Protein)
ggvenn(list("Plasma EV study" = proteins.prev.plasmaEVs$Protein, 
            "Urine EV study" = proteins.prev.urineEVs$Protein),
       stroke_size = 0.1,
       set_name_size = 5,
       text_size = 3,
       fill_color = c("brown1", "navajowhite")) +
  ggtitle("GBM Vs HC significant proteins") +
  labs(caption = paste("Common :", paste(common.prev_study, collapse = ", "))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.5)),
        plot.caption = element_text(hjust = 0.5))
ggsave(paste0(plot_dir_path, "1_PlasmaVsUrine.jpg"))

ggvenn(list("Cohort 1" = rownames(data), 
            "Cohort 2" = rownames(validation_data)),
       stroke_size = 0.1,
       set_name_size = 5,
       text_size = 3,
       fill_color = c("cyan", "orchid")) +
  ggtitle("Proteins identified in our study") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.5)))
ggsave(paste0(plot_dir_path, "2_OurStudy_Cohort1VsCohort2.jpg"))

common.all <- Reduce(intersect, list(proteins.prev.plasmaEVs$Protein, 
                                     proteins.prev.urineEVs$Protein,
                                     proteins.our_study$primary_gene_id))
ggvenn(list("Plasma EV GBM Vs HC" = proteins.prev.plasmaEVs$Protein, 
            "Urine EV GBM Vs HC" = proteins.prev.urineEVs$Protein,
            "Proteins from our study" = proteins.our_study$primary_gene_id),
       stroke_size = 0.1,
       set_name_size = 5,
       text_size = 3,
       fill_color = c("brown1", "navajowhite", "steelblue1")) +
  ggtitle("GBM Vs HC significant proteins & our study proteins") +
  labs(caption = paste("Common :", paste(common.all, collapse = ", "))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.5)),
        plot.caption = element_text(hjust = 0.5))
ggsave(paste0(plot_dir_path, "3_PreviousStudiesVsOurStudy.jpg"))

'BROX' %in% (read.csv("Data/Protein/formatted_data/all_protein_names.csv") %>% dplyr::filter(from_id %in% rownames(data)))[["primary_gene_id"]]
#TRUE
'BROX' %in% (read.csv("Data/Protein/formatted_data/all_protein_names.csv") %>% dplyr::filter(from_id %in% rownames(validation_data)))[["primary_gene_id"]]
#FALSE

#i.e. BROX present in cohort 1 but not in cohort 2, so not present in common set of proteins


dep_sig <- read.table("DE_results_2024/proteomics/3_combat_corrected/p/sig_no_name_PREOPEVsHC.csv", sep = "\t", header = TRUE)
dep_sig.up <- dep_sig %>%
  dplyr::filter(logFC > 0)
dep_sig.down <- dep_sig %>%
  dplyr::filter(logFC < 0)

common.all <- Reduce(intersect, list(proteins.prev.plasmaEVs$Protein, 
                                     proteins.prev.urineEVs$Protein,
                                     dep_sig.up$Molecule))
ggvenn(list("Plasma EV GBM Vs HC" = proteins.prev.plasmaEVs$Protein, 
            "Urine EV GBM Vs HC" = proteins.prev.urineEVs$Protein,
            "Upregulated from our study" = dep_sig.up$Molecule),
       stroke_size = 0.1,
       set_name_size = 5,
       text_size = 3,
       fill_color = c("brown1", "navajowhite", "steelblue1")) +
  ggtitle("GBM Vs HC significant proteins & our study upregulated proteins") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.5)),
        plot.caption = element_text(hjust = 0.5))
ggsave(paste0(plot_dir_path, "4_PreviousStudiesVsOurStudyUp.jpg"), units = "cm", width = 24)


ggvenn(list("Plasma EV GBM Vs HC" = proteins.prev.plasmaEVs$Protein, 
            "Urine EV GBM Vs HC" = proteins.prev.urineEVs$Protein,
            "Downregulated from our study" = dep_sig.down$Molecule),
       stroke_size = 0.1,
       set_name_size = 5,
       text_size = 3,
       fill_color = c("brown1", "navajowhite", "steelblue1")) +
  ggtitle("GBM Vs HC significant proteins & our study downregulated proteins") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.5)),
        plot.caption = element_text(hjust = 0.5))
ggsave(paste0(plot_dir_path, "5_PreviousStudiesVsOurStudyDown.jpg"), units = "cm", width = 24)
