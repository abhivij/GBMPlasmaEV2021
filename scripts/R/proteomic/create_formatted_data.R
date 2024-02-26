library(tidyverse)
library(readxl)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

source("scripts/R/utils.R")
source("scripts/R/plot_data.R")

if(!dir.exists("Data/Protein/formatted_data/")){
  dir.create("Data/Protein/formatted_data/")
}

# input_file_path <- "Data/Protein/norm_output/norm_annotatedQ1-6_NA_FALSE.csv"
# max(formatted_data, na.rm = TRUE)
# 25.10626
# min(formatted_data, na.rm = TRUE)
# -2.664906

# input_file_path <- "Data/Protein/norm_output/norm_annotatedQ1-6_NA_equalizeMedians.csv"
# max(formatted_data, na.rm = TRUE)
# 33.42296
# min(formatted_data, na.rm = TRUE)
# -4.876739

# input_file_path <- "Data/Protein/norm_output/norm_annotatedQ1-6_NA_quantile.csv"
# max(formatted_data, na.rm = TRUE)
# 23.75979
# min(formatted_data, na.rm = TRUE)
# -8.416577


#create a csv file showing results of filtering proteins with specific percentages of NAs
#     so as to determine what filter_na_per to use
compare_filter_na_from_input <- function(input_file_path,
                                         output_filename = "Data/Protein/formatted_data/filter_na_result.csv"){
  protein_data <- read.csv(file = input_file_path)
  
  formatted_data <-  protein_data %>%
    select(-c(GROUP_ORIGINAL)) %>%
    column_to_rownames("SUBJECT_ORIGINAL")
  
  filter_na_result_df <- compare_filter_na(formatted_data)
  
  filepath_parts <- strsplit(input_file_path, split = "/", fixed = TRUE)[[1]]
  filename <- filepath_parts[length(filepath_parts)]
  
  filter_na_result_df[["filename"]] <- filename
  
  write.table(filter_na_result_df, output_filename, sep = ",", 
              row.names = FALSE, append = TRUE,
              col.names = !file.exists(output_filename))
}

compare_filter_na_from_input("Data/Protein/norm_output/norm_annotatedQ1-6_NA_FALSE.csv")
compare_filter_na_from_input("Data/Protein/norm_output/norm_annotatedQ7_NA_FALSE.csv")

# output_file_path <- "Data/Protein/formatted_data/Q1-6_nonorm_formatted.csv"


process_and_format_protein_data("Data/Protein/norm_output/norm_annotatedQ1-6_NA_FALSE.csv",
                                "Data/Protein/formatted_data/Q1-6_nonorm_formatted.csv")
process_and_format_protein_data("Data/Protein/norm_output/norm_annotatedQ1-6_NA_FALSE.csv",
                                "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv",
                                impute = TRUE, filter_na_perc = 50)
process_and_format_protein_data("Data/Protein/norm_output/norm_annotatedQ1-6_NA_FALSE.csv",
                                "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute75fil.csv",
                                impute = TRUE, filter_na_perc = 75)



process_and_format_protein_data("Data/Protein/norm_output/norm_annotatedQ7_NA_FALSE.csv",
                                "Data/Protein/formatted_data/Q7_nonorm_formatted.csv")
process_and_format_protein_data("Data/Protein/norm_output/norm_annotatedQ7_NA_FALSE.csv",
                                "Data/Protein/formatted_data/Q7_nonorm_formatted_impute50fil.csv",
                                impute = TRUE, filter_na_perc = 50)
process_and_format_protein_data("Data/Protein/norm_output/norm_annotatedQ7_NA_FALSE.csv",
                                "Data/Protein/formatted_data/Q7_nonorm_formatted_impute75fil.csv",
                                impute = TRUE, filter_na_perc = 75)


for_data1 <- read.table("Data/Protein/formatted_data/Q1-6_nonorm_formatted.csv", 
                        header=TRUE, sep=",", row.names=1, skip=0,
                   nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
for_data2 <- read.table("Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute75fil.csv", 
                        header=TRUE, sep=",", row.names=1, skip=0,
                        nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
for_data3 <- read.table("Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv", 
                        header=TRUE, sep=",", row.names=1, skip=0,
                        nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")


phenotype_info <- read.table("Data/proteomic_phenotype.txt", header=TRUE, sep="\t") %>%
  filter(GROUP_Q1to6 %in% c("PREOPE", "MET", "HC"))

for_data1 <- for_data1 %>%
  select(phenotype_info$Sample)
for_data2 <- for_data2 %>%
  select(phenotype_info$Sample)
for_data3 <- for_data3 %>%
  select(phenotype_info$Sample)

plot_data(t(for_data1), "static_NA_replace.png", "Static NA replace", 
          colour_label = "Label", groups = phenotype_info$GROUP_Q1to6, 
          shownames = FALSE, text = NA, dim_red = "umap")
plot_data(t(for_data2), "impute_NA_replace_after75fil.png", "75% filter & impute NA replace", 
          colour_label = "Label", groups = phenotype_info$GROUP_Q1to6, 
          shownames = FALSE, text = NA, dim_red = "umap")
plot_data(t(for_data3), "impute_NA_replace_after50fil.png", "50% filter & impute NA replace", 
          colour_label = "Label", groups = phenotype_info$GROUP_Q1to6, 
          shownames = FALSE, text = NA, dim_red = "umap")



################### new cohort
compare_filter_na_from_input("Data/Protein/norm_output/norm__newcohort_processed_NA_FALSE.csv",
                             output_filename = "Data/Protein/formatted_data/newcohort_filter_na_result.csv")


##################################################################
### create data files with common proteins 

data <- read.csv("Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv", 
                 row.names=1)
validation_data <- read.csv("Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil.csv", 
                            row.names = 1)
length(intersect(rownames(data), rownames(validation_data)))
common_proteins <- intersect(rownames(data), rownames(validation_data))

initial_data_common <- data[common_proteins, ]
validation_data_common <- validation_data[common_proteins, ]

write.csv(initial_data_common, 
          "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil_common.csv")
write.csv(validation_data_common, 
          "Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil_common.csv")

##################################################################
###create validation data with proteins from train data

data <- read.csv("Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv", 
                 row.names=1)
validation_data <- read.csv("Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil.csv", 
                            row.names = 1)
length(intersect(rownames(data), rownames(validation_data)))
common_proteins <- intersect(rownames(data), rownames(validation_data))
initial_cohort_specific_proteins <- setdiff(rownames(data), rownames(validation_data))
initial_cohort_proteins <- rownames(data)

validation_data_common <- validation_data[common_proteins, ]
validation_data.2 <-  data.frame(matrix(data = NA, nrow = length(initial_cohort_specific_proteins),
                             ncol = ncol(validation_data),
                             dimnames = list(initial_cohort_specific_proteins,
                                             colnames(validation_data)))
                             )

validation_data_initial_cohort_proteins <- rbind(validation_data_common, validation_data.2)

combined_data <- cbind(validation_data_initial_cohort_proteins[rownames(data), ],
                       data)
write.csv(combined_data, 
          "Data/Protein/formatted_data/combined_data_with_initial_cohort_proteins.csv")

combined_imputed_data <- read.csv("Data/Protein/formatted_data/combined_data_with_initial_cohort_proteins_imputed.csv",
                                  row.names = 1)
validation_cohort_with_initial_cohort_proteins <- combined_imputed_data[, colnames(validation_data)]
write.csv(validation_cohort_with_initial_cohort_proteins, 
          "Data/Protein/formatted_data/validation_cohort_with_initial_cohort_proteins.csv")




##################################################################
# rename samples in validation cohort as done in prediction_pipeline so as to 
#   run this in FEMPipeline

validation_data <- read.csv("Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil_common.csv",
                                   row.names = 1)
colnames(validation_data)[colnames(validation_data) == "SB12_01"] = "SB12"
#use SB22.02
colnames(validation_data)[colnames(validation_data) == "SB22.02"] = "SBtobeused22"
colnames(validation_data)[colnames(validation_data) == "SB22"] = "SB22_dont_include"
colnames(validation_data)[colnames(validation_data) == "SBtobeused22"] = "SB22"

write.csv(validation_data, 
          "Data/Protein/formatted_data/newcohort_common_correctedsamples.csv")




#########################
#create file with protein names - for all proteins in our study
library(UniProt.ws)

data <- read.csv("Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv", 
                 row.names=1)
validation_data <- read.csv("Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil.csv", 
                            row.names = 1)

all_proteins_in_study <- union(rownames(data), rownames(validation_data))
protein_name_df <- data.frame(id = all_proteins_in_study)

up <- UniProt.ws(taxId=9606)
species(up)
keytypes(up)

head(columns(up))

result <- select(
  x = up,
  keys = c(all_proteins_in_study),
  columns = c("protein_name", "gene_primary", "organism_id"),
  keytype = "UniProtKB"
)
colnames(result) <- c("from_id", "uniprot_id", "protein_name", "primary_gene_id", "organism_id")

#just to check if evrything is 9606 - i.e. human -ideally should be the case
result %>% filter(is.na(organism_id) | organism_id != 9606)
# from  entry protein_name primary_gene_name organism_id
# 1 E7EML9 E7EML9      deleted              <NA>        <NA>

result %>% filter(from == "E7EML9")
# from  entry protein_name primary_gene_name organism_id
# 1 E7EML9 E7EML9      deleted              <NA>        <NA>

#so there is some other duplicate entry
result %>% group_by(from) %>% summarize(count = n()) %>% filter(count > 1)
# # A tibble: 1 × 2
# from   count
# <chr>  <int>
#   1 Q6ZMK1     2

result %>% filter(from == "Q6ZMK1")
# from  entry                                                                       protein_name
# 1 Q6ZMK1 P0DTL5                                                          Transmembrane protein 276
# 2 Q6ZMK1 P0DTL6 Zinc finger TRAF-type-containing protein 1 (Cysteine and histidine-rich protein 1)
# primary_gene_name organism_id
# 1           TMEM276        9606
# 2           ZFTRAF1        9606

#not sure what this means - will check if this protein is identified as biomarker

write.csv(result %>% dplyr::select(-c(organism_id)), "Data/Protein/formatted_data/all_protein_names.csv",
          row.names = FALSE)


#######################

#create combined data with PREOPE, MET, HC samples from both cohorts specifically for DE analysis

data_file_path <- "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv"
validation_data_file_path <- "Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil.csv" 
phenotype_file_path <- "Data/proteomic_phenotype_PREOPE_MET_HC_withaddicolumn.txt"
output_dir_path <- "Data/Protein"

data <- read.csv(data_file_path, row.names = 1)
validation_data <- read.csv(validation_data_file_path, row.names = 1)
phenotype <- read.table(phenotype_file_path, header=TRUE, sep="\t")

colnames(validation_data)[colnames(validation_data) == "SB12_01"] = "SB12"
#use SB22.02
colnames(validation_data)[colnames(validation_data) == "SB22.02"] = "SBtobeused22"
colnames(validation_data)[colnames(validation_data) == "SB22"] = "SB22_dont_include"
colnames(validation_data)[colnames(validation_data) == "SBtobeused22"] = "SB22"

phenotype <- phenotype %>%
  dplyr::rename(c("Label" = PREOPE_MET_HC))

output_labels.cohort1 <- phenotype %>%
  dplyr::filter(!is.na(Label)) %>%
  dplyr::select(Sample, Label, data_cohort, Subgroup, Sex, Age) %>%
  dplyr::filter(data_cohort == "initial")
output_labels.cohort2 <- phenotype %>%
  dplyr::filter(!is.na(Label)) %>%
  dplyr::select(Sample, Label, data_cohort, Subgroup, Sex, Age) %>%
  dplyr::filter(data_cohort == "validation")

#currently data format : (transcripts x samples)

data.cohort1 <- data %>% dplyr::select(output_labels.cohort1$Sample)
data.cohort2 <- validation_data %>% dplyr::select(output_labels.cohort2$Sample)

common <- intersect(rownames(data.cohort1), rownames(data.cohort2))  
data.cohort1 <- data.cohort1[common, ]
data.cohort2 <- data.cohort2[common, ]
data <- cbind(data.cohort1, data.cohort2)    

output_labels <- rbind(output_labels.cohort1, output_labels.cohort2)

write.csv(data, "Data/Protein/formatted_data/PREOPE_MET_HC_data.csv")
