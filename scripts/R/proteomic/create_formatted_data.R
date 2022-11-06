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

