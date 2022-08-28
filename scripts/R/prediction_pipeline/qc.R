library(tidyverse)

#train cohort
phenotype <- read.table("Data/transcriptomic_phenotype.txt", header=TRUE, sep="\t")

                    #burden              #resistance            #progression
comparison_vec <- c("PREOPEVsPOSTOPE_TP", "PREOPEVsREC_TP", "POSTOPE_TPVsREC_TP")

count_df <- data.frame(matrix(nrow = 0, ncol = 3, dimnames = list(c(),
                                                                  c("DataCohort", "Class", "Count")
                                                                  )))
count_row <- c()
for(comparison in comparison_vec){
  classes <- strsplit(comparison, split = "Vs")[[1]]
  output_labels <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes) %>%
    dplyr::select(Sample, Label)
  
  count <- summary(factor(output_labels$Label)) 

  count_row <- c("DataCohort" = "Train",
                 "Class" = classes[1],
                 "Count" = count[classes[1]])
  count_df[nrow(count_df) + 1, ] <- count_row
  
  count_row <- c("DataCohort" = "Train",
                 "Class" = classes[2],
                 "Count" = count[classes[2]])
  count_df[nrow(count_df) + 1, ] <- count_row
}

count_df <- count_df %>%
  unique()


#Test cohort (glionet cohort)
test_metadata <- read.table("Data/RNA_validation/metadata_glionet.csv", header = TRUE, sep = ",")
count_rows <- summary(factor(test_metadata$category_old_name))

for(i in c(1:length(count_rows))){
  count_row <- c("DataCohort" = "Test",
                 "Class" = names(count_rows)[i],
                 "Count" = count_rows[i])
  count_df[nrow(count_df) + 1, ] <- count_row
}

count_df <- count_df %>%
  arrange(Class)

write.table(count_df, "Data/validation_prediction_result/summary_count.csv", 
            sep = ",", row.names = FALSE, col.names = TRUE)