library(missForest)

combined_data <- read.csv("Data/Protein/formatted_data/combined_data_with_initial_cohort_proteins.csv",
                          row.names = 1)
combined_data.impute <- t(missForest(t(combined_data))$ximp)
write.csv(combined_data.impute,
          "Data/Protein/formatted_data/combined_data_with_initial_cohort_proteins_imputed.csv")

