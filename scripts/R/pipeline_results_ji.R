library(tidyverse)
source("dataset_pipeline_arguments.R")
setwd("../../fem_pipeline_results/")

compute_jaccard_index_pairwise <- function(fsm1, fsm2, features_info, total_iter = 30){
  features_info_subset <- features_info %>%
    filter(FSM %in% c(fsm1, fsm2)) %>%
    select(-c(1,2))
  if (fsm1 == fsm2) {
    total_ji <- 0
    count <- 0
    # averaging over all possible pairs of iterations except same iteration to itself
    # not considering same iteration to itself when fsm1=fsm2, since this would be same value
    for (i in c(1:(total_iter-1))){
      for (j in c((i+1):total_iter)){
        sums <- colSums(features_info_subset[c(i, j), ])
        total_ji <- total_ji +
          (sum(sums == 2) / sum(sums != 0))
        count <- count + 1
      }
    } 
    ji <- total_ji / count
  }
  else {
    total_ji <- 0
    count <- 0
    # averaging over all possible pairs of iterations
    for (i in c(1:total_iter)){
      for (j in c(i:total_iter)){
        sums <- colSums(features_info_subset[c(i, j), ])
        total_ji <- total_ji +
          (sum(sums == 2) / sum(sums != 0))
        count <- count + 1
      }
    } 
    ji <- total_ji / count
  }
  return (ji)
}


compute_all_jaccard_index <- function(fsm_vector, features_info){
  ji_df <- data.frame()
  for(i in c(1:length(fsm_vector))){
    ji <- compute_jaccard_index_pairwise(fsm_vector[i], fsm_vector[i], features_info)
    ji_df <- rbind(ji_df, 
                       data.frame(FSM1 = fsm_vector[i], FSM2 = fsm_vector[i], JI = ji))        
  }
  return (ji_df)
}

print("Computing JI.....")

dparg_id_vec <- c(41:43, 125:127, 131, 133, 135, 137)

fsm_vector <- c("all", 
                "t-test", "t-test_BH",
                "t-test_pval_0.025", "t-test_pval_0.01", "t-test_pval_0.005",
                "wilcoxontest", "wilcoxontest_BH",
                "wilcoxontest_pval_0.025", "wilcoxontest_pval_0.001", "wilcoxontest_pval_0.005",
                "ranger_impu_cor", 
                "mrmr10", "mrmr20",
                "mrmr30", "mrmr50", 
                "mrmr75", "mrmr100",
                "RF_RFE", "ga_rf")

for(arg in dparg_id_vec){
  ds <- dataset_pipeline_arguments[[arg]]
  dataset_id <- paste(ds$dataset_id, ds$classification_criteria, sep = "_")
  print(dataset_id)  
  
  features_file <- paste(dataset_id, "features.csv", sep = "_")
  
  tryCatch({
    features_info <- read.table(features_file, sep = ',', header = TRUE)
    
    features_info <- features_info %>%
      filter(FSM %in% fsm_vector)
    
    ji_df <- compute_all_jaccard_index(fsm_vector, features_info)
    
    ji_df <- cbind(DataSetId = dataset_id, ji_df)
    
    dir_path <- 'JI'
    if (!dir.exists(dir_path)){
      dir.create(dir_path)
    }
    ji_data_file_name <- "all_ji.csv"
    file_path <- paste(dir_path, ji_data_file_name, sep = "/")
    write.table(ji_df, file = file_path,
                quote = FALSE, sep = ",", row.names = FALSE, append = TRUE,
                col.names = !file.exists(file_path))    
  },
  error = function(cond){
    print(cond)
  })
    
}

print("JI computation complete")

