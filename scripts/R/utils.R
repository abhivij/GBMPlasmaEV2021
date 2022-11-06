library(missForest)

append_path <- function(dirname, filename){
  if (dirname != "") {
    filename <- paste(dirname, filename, sep = "/")
  }
  return (filename)
}

# filter_na_perc <- 50
# data <- formatted_data
impute_data <- function(data, filter_na_perc = 75, impute = TRUE){
  data <- data[, colSums(is.na(data)) < (filter_na_perc/100*nrow(data))]
  print(sum(is.na(data)))
  if (impute) {
    data <- missForest(data)$ximp
  }
  return (data)
}

compare_filter_na <- function(data){
  filter_perc_vec <- c(100, 90, 80, 75, 50, 25, 20, 10, 5)
  info_df <- data.frame()
  row_info <- c(dim(data)[1], "-", dim(data)[2], 100)
  info_df <- rbind(info_df, row_info)
  colnames(info_df) <- c("num_samples", "NA_filter_percent", "num_cols", "col_percent")
  for(perc in filter_perc_vec){
    print(perc)
    filtered_data <- data[, colSums(is.na(data)) < (perc/100*nrow(data))]
    row_info <- c(dim(filtered_data)[1], perc, dim(filtered_data)[2], 
                  (dim(filtered_data)[2]*100/dim(data)[2]))
    info_df <- rbind(info_df, row_info)
  }
  return (info_df)
}




# comparison_list <- list(c("PREOPE", "MET"), c("PREOPE", "HC"), c("MET", "HC"))
# class_column_name <- "GROUP_Q1to6"
insert_comparison_columns <- function(phenotype_info, comparison_list, class_column_name){
  #below code is to generalize
  # mutate(PREOPEVsMET = case_when(
  #   GROUP_Q1to6 == "PREOPE" ~ "PREOPE",
  #   GROUP_Q1to6 == "MET" ~ "MET",
  #   TRUE ~ NA_character_))     
  for(comparison in comparison_list){
    comparison_column_name <- paste0(comparison[1], "Vs", comparison[2])
    phenotype_info <- phenotype_info %>%
      mutate("{comparison_column_name}" := case_when(
        !!sym(class_column_name) == comparison[1] ~ comparison[1],
        !!sym(class_column_name) == comparison[2] ~ comparison[2],
        TRUE ~ NA_character_
      ))
  }
  phenotype_info
}


# fsm_vec <- c("all", 
#              "t-test", "wilcoxontest",
#              "ranger_impu_cor", 
#              "mrmr10", "mrmr20",
#              "mrmr30", "mrmr50", 
#              "mrmr75", "mrmr100",
#              "RF_RFE", "ga_rf")

fsm_vector <- c("all", 
                "t-test", "wilcoxontest",
                "ranger_impu_cor",
                "ranger_pos_impu_cor",
                "mrmr10", "mrmr20",
                "mrmr30", "mrmr50", 
                "mrmr75", "mrmr100",
                "mrmr_perc50",
                "RF_RFE",
                "ga_rf")



process_and_format_protein_data <- function(input_file_path, output_file_path,
                                            impute = FALSE, filter_na_perc = 75){
  protein_data <- read.csv(file = input_file_path)
  
  protein_data <- protein_data %>%
    arrange(SUBJECT_ORIGINAL)
  
  formatted_data <-  protein_data %>%
    select(-c(GROUP_ORIGINAL)) %>%
    column_to_rownames("SUBJECT_ORIGINAL")
  
  print(max(formatted_data, na.rm = TRUE))
  print(min(formatted_data, na.rm = TRUE))
  
  print(sum(is.na(formatted_data)))
  if(impute){
    #impute = TRUE does imputation using missForest
    formatted_data <- impute_data(formatted_data, 
                                  filter_na_perc = filter_na_perc,
                                  impute = TRUE)
  }else{
    #impute = FALSE replaces NA with (min - 25th_quantile)
    quantiles <- quantile(formatted_data, na.rm = TRUE)
    na_repl_value <- quantiles["0%"] - quantiles["25%"]
    formatted_data[is.na(formatted_data)] <- na_repl_value
  }
  print(sum(is.na(formatted_data)))
  
  
  formatted_data <- t(formatted_data)
  write.csv(formatted_data, output_file_path)
}
