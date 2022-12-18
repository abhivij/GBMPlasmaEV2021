library(missForest)
library(tidyverse)
library(ggrepel)
library(umap)
library(ggvenn)
library(sva)

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

comparison = "PREOPEVsPOSTOPE_TP"
omics_type = "proteomic"
classes = c("POSTOPE_TP", "PREOPE")
data_to_show = c("initial", "validation", "both")
data_to_show = "initial"
show_only_common = TRUE
show_imputed_based_on_initial = FALSE
perform_filter = FALSE
dim_red = "UMAP"
norm = "quantile_train_param"
norm <- "none"

create_dim_red_plots <- function(comparison, classes,
                                 omics_type = "proteomics",
                                 dim_red, norm,
                                 data_to_show,
                                 show_only_common = FALSE,
                                 show_imputed_based_on_initial = FALSE,
                                 shownames = FALSE,
                                 perform_filter = TRUE){
  if(omics_type == "proteomics"){
    data_file_path <- "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv"
    phenotype_file_path <- "Data/proteomic_phenotype.txt"
    
    if(show_imputed_based_on_initial){
      validation_data_file_path <- "Data/Protein/formatted_data/validation_cohort_with_initial_cohort_proteins.csv"
    } else{
      validation_data_file_path <- "Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil.csv"      
    }
    validation_metadata_file_path <- "Data/RNA_validation/metadata_glionet.csv"
    
  } else if(omics_type == "transcriptomics"){
    #to add later
    data_file_path <- ""
  }

  data <- read.csv(data_file_path, row.names = 1)
  phenotype <- read.table(phenotype_file_path, header=TRUE, sep="\t")
  
  validation_data <- read.csv(validation_data_file_path, row.names = 1)
  validation_metadata <- read.csv(validation_metadata_file_path) %>%
    rename("Label" = "category_old_name", "Sample" = "sample_id") %>%
    select(Sample, age, gender, Label)
  
  colnames(validation_data)[colnames(validation_data) == "SB12_01"] = "SB12"
  
  #use SB22.02
  colnames(validation_data)[colnames(validation_data) == "SB22.02"] = "SBtobeused22"
  colnames(validation_data)[colnames(validation_data) == "SB22"] = "SB22_dont_include"
  colnames(validation_data)[colnames(validation_data) == "SBtobeused22"] = "SB22"
  
  validation_metadata <- validation_metadata %>%
    filter(Sample != "SB7")
  
  title <- paste0(dim_red, " plot of ", paste(classes, collapse = ", "), " samples ", norm, " ",
                  data_to_show, " data")
  if(show_only_common){
    title <- paste(title, "common")
  }
  if(show_imputed_based_on_initial){
    title <- paste(title, "imputed with initial prot")
  }
  
  phenotype <- phenotype %>%
    mutate(cohort = "initial")
  validation_metadata <- validation_metadata %>%
    mutate(cohort = "validation")
  
  output_labels.initial <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes) %>%
    dplyr::select(Sample, Label, cohort)
  output_labels.validation <- validation_metadata %>%
    filter(Label %in% classes) %>%
    dplyr::select(Sample, Label, cohort)
  
  if(show_only_common){
    common_proteins <- intersect(rownames(data), rownames(validation_data))  
    data <- data[common_proteins, ]
    validation_data <- validation_data[common_proteins, ]
  }
  if(data_to_show == "initial"){
    data_of_interest <- data
    output_labels <- output_labels.initial
  } else if(data_to_show == "validation"){
    data_of_interest <- validation_data
    output_labels <- output_labels.validation
  } else if(data_to_show == "both"){
    data_of_interest <- cbind(data, validation_data)
    output_labels <- rbind(output_labels.initial, output_labels.validation)
  }
  output_labels <- output_labels %>%
    mutate(Label = factor(Label), cohort = factor(cohort))
  
  group_counts <- output_labels %>%
    dplyr::mutate(Label = paste(cohort, Label, sep = "_")) %>%
    group_by(Label) %>%
    summarise(n = n())
  
  group_counts_text <- paste(apply(group_counts, MARGIN = 1, FUN = function(x){paste(x[1], x[2], sep = ":")}),
                             collapse = "  ")
  data_of_interest <- data_of_interest[, output_labels$Sample]
  data <- data_of_interest
  
  #currently data format : (transcripts x samples)
  if(perform_filter){
    keep <- edgeR::filterByExpr(data, group = output_labels$Label)
    data <- data[keep, ]
  }

  
  if(norm == "quantile_train_param"){
    #adapted from https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
    data.rank <- apply(data, 2, rank, ties.method="average")
    data.sorted <- data.frame(apply(data, 2, sort))
    data.mean <- apply(data.sorted, 1, mean)
    index_to_mean <- function(index, data_mean){
      #index can be int or int+0.5
      #if int+0.5, take average of the numbers in those positions
      int.result <- data_mean[index]
      index.int <- floor(index)
      #some of the values in point5.result might be NA
      #but they won't be chosen
      point5.result <- (data_mean[index.int] + data_mean[index.int+1])/2
      point5.indices <- index%%1 != 0
      result <- int.result
      result[point5.indices] <- point5.result[point5.indices]
      return (result)
    }
    data.norm <- apply(data.rank, 2, index_to_mean, data_mean = data.mean)
    rownames(data.norm) <- rownames(data)
    data <- data.norm
  }
  data <- as.data.frame(t(as.matrix(data)))
  
  if(shownames){
    text <- rownames(data)
  } else{
    text <- ""
  }
  set.seed(1)
  if(dim_red == "PCA"){
    result <- prcomp(data)
    dim_red_df <- data.frame(x = result$x[,1], y = result$x[,2])    
    xlab <- "PCA 1"
    ylab <- "PCA 2"  
  } else if(dim_red == "UMAP"){
    result <- umap(data)
    dim_red_df <- data.frame(x = result$layout[,1], y = result$layout[,2])  
    xlab <- "UMAP 1"
    ylab <- "UMAP 2"
  }
  
  ggplot2::ggplot(dim_red_df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(colour = output_labels$Label,
                                     shape = output_labels$cohort), size = 3) +
    # geom_text_repel(aes(label = text)) +
    ggplot2::labs(title = title) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    labs(caption = paste(paste("Data dimension :", paste(dim(data), collapse = "x")), "\n",
                         group_counts_text),
         colour = "Label",
         shape = "Cohort")
  
  dir_path <- "plots/qc/dim_red/"
  if(!dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE)
  }
  file_name <- paste0(gsub(title, pattern = " |,", replacement = "-"), ".jpg")
  file_path <- paste(dir_path, file_name, sep = "/")
  ggplot2::ggsave(file_path, units = "cm", width = 30)
}
