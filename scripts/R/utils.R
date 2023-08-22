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


show_imputed_based_on_initial = FALSE
shownames = FALSE
perform_filter = TRUE
batch_effect_correction = "none"
plot_dir_path = "plots/qc/dim_red/"
best_features_file_path = NA
dataset_replace_str = NA

comparison = "POSTOPE_TPVsREC_TP"
omics_type = "transcriptomics"
classes = c("POSTOPE_TP", "REC_TP")
data_to_show = "both"
show_only_common = TRUE
show_imputed_based_on_initial = FALSE
perform_filter = FALSE
dim_red = "UMAP"
norm = "log_cpm"
batch_effect_correction = "combat"
plot_dir_path = "plots/qc/dim_red/tr_combat/"

create_dim_red_plots <- function(comparison, classes,
                                 omics_type = "proteomics",
                                 dim_red, norm,
                                 data_to_show,
                                 show_only_common = FALSE,
                                 show_imputed_based_on_initial = FALSE,
                                 shownames = FALSE,
                                 perform_filter = TRUE,
                                 batch_effect_correction = "none",
                                 plot_dir_path = "plots/qc/dim_red/",
                                 best_features_file_path = NA,
                                 dataset_replace_str = NA
                                 ){
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
    data_file_path <- "Data/RNA/umi_counts_initial_cohort.csv"
    phenotype_file_path <- "Data/transcriptomic_phenotype.txt"
    validation_data_file_path <- "Data/RNA/umi_counts_validation_cohort.csv"      
    validation_metadata_file_path <- "Data/RNA_validation/metadata_glionet.csv"

    # validation_metadata <- read.csv("Data/RNA_validation/metadata_glionet.csv") %>%
    #   mutate(sample_category = factor(sample_category)) %>%
    #   mutate(sample_category = recode_factor(sample_category, "PRE-OP" = "PREOPE",
    #                                          "POST-OP" = "POSTOPE_TP",
    #                                          "RECURRENCE" = "REC_TP"))
  }

  data <- read.csv(data_file_path, row.names = 1)
  phenotype <- read.table(phenotype_file_path, header=TRUE, sep="\t")
  
  validation_data <- read.csv(validation_data_file_path, row.names = 1)
  validation_metadata <- read.csv(validation_metadata_file_path) %>%
    rename("Label" = "category_old_name", "Sample" = "sample_id") %>%
    select(Sample, age, gender, Label)
  
  if(omics_type == "proteomics"){
    colnames(validation_data)[colnames(validation_data) == "SB12_01"] = "SB12"
    #use SB22.02
    colnames(validation_data)[colnames(validation_data) == "SB22.02"] = "SBtobeused22"
    colnames(validation_data)[colnames(validation_data) == "SB22"] = "SB22_dont_include"
    colnames(validation_data)[colnames(validation_data) == "SBtobeused22"] = "SB22"
    
    validation_metadata <- validation_metadata %>%
      filter(Sample != "SB7")    
  } else if(omics_type == "transcriptomics"){
    colnames(validation_data) <- paste0("S", colnames(validation_data))
  }
  
  title <- paste0(dim_red, " plot of ", omics_type, " ",
                  paste(classes, collapse = ", "), " samples ", norm, " ",
                  data_to_show, " data")
  if(show_only_common){
    title <- paste(title, "common")
  }
  if(show_imputed_based_on_initial){
    title <- paste(title, "imputed with initial prot")
  }
  title <- paste(title, batch_effect_correction)
  
  phenotype <- phenotype %>%
    mutate(cohort = "initial")
  validation_metadata <- validation_metadata %>%
    mutate(cohort = "validation")
  
  output_labels.train <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes) %>%
    dplyr::select(Sample, Label, cohort)
  output_labels.test <- validation_metadata %>%
    filter(Label %in% classes) %>%
    dplyr::select(Sample, Label, cohort)
  
  
  #currently data format : (transcripts x samples)
  
  data.train <- data %>% dplyr::select(output_labels.train$Sample)
  data.test <- validation_data %>% dplyr::select(output_labels.test$Sample)
  
  if(show_only_common){
    common <- intersect(rownames(data.train), rownames(data.test))  
    data.train <- data.train[common, ]
    data.test <- data.test[common, ]
  }
  
  ################obtain best biomarkers
  if(!is.na(best_features_file_path) & !is.na(dataset_replace_str)){
    
    #note : performing t(), selecting based on columns, and then again t()
    #       is done instead of just slecting based on rows, so as to keep the code same
    #       as in pipeline_validation_data.R
    
    data.train <- as.data.frame(t(as.matrix(data.train)))
    data.test <- as.data.frame(t(as.matrix(data.test)))
    
    best_features <- read.csv(best_features_file_path)  
    best_features_sub <- best_features %>%
      mutate(dataset_id = gsub(dataset_replace_str, "", dataset_id)) %>%
      filter(is_best == 1, dataset_id == comparison)
    biomarkers <- strsplit(best_features_sub$biomarkers, split = "|", fixed = TRUE)[[1]]  
    
    data.train <- data.frame(data.train)[, biomarkers]  #data.frame() replaces - in colnames to .
    colnames(data.test) <- gsub("-", ".", colnames(data.test))
    sum(biomarkers %in% colnames(data.test))
    
    available_biomarkers <- c(biomarkers[biomarkers %in% colnames(data.test)])
    print(available_biomarkers)
    non_available_biomarkers <- biomarkers[!biomarkers %in% colnames(data.test)] 
    print(non_available_biomarkers)
    
    print(length(available_biomarkers))
    print(length(non_available_biomarkers))
    
    data.test <- data.test %>%
      select(available_biomarkers)
    for(nab in non_available_biomarkers){
      data.test[[nab]] <- 0
    }    
    data.train <- as.data.frame(t(as.matrix(data.train)))
    data.test <- as.data.frame(t(as.matrix(data.test)))
  }
  ################obtain best biomarkers end
  
  if(perform_filter){
    keep <- edgeR::filterByExpr(data.train, group = output_labels.train$Label)
    data.train <- data.train[keep, ]
    data.test <- data.test[keep, ]  
  }

  if(norm == "quantile_train_param"){
    #adapted from https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
    data.train.rank <- apply(data.train, 2, rank, ties.method="average")
    data.train.sorted <- data.frame(apply(data.train, 2, sort))
    data.train.mean <- apply(data.train.sorted, 1, mean)
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
    data.train.norm <- apply(data.train.rank, 2, index_to_mean, data_mean = data.train.mean)
    rownames(data.train.norm) <- rownames(data.train)
    data.train <- data.train.norm
    
    data.test.rank <- apply(data.test, 2, rank, ties.method="average")
    #use params i.e. mean values of rows, from training data
    data.test.norm <- apply(data.test.rank, 2, index_to_mean, data_mean = data.train.mean)
    rownames(data.test.norm) <- rownames(data.test)
    data.test <- data.test.norm
  } else if(norm == "log_cpm"){
    #calculating norm log cpm
    data.train <- edgeR::cpm(data.train, log=TRUE)
    data.test <- edgeR::cpm(data.test, log=TRUE)
  }
  data.train <- as.data.frame(t(as.matrix(data.train)))
  data.test <- as.data.frame(t(as.matrix(data.test)))
  
  if(data_to_show == "initial"){
    data_of_interest <- data.train
    output_labels <- output_labels.train
  } else if(data_to_show == "validation"){
    data_of_interest <- data.test
    output_labels <- output_labels.test
  } else if(data_to_show == "both"){
    data_of_interest <- rbind(data.train, data.test)
    output_labels <- rbind(output_labels.train, output_labels.test)
  }
  output_labels <- output_labels %>%
    mutate(Label = factor(Label), cohort = factor(cohort))

  group_counts <- output_labels %>%
    dplyr::mutate(Label = paste(cohort, Label, sep = "_")) %>%
    group_by(Label) %>%
    summarise(n = n())

  group_counts_text <- paste(apply(group_counts, MARGIN = 1, FUN = function(x){paste(x[1], x[2], sep = ":")}),
                             collapse = "  ")
  
  all.equal(rownames(data_of_interest), output_labels$Sample)
  data <- data_of_interest
  
  if(batch_effect_correction == "combat"){
    data <- as.data.frame(t(as.matrix(data)))
    data.combat = ComBat(dat=data, batch=output_labels$cohort)
    data.combat <- as.data.frame(t(as.matrix(data.combat)))
    data <- data.combat
  } else if(batch_effect_correction == "combat_ref"){
    data <- as.data.frame(t(as.matrix(data)))
    data.combat = ComBat(dat=data, batch=output_labels$cohort, ref.batch = 'initial')
    data.combat <- as.data.frame(t(as.matrix(data.combat)))
    
    data <- as.data.frame(t(as.matrix(data)))
    data.train <- data[output_labels.train$Sample, ]
    data.train.combat <- data.combat[output_labels.train$Sample, ]
    print(all.equal(data.train, data.train.combat))
    # TRUE
    # 
    # data.test <- data[output_labels.test$Sample, ]
    # data.test.combat <- data.combat[output_labels.test$Sample, ]
    # print(all.equal(data.test, data.test.combat))
    # FALSE
    data <- data.combat
  }
  
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
    geom_text_repel(aes(label = text)) +
    ggplot2::labs(title = title) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    labs(caption = paste(paste("Data dimension :", paste(dim(data), collapse = "x")), "\n",
                         group_counts_text),
         colour = "Label",
         shape = "Cohort")
  
  if(!dir.exists(plot_dir_path)){
    dir.create(plot_dir_path, recursive = TRUE)
  }
  file_name <- paste0(gsub(title, pattern = " |,", replacement = "-"), ".jpg")
  file_path <- paste(plot_dir_path, file_name, sep = "/")
  ggplot2::ggsave(file_path, units = "cm", width = 30)
}



#create dim red for PREOPE MET HC

dim_red = "UMAP"
shownames = FALSE
perform_filter = TRUE
batch_effect_correction = "none"
plot_dir_path = "plots_PREOPE_MET_HC/qc/dim_red/"
best_features_file_path = NA
dataset_replace_str = NA

comparison = "PREOPEVsMET"
classes = c("MET", "PREOPE")
omics_type = "transcriptomics"
norm = "log_cpm"
dim_red = "UMAP"
shownames = FALSE
perform_filter = TRUE
batch_effect_correction = "combat"
plot_dir_path = "plots_comparison_set2/qc/dim_red/"
best_features_file_path = NA
dataset_replace_str = NA
file_name_prefix = 11

best_features_file_path = "Data/selected_features/best_features_with_add_col.csv"
dataset_replace_str = "GBM_combined_transcriptomic_combat_compset2_"

boxplot_dir_path = "plots_comparison_set2/qc/boxplot/"


comparison = "PREOPEVsMET"
classes = c("MET", "PREOPE")
omics_type = "proteomics"
norm = "quantile_train_param"
dim_red = "UMAP"
shownames = FALSE
perform_filter = FALSE
batch_effect_correction = "combat"
plot_dir_path = "plots_comparison_set2/qc/dim_red_best/"
file_name_prefix = 1
best_features_file_path = "Data/selected_features/best_features_with_add_col.csv"
dataset_replace_str = "GBM_combined_proteomic_combat_compset2_"
boxplot_dir_path = "plots_comparison_set2/qc/boxplot/"

create_dim_red_plots_PMH <- function(comparison, classes,
                                     omics_type, norm,
                                     dim_red = "UMAP",
                                     shownames = FALSE,
                                     perform_filter = TRUE,
                                     batch_effect_correction = "none",
                                     plot_dir_path = "plots_PREOPE_MET_HC/qc/dim_red/",
                                     boxplot_dir_path = "plots_PREOPE_MET_HC/qc/boxplot/",
                                     file_name_prefix = "",
                                     best_features_file_path = NA,
                                     dataset_replace_str = NA){
  if(omics_type == "proteomics"){
    data_file_path <- "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv"
    validation_data_file_path <- "Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil.csv" 
    phenotype_file_path <- "Data/proteomic_phenotype_PREOPE_MET_HC.txt"
  } else if(omics_type == "transcriptomics"){
    data_file_path <- "Data/RNA/umi_counts_initial_cohort.csv"
    validation_data_file_path <- "Data/RNA/umi_counts_validation_cohort.csv"      
    phenotype_file_path <- "Data/transcriptomic_phenotype_PREOPE_MET_HC.txt"
    
    # validation_metadata <- read.csv("Data/RNA_validation/metadata_glionet.csv") %>%
    #   mutate(sample_category = factor(sample_category)) %>%
    #   mutate(sample_category = recode_factor(sample_category, "PRE-OP" = "PREOPE",
    #                                          "POST-OP" = "POSTOPE_TP",
    #                                          "RECURRENCE" = "REC_TP"))
  }
  
  data <- read.csv(data_file_path, row.names = 1)
  validation_data <- read.csv(validation_data_file_path, row.names = 1)
  
  phenotype <- read.table(phenotype_file_path, header=TRUE, sep="\t")
  
  if(omics_type == "proteomics"){
    colnames(validation_data)[colnames(validation_data) == "SB12_01"] = "SB12"
    #use SB22.02
    colnames(validation_data)[colnames(validation_data) == "SB22.02"] = "SBtobeused22"
    colnames(validation_data)[colnames(validation_data) == "SB22"] = "SB22_dont_include"
    colnames(validation_data)[colnames(validation_data) == "SBtobeused22"] = "SB22"

    #SB7 sample not required for this analysis
    # so no need to filter out
    # validation_metadata <- validation_metadata %>%
    #   filter(Sample != "SB7")    
    
  } else if(omics_type == "transcriptomics"){
    colnames(validation_data) <- paste0("S", colnames(validation_data))
  }
  
  title <- paste0(dim_red, " plot of ", omics_type, " ",
                  paste(rev(classes), collapse = ", "), " samples ", 
                  norm, " data")
  
  title <- paste0(title, " BEC-", batch_effect_correction)
  
  output_labels.cohort1 <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes) %>%
    dplyr::select(Sample, Label, data_cohort, Subgroup, Sex, Age) %>%
    filter(data_cohort == "initial")
  output_labels.cohort2 <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes) %>%
    dplyr::select(Sample, Label, data_cohort, Subgroup, Sex, Age) %>%
    filter(data_cohort == "validation")
  
  #currently data format : (transcripts x samples)
  
  data.cohort1 <- data %>% dplyr::select(output_labels.cohort1$Sample)
  data.cohort2 <- validation_data %>% dplyr::select(output_labels.cohort2$Sample)
  
  #take common only if both cohorts contain samples
  if(ncol(data.cohort1) > 0 & ncol(data.cohort2) > 0){
    common <- intersect(rownames(data.cohort1), rownames(data.cohort2))  
    data.cohort1 <- data.cohort1[common, ]
    data.cohort2 <- data.cohort2[common, ]
    data <- cbind(data.cohort1, data.cohort2)    
  } else if(ncol(data.cohort1) > 0){
    data <- data.cohort1
  } else if(ncol(data.cohort2) > 0){
    data <- data.cohort2
  } else{
    print("no samples")
  }

  output_labels <- rbind(output_labels.cohort1, output_labels.cohort2)
  
  if(perform_filter){
    keep <- edgeR::filterByExpr(data, group = output_labels$Label)
    data_left_out <- data[!keep, ]
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
    
  } else if(norm == "log_cpm"){
    #calculating norm log cpm
    data <- edgeR::cpm(data, log=TRUE)
  }
  
  data <- as.data.frame(t(as.matrix(data)))
  
  output_labels <- output_labels %>%
    mutate(Label = factor(Label), data_cohort = factor(data_cohort))
  
  group_counts <- output_labels %>%
    dplyr::mutate(Label = paste(data_cohort, Label, sep = "_")) %>%
    group_by(Label) %>%
    summarise(n = n())
  
  group_counts_text <- paste(apply(group_counts, MARGIN = 1, FUN = function(x){paste(x[1], x[2], sep = ":")}),
                             collapse = "  ")
  
  all.equal(rownames(data), output_labels$Sample)
  
  if(batch_effect_correction == "combat"){
    data <- as.data.frame(t(as.matrix(data)))
    data.combat = ComBat(dat=data, batch=output_labels$data_cohort)
    data.combat <- as.data.frame(t(as.matrix(data.combat)))
    data <- data.combat
  } else if(batch_effect_correction == "combat_ref"){
    data <- as.data.frame(t(as.matrix(data)))
    data.combat = ComBat(dat=data, batch=output_labels$cohort, ref.batch = 'initial')
    data.combat <- as.data.frame(t(as.matrix(data.combat)))
    
    data <- as.data.frame(t(as.matrix(data)))
    data.cohort1 <- data[output_labels.cohort1$Sample, ]
    data.cohort1.combat <- data.combat[output_labels.cohort1$Sample, ]
    print(all.equal(data.cohort1, data.cohort1.combat))
    # TRUE
    # 
    # data.cohort2 <- data[output_labels.cohort2$Sample, ]
    # data.cohort2.combat <- data.combat[output_labels.cohort2$Sample, ]
    # print(all.equal(data.cohort2, data.cohort2.combat))
    # FALSE
    data <- data.combat
  }
  
  
  ################obtain best biomarkers
  
  if(!is.na(best_features_file_path) & !is.na(dataset_replace_str)){

    best_features <- read.csv(best_features_file_path)  
    best_features_sub <- best_features %>%
      mutate(dataset_id = gsub(dataset_replace_str, "", dataset_id)) %>%
      filter(is_best > 0, dataset_id == comparison)
    
    if(!dir.exists(boxplot_dir_path)){
      dir.create(boxplot_dir_path, recursive = TRUE)
    }
    all_protein_names <- read.csv("Data/Protein/formatted_data/all_protein_names.csv")
    
    #plot biomarker boxplot for the multiple biomarker sets identified
    for(i in c(1:nrow(best_features_sub))){
      # i <- 1
      biomarkers <- strsplit(best_features_sub[i, "biomarkers"], split = "|", fixed = TRUE)[[1]] 
      data_sub <- data[, biomarkers]
      
      #plot biomarker boxplot
      
      #args - data, biomarkers, classes, output_labels, omics_type
      
      data_to_plot <- data_sub %>%
        rownames_to_column(var = "Sample") %>%
        pivot_longer(cols = !Sample, names_to = "biomarker", values_to = "norm_expr") %>%
        inner_join(output_labels %>%
                     dplyr::select(c(Sample, Label)))
    
      if(omics_type == "transcriptomics"){
        x_lab <- "transcripts"
        y_lab <- "Log CPM expression"
      } else{
        x_lab <- "proteins"
        y_lab <- "Quantile normalized expression"
        
        data_to_plot <- data_to_plot %>%
          left_join(all_protein_names %>% dplyr::select(c(from_id, primary_gene_id)), 
                    by = c("biomarker" = "from_id")) %>%
          dplyr::select(-c(biomarker)) %>%
          dplyr::rename(c("biomarker" = "primary_gene_id"))
      }
      x_lab <- paste0(x_lab, "(", length(biomarkers), ")")
      
      biomarker_agg <- data_to_plot %>%
        group_by(biomarker) %>%
        summarize(med_expr = median(norm_expr)) %>%
        arrange(desc(med_expr))
      
      data_to_plot <- data_to_plot %>%
        mutate(Label = factor(Label, levels = rev(classes)),
               biomarker = factor(biomarker, biomarker_agg$biomarker)) 
      
      ggplot(data_to_plot, aes(x = biomarker, 
                               y = norm_expr,
                               fill = Label)) +
        geom_boxplot(size = 0.2, alpha = 0.5) +
        ggtitle(paste(sub("Vs", " Vs ", comparison), omics_type)) +
        xlab(x_lab) +
        ylab(y_lab) +
        guides(fill = guide_legend(title = "Condition")) +
        theme(axis.text.x = element_text(size=rel(1.2), angle = 90),
              axis.text.y = element_text(size=rel(1.2)),
              axis.title.x = element_text(size=rel(1.3)),
              axis.title.y = element_text(size=rel(1.3)),
              plot.title  = element_text(size=rel(1.5)))  
      
      file_name <- paste(comparison, omics_type, i, ".jpg",
                         sep = "_")
      file_path <- paste(boxplot_dir_path, file_name, sep = "/")
      ggplot2::ggsave(file_path, units = "cm", width = 30)
    }
    
    
    #plot dim red plot just for the best biomarker set chosen
    best_features_sub <- best_features %>%
      mutate(dataset_id = gsub(dataset_replace_str, "", dataset_id)) %>%
      filter(is_best == 1, dataset_id == comparison)
    biomarkers <- strsplit(best_features_sub[1, "biomarkers"], split = "|", fixed = TRUE)[[1]] 
    data <- data[, biomarkers]
    
    title <- paste(title, "best biomarkers")
    
    # 
    # data.cohort1 <- data.frame(data.cohort1)[, biomarkers]  #data.frame() replaces - in colnames to .
    # colnames(data.cohort2) <- gsub("-", ".", colnames(data.cohort2))

  }
  ################obtain best biomarkers end
  
  
  
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
    ggplot2::geom_point(ggplot2::aes(shape = output_labels$Label,
                                     fill = output_labels$Subgroup,
                                     colour = output_labels$data_cohort), size = 3, stroke = 1) +
    geom_text_repel(aes(label = text), size = 3, colour = "grey") +
    ggplot2::labs(title = title) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::scale_shape_manual(name = "Label (Point Shape)", values = c(21, 22)) +
    ggplot2::scale_colour_manual(name = "Data Cohort (Point Border)", values = c("black", "grey")) +
    ggplot2::guides(fill = guide_legend(override.aes = list(shape = 21, colour = NA))) +
    ggplot2::guides(colour = guide_legend(override.aes = list(shape = 1))) +
    labs(caption = paste(paste("Data dimension :", paste(dim(data), collapse = "x")), "\n",
                         group_counts_text),
         fill = "Subgroup (Point Fill)")
  
  if(!dir.exists(plot_dir_path)){
    dir.create(plot_dir_path, recursive = TRUE)
  }
  file_name <- paste(file_name_prefix, paste0(gsub(title, pattern = " |,", replacement = "-"), ".jpg"),
                     sep = "_")
  file_path <- paste(plot_dir_path, file_name, sep = "/")
  ggplot2::ggsave(file_path, units = "cm", width = 30)
}


