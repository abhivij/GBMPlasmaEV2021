library(diann)
library(tidyverse)
library(readxl)
library(ggvenn)

df <- diann_load("Data/Protein/SBGN Plasma-EV Analysis/SBGN Plasma-EV DIA.tsv")
protein.groups <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], 
                               group.header="Protein.Group", 
                               id.header = "Precursor.Id", 
                               quantity.header = "Precursor.Normalised")


sample_info <- read_xlsx("Data/Protein/SBGN Plasma-EV Analysis/Samples Input to DIA-NN_modified.xlsx")
colnames(sample_info) <- c("file_name", "sample_name", "remove", "comment")

sample_info <- sample_info %>%
  filter(is.na(remove))
sample_info <- sample_info %>%
  mutate(sample_name = sub("-", "", sample_name, fixed = TRUE))

colnames(protein.groups) <- sub("R:\\PRJ-DIA_SWATH_2022\\DIA Eclipse RAW data\\",
                                "",
                                colnames(protein.groups), fixed = TRUE)
colnames(protein.groups)

#samples of interest
protein.groups.soi <- protein.groups[, sample_info$file_name]
dim(protein.groups.soi)
# [1] 6764   55

colnames(protein.groups.soi) <- sample_info$sample_name
rownames(protein.groups.soi) <- sub("_HUMAN", "", rownames(protein.groups.soi), fixed = TRUE)
rownames(protein.groups.soi) <- sapply(strsplit(rownames(protein.groups.soi), split = "|", fixed = TRUE), function(x){x[2]})

protein.groups.soi[is.na(protein.groups.soi)] <- 0


initial_cohort <- read.csv("Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv", 
                           row.names = 1)


ggvenn(list("Initial cohort" = rownames(initial_cohort), 
            "Validation Cohort" = rownames(protein.groups.soi)),
       fill_color = c("green", "orange"),
       stroke_size = 0.1,
       set_name_size = 5,
       text_size = 3)
ggsave(paste0("plots/proteins/common_proteins_bw_cohorts.png"))

protein.groups.soi <- log2(protein.groups.soi)
protein.groups.soi[protein.groups.soi == -Inf] <- 0



validation_metadata <- read.csv("Data/RNA_validation/metadata_glionet.csv")



data <- data.frame(t(protein.groups.soi[, validation_metadata$sample_id]))
groups <- validation_metadata$category_old_name
title_prefix <- "New cohort no normalization"
create_box_plot(data, groups, title_prefix)


data <- data.frame(protein.groups.soi[, validation_metadata$sample_id])
norm_data <- preprocessCore::normalize.quantiles(as.matrix(data))
norm_data <- data.frame(norm_data, row.names = rownames(data))
colnames(norm_data) <- colnames(data)
data <- as.data.frame(t(as.matrix(norm_data)))
groups <- validation_metadata$category_old_name
title_prefix <- "New cohort quantile normalization"
create_box_plot(data, groups, title_prefix)



initial_cohort <- read.csv("Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv", 
                           row.names = 1)
phenotype <- read.table("Data/proteomic_phenotype.txt", header=TRUE, sep="\t")
phenotype <- phenotype %>%
  mutate(group = case_when(
    PREOPE_POSTOPE_TP_PREREC_REC_TP == "PREREC" ~ NA_character_,
    TRUE ~ PREOPE_POSTOPE_TP_PREREC_REC_TP
    ))
phenotype <- phenotype %>%
  filter(!is.na(group))

data <- data.frame(t(initial_cohort[, phenotype$Sample]))
groups <- phenotype$group  
title_prefix <- "Initial cohort no normalization"
create_box_plot(data, groups, title_prefix,
                ylabel = "Expression across proteins (log2 intensity)")





data <- data.frame(initial_cohort[, phenotype$Sample])
norm_data <- preprocessCore::normalize.quantiles(as.matrix(data))
norm_data <- data.frame(norm_data, row.names = rownames(data))
colnames(norm_data) <- colnames(data)
data <- as.data.frame(t(as.matrix(norm_data)))
groups <- phenotype$group  
title_prefix <- "Initial cohort quantile normalization"
create_box_plot(data, groups, title_prefix,
                ylabel = "Expression across proteins (log2 intensity)")


create_box_plot <- function(data, groups,
                            title_prefix = "",
                            ylabel = "Expression across proteins (log2 LFQ)"){
  data_to_plot <- cbind(data, "label" = groups) %>%
    rownames_to_column("sample_name")
  data_to_plot <- data_to_plot %>%
    pivot_longer(!c(sample_name, label), names_to = "molecules") %>%
    arrange(label)
  data_to_plot <- data_to_plot %>%
    mutate(sample_name = factor(sample_name, levels = unique(data_to_plot$sample_name)))
  
  ggplot(data_to_plot, aes(x = sample_name, y = value)) +
    geom_boxplot(aes(fill = label)) +
    xlab("Sample Name") +
    ylab(ylabel) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(title_prefix)
  
  title <- paste(title_prefix, "boxplot")
  
  dir_path <- "plots/proteins/"
  file_name <- paste0(gsub(title, pattern = " ", replacement = "-"), ".png")
  file_path <- paste(dir_path, file_name, sep = "/")
  ggplot2::ggsave(file_path, units = "cm", width = 20)
  
}



plot_data <- function(comparison, classes, 
                      best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                      use_best_transcripts = TRUE,
                      perform_filter = TRUE, norm = "norm_log_tmm", use_train_param = TRUE){
  
  data <- read.table("data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                     nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
  phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t")
  
  output_labels.train <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes, country == "AU") %>%
    dplyr::select(Sample, Label, age, age_group, sex, FEV1) %>%
    dplyr::mutate(Label = factor(Label), age = as.numeric(age), 
                  age_group = factor(age_group), sex = factor(sex),
                  FEV1 = as.numeric(FEV1)) %>%
    arrange(Label, Sample)
  data.train <- data[, output_labels.train$Sample]
  print(summary(output_labels.train))
  
  output_labels.test <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes, country == "DK") %>%
    dplyr::select(Sample, Label, age, age_group, sex, FEV1) %>%
    dplyr::mutate(Label = factor(Label), age = as.numeric(age), 
                  age_group = factor(age_group), sex = factor(sex),
                  FEV1 = as.numeric(FEV1)) %>%    
    arrange(Label, Sample)
  data.test <- data[, output_labels.test$Sample]
  print(summary(output_labels.test))
  
  #currently data.train, data.test format : (transcripts x samples)
  
  if(perform_filter){
    keep <- edgeR::filterByExpr(data.train, group = output_labels.train$Label)
    data.train <- data.train[keep, ]
    if(use_train_param){
      data.test <- data.test[keep, ]  
    } else{
      keep <- edgeR::filterByExpr(data.test, group = output_labels.test$Label)  
      data.test <- data.test[keep, ] 
    }
  }
  
  
  if(norm == "norm_log_tmm"){
    #calculating norm log tmm
    dge <- edgeR::DGEList(counts = data.train, group = output_labels.train$Label)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    tmm <- edgeR::cpm(dge, log = TRUE)
    data.train <- tmm
    
    dge <- edgeR::DGEList(counts = data.test, group = output_labels.test$Label)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    tmm <- edgeR::cpm(dge, log = TRUE)
    data.test <- tmm
    
    data.train <- as.data.frame(t(as.matrix(data.train)))
    data.test <- as.data.frame(t(as.matrix(data.test)))  
    
    #normalizing the data
    normparam <- caret::preProcess(data.train) 
    data.train <- predict(normparam, data.train)
    
    if(use_train_param){
      data.test <- predict(normparam, data.test) #normalizing test data using params from train data       
    } else{
      normparam <- caret::preProcess(data.test)
      data.test <- predict(normparam, data.test) #normalizing test data using params from train data 
    }
  }
  
  #now data.train, data.test format : (samples x transcripts)
  
  #get best biomarkers only
  if(use_best_transcripts){
    best_features <- read.csv(best_features_file_path)  
    best_features_sub <- best_features %>%
      mutate(dataset_id = gsub("CF_EV_AU_zlogtmm_", "", dataset_id)) %>%
      filter(is_best == 1, dataset_id == comparison)
    
    biomarkers <- strsplit(best_features_sub$biomarkers, split = "|", fixed = TRUE)[[1]]  
    
    features_with_slash <- colnames(data.train)[grepl("/", colnames(data.train), fixed = TRUE)] 
    for(f in features_with_slash){
      f_replaced <- gsub("/|-", ".", f) 
      if(f_replaced %in% biomarkers){
        biomarkers[biomarkers == f_replaced] = f
      }
    }
    biomarkers <- gsub(".", "-", biomarkers, fixed = TRUE)
    
    
    data.train <- data.train[, biomarkers]
    data.test <- data.test[, biomarkers]    
  }
  
  if(perform_filter && norm == "norm_log_tmm"){
    pp = TRUE
  } else{
    pp = FALSE
  }
  
  # create_dim_red_plots(data = data.train, groups = output_labels.train$Label, 
  #                      title_prefix = paste(comparison, 
  #                                           "best_transcripts", use_best_transcripts,
  #                                           "pp", pp, "train_params", use_train_param,
  #                                           "train"), 
  #                      dim_red = "UMAP")
  # create_dim_red_plots(data = data.test, groups = output_labels.test$Label, 
  #                      title_prefix = paste(comparison, 
  #                                           "best_transcripts", use_best_transcripts,
  #                                           "pp", pp, "train_params", use_train_param,
  #                                           "test"), 
  #                      dim_red = "UMAP")  
  # 
  # 
  # create_dim_red_plots(data = data.train, groups = output_labels.train$Label, 
  #                      title_prefix = paste(comparison, 
  #                                           "best_transcripts", use_best_transcripts,
  #                                           "pp", pp, "train_params", use_train_param,
  #                                           "train"), 
  #                      dim_red = "tSNE")
  # create_dim_red_plots(data = data.test, groups = output_labels.test$Label, 
  #                      title_prefix = paste(comparison, 
  #                                           "best_transcripts", use_best_transcripts,
  #                                           "pp", pp, "train_params", use_train_param,
  #                                           "test"), 
  #                      dim_red = "tSNE")  
  # 
  
  create_box_plot(data = data.train, groups = output_labels.train$Label, 
                  title_prefix = paste(comparison, 
                                       "best_transcripts", use_best_transcripts,
                                       "pp", pp, "train_params", use_train_param,
                                       "train"))
  create_box_plot(data = data.test, groups = output_labels.test$Label, 
                  title_prefix = paste(comparison, 
                                       "best_transcripts", use_best_transcripts,
                                       "pp", pp, "train_params", use_train_param,
                                       "test"))
  
}
