library(diann)
library(tidyverse)
library(readxl)
library(ggvenn)

source("scripts/R/utils.R")

df <- diann_load("Data/Protein/SBGN Plasma-EV Analysis/SBGN Plasma-EV DIA.tsv")
protein.groups <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], 
                               group.header="Protein.Group", 
                               id.header = "Precursor.Id", 
                               quantity.header = "Precursor.Normalised")

# protein.groups <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], 
#                                group.header="Protein.Group", 
#                                id.header = "Precursor.Id", 
#                                quantity.header = "PG.Normalised")


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

protein.groups.soi <- data.frame(protein.groups.soi)
protein.groups.soi <- protein.groups.soi %>%
  mutate(X = rownames(.), .before = "SB1")
row.names(protein.groups.soi) <- c()

write.csv(protein.groups.soi, "Data/Protein/formatted_data/new_cohort_samples_of_interest.csv",
          row.names = FALSE)

protein.groups.soi2 <- read.csv("Data/Protein/formatted_data/new_cohort_samples_of_interest.csv")
all.equal(protein.groups.soi, protein.groups.soi2)


strange_proteins <- protein.groups.soi %>%
  filter(grepl('.', x = X, fixed = TRUE)) %>%
  select(X)

protein.groups.soi %>%
  filter(X %in% sub(".1", "", strange_proteins$X[!strange_proteins$X %in% 'NA.']) ) %>%
  select(X)



#############
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
title_prefix <- "second_New cohort quantile normalization"
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



##############


#imputation


protein.groups.soi.filtered <- protein.groups.soi %>%
  filter(!grepl('.', x = X, fixed = TRUE)) %>%
  column_to_rownames("X")
protein.groups.soi.filtered <- log2(protein.groups.soi.filtered)


na_perc_df <- data.frame(matrix(nrow = 0, ncol = 2, dimnames = list(c(), c("perc", "protein_count"))))
for(filter_na_perc in c(25, 30, 40, 50, 60, 70, 75, 80, 90)){
  data <- t(protein.groups.soi.filtered)
  # filter_na_perc = 50
  data <- data.frame(data[, colSums(is.na(data)) < (filter_na_perc/100*nrow(data))])
  na_perc_df <- rbind(na_perc_df, data.frame(perc = filter_na_perc, protein_count = dim(data)[2]))
  print(paste(filter_na_perc, paste0(dim(data), collapse = "_")))  
}
write.csv(na_perc_df, "Data/Protein/formatted_data/filterna_validationcohort_diann.csv", row.names = FALSE)

new_cohort.imputed <- impute_data(t(protein.groups.soi.filtered), 
                              filter_na_perc = 50,
                              impute = TRUE)
dim(new_cohort.imputed)
new_cohort.imputed <- data.frame(new_cohort.imputed)
new_cohort.imputed <- t(new_cohort.imputed)
new_cohort.imputed <- data.frame(new_cohort.imputed)

ggvenn(list("Initial cohort" = rownames(initial_cohort), 
            "Validation Cohort filtered 50" = rownames(new_cohort.imputed)),
       fill_color = c("green", "orange"),
       stroke_size = 0.1,
       set_name_size = 5,
       text_size = 3)
ggsave(paste0("plots/proteins/common_proteins_bw_cohorts_50filtered.png"))



validation_metadata <- read.csv("Data/RNA_validation/metadata_glionet.csv")



data <- data.frame(t(new_cohort.imputed[, validation_metadata$sample_id]))
groups <- validation_metadata$category_old_name
title_prefix <- "New cohort imputed no normalization"
create_box_plot(data, groups, title_prefix)


data <- data.frame(new_cohort.imputed[, validation_metadata$sample_id])
norm_data <- preprocessCore::normalize.quantiles(as.matrix(data))
norm_data <- data.frame(norm_data, row.names = rownames(data))
colnames(norm_data) <- colnames(data)
data <- as.data.frame(t(as.matrix(norm_data)))
groups <- validation_metadata$category_old_name
title_prefix <- "New cohort imputed quantile normalization"
create_box_plot(data, groups, title_prefix)




#############

df <- diann_load("Data/Protein/Less strict/HBB Plasma-EV.tsv")
protein.groups <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], 
                               group.header="Protein.Group", 
                               id.header = "Precursor.Id", 
                               quantity.header = "Precursor.Normalised")

dim(protein.groups)


colnames(protein.groups) <- sub("R:\\PRJ-swath_2021\\HBB GBM-EV SWATH 2020 Study\\RAW SWATH Data - TripleTOF\\SWATH MS files\\SWATH data\\MzML format VRDprot\\",
                                "",
                                colnames(protein.groups), fixed = TRUE)
colnames(protein.groups) <- sapply(colnames(protein.groups), FUN = function(x){
  return(strsplit(x, split = "_", fixed = TRUE)[[1]][3])
}) 

colnames(protein.groups) <- gsub("rep", "", colnames(protein.groups))

phenotype <- read.table("Data/proteomic_phenotype.txt", header=TRUE, sep="\t")
phenotype <- phenotype %>%
  mutate(group = case_when(
    PREOPE_POSTOPE_TP_PREREC_REC_TP == "PREREC" ~ NA_character_,
    TRUE ~ PREOPE_POSTOPE_TP_PREREC_REC_TP
  ))
phenotype <- phenotype %>%
  filter(!is.na(group))

protein.groups[, phenotype$Sample]


phenotype$Sample[!phenotype$Sample %in% colnames(protein.groups)]

sum(phenotype$Sample %in% colnames(protein.groups))


phenotype <- phenotype %>%
  filter(Sample != "HB27")

rownames(protein.groups) <- sapply(strsplit(rownames(protein.groups), split = "|", fixed = TRUE), function(x){x[2]})
data <- data.frame(protein.groups[, phenotype$Sample])

ggvenn(list("Initial cohort original" = rownames(initial_cohort), 
            "Initial cohort diann" = rownames(protein.groups)),
       fill_color = c("green", "orange"),
       stroke_size = 0.1,
       set_name_size = 5,
       text_size = 3)
ggsave(paste0("plots/proteins/common_proteins_initial_cohort.png"))


# filter_na_perc = 50
data <- data.frame(t(data))
data <- data.frame(data[, colSums(is.na(data)) < (filter_na_perc/100*nrow(data))])

ggvenn(list("original" = rownames(initial_cohort), 
            "diann filtered 50" = colnames(data)),
       fill_color = c("green", "orange"),
       stroke_size = 0.1,
       set_name_size = 5,
       text_size = 3)
ggsave(paste0("plots/proteins/common_proteins_initial_cohort_diannfiltered.png"))



###########################


#new cohort from skyline

new_cohort <- read.csv("Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil.csv", 
                       row.names = 1)
initial_cohort <- read.csv("Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv", 
                           row.names = 1)

ggvenn(list("Initial cohort" = rownames(initial_cohort), 
            "Validation Cohort" = rownames(new_cohort)),
       fill_color = c("green", "orange"),
       stroke_size = 0.1,
       set_name_size = 5,
       text_size = 3)
ggsave(paste0("plots/proteins/common_proteins_bw_cohorts_skyline.png"))


validation_metadata <- read.csv("Data/RNA_validation/metadata_glionet.csv")

validation_metadata$sample_id
validation_metadata$sample_id %in% colnames(new_cohort)

qc_samples <- c("SB12_01", "SB12_02", "SB12_03", "SB12_04", "SB12_05", "SB12_06", "SB12_200ng", "SB12_300ng",
                "SB22", "SB22.02",
                "SB53", "SB53_2",
                "SB54", "SB54_2")

data <- data.frame(t(new_cohort[, qc_samples]))
groups <- rep("QC", length(qc_samples))
title_prefix <- "New cohort imputed QC samples"
create_box_plot(data, groups, title_prefix)

#sb12 1, 2, 3, 4, 5, 6 are similar
#use 1
colnames(new_cohort)[colnames(new_cohort) == "SB12_01"] = "SB12"

#use SB22.02
colnames(new_cohort)[colnames(new_cohort) == "SB22.02"] = "SBtobeused22"
colnames(new_cohort)[colnames(new_cohort) == "SB22"] = "SB22_dont_include"
colnames(new_cohort)[colnames(new_cohort) == "SBtobeused22"] = "SB22"

validation_metadata <- validation_metadata %>%
  filter(sample_id != "SB7")

data <- data.frame(t(new_cohort[, validation_metadata$sample_id]))
groups <- validation_metadata$category_old_name
title_prefix <- "New cohort imputed no normalization"
create_box_plot(data, groups, title_prefix,
                ylabel = "Expression across proteins (log2 intensity)")




data <- data.frame(new_cohort[, validation_metadata$sample_id])
norm_data <- preprocessCore::normalize.quantiles(as.matrix(data))
norm_data <- data.frame(norm_data, row.names = rownames(data))
colnames(norm_data) <- colnames(data)
data <- as.data.frame(t(as.matrix(norm_data)))
groups <- validation_metadata$category_old_name
title_prefix <- "New cohort imputed quantile normalization"
create_box_plot(data, groups, title_prefix,
                ylabel = "Expression across proteins (log2 intensity)")



#new cohort protein count
protein_data <- read.csv("Data/Protein/norm_output/norm__newcohort_processed_NA_FALSE.csv")

protein_data <- protein_data %>%
  arrange(SUBJECT_ORIGINAL) %>%
  mutate(SUBJECT_ORIGINAL = sub("SB12_01", "SB12", SUBJECT_ORIGINAL)) %>%
  mutate(SUBJECT_ORIGINAL = sub("SB22-02", "SBcorrect", SUBJECT_ORIGINAL)) %>%
  mutate(SUBJECT_ORIGINAL = sub("SB22", "SB-ignore-22", SUBJECT_ORIGINAL)) %>%
  mutate(SUBJECT_ORIGINAL = sub("SBcorrect", "SB22", SUBJECT_ORIGINAL))

formatted_data <-  protein_data %>%
  select(-c(GROUP_ORIGINAL)) %>%
  column_to_rownames("SUBJECT_ORIGINAL")

print(max(formatted_data, na.rm = TRUE))
print(min(formatted_data, na.rm = TRUE))

print(sum(is.na(formatted_data)))

protein_count <- rowSums(!is.na(formatted_data))

summary(protein_count)
protein_count <- data.frame(protein_count) %>%
  rownames_to_column("sample_id")

data_to_plot <- protein_count %>%
  inner_join(validation_metadata %>% select(sample_id, category_old_name))

data_to_plot <- data_to_plot %>%
  separate(sample_id, into = c(NA, "suf"), sep = "B", remove = FALSE, convert = TRUE) %>%
  arrange(category_old_name, suf) %>%
  select(-c(suf)) 

data_to_plot <- data_to_plot %>%
  mutate(sample_id = factor(sample_id, levels = unique(data_to_plot$sample_id)))

ggplot(data_to_plot) +
  geom_bar(aes(x = sample_id, y = protein_count, fill = category_old_name), stat = "identity") +
  xlab("Sample") +
  ylab("Protein Count") +
  labs(fill = "category") +
  ggtitle("Protein count validation cohort (skyline)") +
  scale_y_continuous(breaks = seq(0, max(data_to_plot$protein_count), by = 200)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/proteins/new_cohort_protein_count.png")

summary(data_to_plot$protein_count)



##initial cohort quantile norm
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

x.train <- data.frame(initial_cohort[, phenotype$Sample])
x.test <- data.frame(new_cohort[, validation_metadata$sample_id])


x.train.rank <- apply(x.train, 2, rank, ties.method="average")
x.train.sorted <- data.frame(apply(x.train, 2, sort))
x.train.mean <- apply(x.train.sorted, 1, mean)
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
x.train.norm <- apply(x.train.rank, 2, index_to_mean, data_mean = x.train.mean)
rownames(x.train.norm) <- rownames(x.train)
x.train <- x.train.norm

x.test.rank <- apply(x.test, 2, rank, ties.method="average")
#use params i.e. mean values of rows, from training data
x.test.norm <- apply(x.test.rank, 2, index_to_mean, data_mean = x.train.mean)
rownames(x.test.norm) <- rownames(x.test)
x.test <- x.test.norm

x.train <- as.data.frame(t(as.matrix(x.train)))
x.test <- as.data.frame(t(as.matrix(x.test))) 

groups <- phenotype$group  
title_prefix <- "Initial cohort quantile normalization new"
create_box_plot(x.train, groups, title_prefix,
                ylabel = "Expression across proteins (log2 intensity)")

groups <- validation_metadata$category_old_name
title_prefix <- "New cohort imputed quantile normalization with train param"
create_box_plot(x.test, groups, title_prefix,
                ylabel = "Expression across proteins (log2 intensity)")
