#PREOPE Vs POSTOPE_T proteomic 
#  and PREOPE Vs REC_TP transcriptomic showed no variation in the boxplots
#  shown in https://unsw.sharepoint.com/:p:/s/GBMPlasmaEV/EU1ARRr9SjNAqCtrI6ks__gBloWCQhelasTQI66hmO581w?e=Z4kNnN

#checking if those biomarkers obtained constant value due to imputation step

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

library(tidyverse)


normalize_data <- function(x, norm){
  #x : transcripts x samples
  #output x.norm : transcripts x samples
  if(norm == "norm_log_cpm_simple"){
    x.norm <- edgeR::cpm(x, log=TRUE)
    x.norm <- as.data.frame(t(as.matrix(x.norm)))    
    normparam <- caret::preProcess(x.norm) 
    x.norm <- predict(normparam, x.norm)
    x.norm <- as.data.frame(t(as.matrix(x.norm)))
  } else if(norm == "quantile"){
    x.norm <- preprocessCore::normalize.quantiles(as.matrix(x))
    x.norm <- data.frame(x.norm, row.names = rownames(x))
    colnames(x.norm) <- colnames(x)
  } else{
    print("norm method not supported")
    return
  }
  return (x.norm)
}


best_features_file_path = "Data/selected_features/best_features_with_add_col.csv"
dataset_id <- "GBMPlasmaEV_proteomic_impute50fil_quantile_PREOPEVsPOSTOPE_T"
phenotype_file_path <- "Data/proteomic_phenotype.txt"
plot_dir_path = "plots/FEMPipeline/selected_features_boxplot/issue"


data <- read.table("Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv", header=TRUE, sep=",", row.names=1, skip=0,
                   nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")

input_data <- read.csv(file = "Data/Protein/norm_output/norm_annotatedQ1-6_NA_FALSE.csv")
formatted_data <-  input_data %>%
  select(-c(GROUP_ORIGINAL)) %>%
  column_to_rownames("SUBJECT_ORIGINAL")
print(max(formatted_data, na.rm = TRUE))
print(min(formatted_data, na.rm = TRUE))
print(sum(is.na(formatted_data)))
formatted_data <- t(formatted_data)
input_data <- data.frame(formatted_data)

sum(is.na(data))
data["P0DOX6",]

sum(is.na(input_data))
input_data["P0DOX6",]


data_to_plot <- data.frame("original" = t(input_data["Q5RKV6",]),
                           "imputed" = t(data["Q5RKV6",]))
colnames(data_to_plot) <- c("original", "imputed")

data_to_plot <- data_to_plot %>%
  pivot_longer(c("original", "imputed"), names_to = "data_type", values_to = "expr_value")
ggplot(data_to_plot, aes(x = data_type, 
                         y = expr_value,
                         fill = data_type)) +
  geom_boxplot(size = 0.2, alpha = 0.5) +
  ggtitle("Expression of Q5RKV6 in all samples") +
  theme(axis.text.x = element_text(size=rel(1.2), angle = 45, hjust = 1),
        axis.text.y = element_text(size=rel(1.2)),
        axis.title.x = element_text(size=rel(1.5)),
        axis.title.y = element_text(size=rel(1.5)),
        plot.title  = element_text(size=rel(1.5)))
ggsave("plots/FEMPipeline/selected_features_boxplot/issue/Q5RKV6_all_samples.png")



data_to_plot <- data.frame("original" = t(input_data["P0DOX6",]),
                           "imputed" = t(data["P0DOX6",]))
colnames(data_to_plot) <- c("original", "imputed")
boxplot(data_to_plot)
data_to_plot <- data_to_plot %>%
  pivot_longer(c("original", "imputed"), names_to = "data_type", values_to = "expr_value")
ggplot(data_to_plot, aes(x = data_type, 
                         y = expr_value,
                         fill = data_type)) +
  geom_boxplot(size = 0.2, alpha = 0.5) +
  ggtitle("Expression of P0DOX6 in PREOPE, POSTOPE_T samples") +
  theme(axis.text.x = element_text(size=rel(1.2), angle = 45, hjust = 1),
        axis.text.y = element_text(size=rel(1.2)),
        axis.title.x = element_text(size=rel(1.5)),
        axis.title.y = element_text(size=rel(1.5)),
        plot.title  = element_text(size=rel(1.5)))
ggsave("plots/FEMPipeline/selected_features_boxplot/issue/p0dox6_PREOPEVsPOSTOPE_T_samples.png")

classification_criteria <- strsplit(dataset_id, split = split_str)[[1]]
classification_criteria <- classification_criteria[length(classification_criteria)]
extracted_samples <- phenotype[!is.na(phenotype[classification_criteria]), ]
data <- data %>% select(extracted_samples$Sample) 
input_data <- input_data %>% select(extracted_samples$Sample) 
data_to_plot <- data.frame("original" = t(input_data["P0DOX6",]),
                           "imputed" = t(data["P0DOX6",]))
colnames(data_to_plot) <- c("original", "imputed")
boxplot(data_to_plot)
data_to_plot <- data_to_plot %>%
  pivot_longer(c("original", "imputed"), names_to = "data_type", values_to = "expr_value")
ggplot(data_to_plot, aes(x = data_type, 
                         y = expr_value,
                         fill = data_type)) +
  geom_boxplot(size = 0.2, alpha = 0.5) +
  ggtitle("Expression of P0DOX6 in all samples") +
  theme(axis.text.x = element_text(size=rel(1.2), angle = 45, hjust = 1),
        axis.text.y = element_text(size=rel(1.2)),
        axis.title.x = element_text(size=rel(1.5)),
        axis.title.y = element_text(size=rel(1.5)),
        plot.title  = element_text(size=rel(1.5)))
ggsave("plots/FEMPipeline/selected_features_boxplot/issue/p0dox6_all_samples.png")

