library(tidyverse)

transcriptomics <- read.csv("Data/prediction_result/transcriptomics.csv")

proteomics <- read.csv("Data/prediction_result/proteomics.csv")

t.t <- transcriptomics %>%
  filter(TestDataClassName != "PREREC")
mean(t.t$TestDataClassId == t.t$prediction_with_cutoff_0.37)
t.v <- transcriptomics %>%
  filter(TestDataClassName == "PREREC") 
mean(t.v$TestDataClassId == t.v$prediction_with_cutoff_0.37)

p.t <- proteomics %>%
  filter(TestDataClassName != "PREREC")
mean(p.t$TestDataExpectedClassName == p.t$Prediction)
p.v <- proteomics %>%
  filter(TestDataClassName == "PREREC") 
mean(p.v$TestDataExpectedClassName == p.v$Prediction)


result_file_name = "Data/prediction_result/transcriptomics_POSTOPE_TPVsREC_TP_with_new_validation_data.csv"
comparison = "POSTOPE_TPVsREC_TP"

create_result_heatmap <- function(result_file_name, comparison){
  result_df <- read.csv(result_file_name) %>%
    select(-c(3))
  cutoff <- strsplit(colnames(result_df)[4], split = "cutoff_")[[1]][2]
  
  colnames(result_df) <- c("sample", "true_label", "pred_prob", "pred", "type")

  
  train_result <- result_df %>%
    filter(type == "train")
  sum(train_result$true_label == train_result$pred)
  mean(train_result$true_label == train_result$pred)
  
  caret::confusionMatrix(data = factor(train_result$pred), 
                         reference = factor(train_result$true_label))
  
  
  test_result <- result_df %>%
    filter(type == "test", !is.na(true_label))
  sum(test_result$true_label == test_result$pred)
  mean(test_result$true_label == test_result$pred)
  
  caret::confusionMatrix(data = factor(test_result$pred), 
                         reference = factor(test_result$true_label))
  
}