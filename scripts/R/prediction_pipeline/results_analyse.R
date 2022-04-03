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
