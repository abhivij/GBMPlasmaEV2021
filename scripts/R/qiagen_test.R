library(tidyverse)
library(readxl)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

old_output <- read_excel("Data/RNA/158629.all_samples.summary.xlsx", 
                         sheet = "summary")
new_output <- read_excel("Data/qiagen_results_test/GBMPlasmaEV_initial_cohort_RNAData/219155.all_samples.summary.xlsx",
                        sheet = "summary")

samples_old <- sort(colnames(old_output)[-1])
samples_new <- sort(colnames(new_output)[-1])

all.equal(samples_old, samples_new)
#TRUE

old_output.sort <- old_output[, c(colnames(old_output)[1], samples_old)]
new_output.sort <- new_output[, c(colnames(old_output)[1], samples_old)]

colSums(old_output.sort[-1,-1])

old_output.sort.ex3 <- old_output.sort[-c(5, 12, 13), ]
new_output.sort.ex3 <- new_output.sort[-c(5, 12, 13), ]


all.equal(old_output.sort.ex3, new_output.sort.ex3)



old_output.sort.ex3.t <- data.frame(t(old_output.sort.ex3 %>%
                             column_to_rownames("read set"))) 

old_output.sort.ex3.t <- old_output.sort.ex3.t %>%
  mutate(short_reads_perc = too_short_reads * 100 / total_reads)

summary(old_output.sort.ex3.t$short_reads_perc)



old_output.sort.t <- data.frame(t(old_output.sort %>%
                                    column_to_rownames("read set")))
old_output.sort.t <- old_output.sort.t %>%
  mutate(short_reads_perc = too_short_reads * 100 / total_reads) %>%
  mutate(mir_reads_perc = miRNA_Reads * 100 / total_reads)
summary(old_output.sort.t$short_reads_perc)
summary(old_output.sort.t$mir_reads_perc)



new_output.sort.t <- data.frame(t(new_output.sort %>%
                                    column_to_rownames("read set")))
new_output.sort.t <- new_output.sort.t %>%
  mutate(short_reads_perc = too_short_reads * 100 / total_reads) %>%
  mutate(mir_reads_perc = miRNA_Reads * 100 / total_reads)
summary(new_output.sort.t$short_reads_perc)
summary(new_output.sort.t$mir_reads_perc)
