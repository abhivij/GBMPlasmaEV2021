library(tidyverse)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)


x <- read.table("Data/RNA/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")

metadata <- read.csv("Data/transcriptomic_sample_metadata.csv")

metadata <- metadata %>%
  select(c(1,2))
colnames(metadata) <- c("sample", "condition")
metadata <- metadata %>%
  filter(condition != "OUT")

x <- x[, metadata$sample]
#   x format (transcripts x samples)

x.norm <- edgeR::cpm(x, log=TRUE)
x.norm <- as.data.frame(t(as.matrix(x.norm)))

#normalizing the data
normparam <- caret::preProcess(x.norm)
x.norm <- predict(normparam, x.norm)
#after this step, sum of columns (i.e. transcripts) will be 0



data <- x.norm %>%
  rownames_to_column("sample") 


data <- metadata %>%
  inner_join(data) %>%
  arrange(condition) %>%
  select(-c(sample))


result <- data %>%
  group_by(condition) %>%
  summarise(across(where(is.numeric), mean))

result_var <- result %>%
  summarise(across(where(is.numeric), var))
result_var <- as.data.frame(t(result_var))
colnames(result_var) <- "variance_across_conditions"

result_var <- result_var %>%
  rownames_to_column("transcripts") %>%
  arrange(variance_across_conditions) 
var_threshold <- quantile(result_var[, 2], )["25%"]

filtered_result_var <- result_var %>%
  filter(variance_across_conditions <= var_threshold)


result_mean <- result %>%
  summarise(across(where(is.numeric), mean))
result_mean <- as.data.frame(t(result_mean))
colnames(result_mean) <- "mean_across_conditions"
result_mean <- result_mean %>%
  rownames_to_column("transcripts") %>%
  arrange(desc(mean_across_conditions)) 


var_and_mean <- filtered_result_var %>%
  inner_join(result_mean)
var_and_mean <- var_and_mean %>%
  arrange(desc(mean_across_conditions))


best_ones <- var_and_mean[1:5, ]

write.csv(best_ones, "candidate_housekeeping_genes.csv", row.names = FALSE)



###############




logcpm <- edgeR::cpm(x, log=TRUE)
logcpm <- as.data.frame(t(as.matrix(logcpm)))

logcpm <- logcpm[,best_ones$transcripts]



data_to_plot <- data.matrix(logcpm)
Heatmap(data_to_plot, name = "Expression logCPM",
             rect_gp = gpar(col = "black", lwd = 1),
             border = TRUE,
             row_names_side = "left", 
             show_column_dend = FALSE,
             show_row_dend = FALSE,
             show_row_names = FALSE,
             row_split = metadata$condition,
             left_annotation = HeatmapAnnotation(
               "condition" = metadata$condition,
               which = "row",
               show_annotation_name = FALSE,
               border = TRUE,
               gp = gpar(col = "black", lwd = 1),
               gap = unit(1, units = "mm")
             ))

