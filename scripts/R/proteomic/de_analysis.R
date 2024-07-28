#DE analysis for all proteomic data - PREOPE, MET, HC - 2024 Jan 

library(edgeR)
library(tidyverse)
library(sva)

data <- read.csv("Data/Protein/formatted_data/PREOPE_MET_HC_data.csv", row.names = 1)
phenotype <- read.table("Data/proteomic_phenotype_PREOPE_MET_HC_withaddicolumn.txt", header=TRUE, sep="\t")

data <- data[, phenotype$Sample]
all.equal(phenotype$Sample, colnames(data))

#quantile norm

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
data <- as.data.frame(data.norm)

# plotMDS(data, labels = phenotype$PREOPE_MET_HC)

summary(factor(paste(phenotype$data_cohort, phenotype$PREOPE_MET_HC)))
# 
# initial HC       initial MET    initial PREOPE validation PREOPE 
# 21                21                11                15

model_matrix <- model.matrix(~ 0 + phenotype$PREOPE_MET_HC)
colnames(model_matrix) <- sub("phenotype$PREOPE_MET_HC", "", colnames(model_matrix), fixed = TRUE)


contr_matrix <- makeContrasts(contrasts = "PREOPE - MET",
                              levels = colnames(model_matrix))
fit <- lmFit(data, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("protein")
result <- top.table %>%
  dplyr::select(protein, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = protein, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)


# 'B1AJZ1' %in% rownames(data)
# 'Q674X7' %in% rownames(data)
# 
# 'Q6ZMK1' %in% rownames(data)
# 'P0DTL5' %in% rownames(data)
# 
# 'P13746' %in% rownames(data)
# 'P04439' %in% rownames(data)

plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs MET",
                         output_dir_path = "DE_results_2024/proteomics/1_only_condition/p/",
                         plot_file_name = "PREOPEVsMET.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         molecule_names_file_path = "Data/Protein/formatted_data/all_protein_names.csv",
                         molecule_names_file_columns = c(1, 3, 4),
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs MET",
                         output_dir_path = "DE_results_2024/proteomics/1_only_condition/padj/",
                         plot_file_name = "PREOPEVsMET.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         molecule_names_file_path = "Data/Protein/formatted_data/all_protein_names.csv",
                         molecule_names_file_columns = c(1, 3, 4),
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "PREOPE - HC",
                              levels = colnames(model_matrix))
fit <- lmFit(data, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("protein")
result <- top.table %>%
  dplyr::select(protein, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = protein, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)
plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs HC",
                         output_dir_path = "DE_results_2024/proteomics/1_only_condition/p/",
                         plot_file_name = "PREOPEVsHC.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         molecule_names_file_path = "Data/Protein/formatted_data/all_protein_names.csv",
                         molecule_names_file_columns = c(1, 3, 4),
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs HC",
                         output_dir_path = "DE_results_2024/proteomics/1_only_condition/padj/",
                         plot_file_name = "PREOPEVsHC.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         molecule_names_file_path = "Data/Protein/formatted_data/all_protein_names.csv",
                         molecule_names_file_columns = c(1, 3, 4),
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "MET - HC",
                              levels = colnames(model_matrix))
fit <- lmFit(data, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("protein")
result <- top.table %>%
  dplyr::select(protein, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = protein, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)
plot_volcano_and_save_DE(result, plot_title = "MET Vs HC",
                         output_dir_path = "DE_results_2024/proteomics/1_only_condition/p/",
                         plot_file_name = "METVsHC.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         molecule_names_file_path = "Data/Protein/formatted_data/all_protein_names.csv",
                         molecule_names_file_columns = c(1, 3, 4),
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "MET Vs HC",
                         output_dir_path = "DE_results_2024/proteomics/1_only_condition/padj/",
                         plot_file_name = "METVsHC.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         molecule_names_file_path = "Data/Protein/formatted_data/all_protein_names.csv",
                         molecule_names_file_columns = c(1, 3, 4),
                         plot_width_cm = 25)



##############################

#now with batch in limma

phenotype <- phenotype %>%
  mutate(data_cohort = case_when(data_cohort == "initial" ~ "Cohort1",
                                 TRUE ~ 'Cohort2'))

model_matrix <- model.matrix(~ 0 + phenotype$PREOPE_MET_HC + phenotype$data_cohort)
colnames(model_matrix) <- sub("phenotype$PREOPE_MET_HC", "", colnames(model_matrix), fixed = TRUE)
colnames(model_matrix) <- sub("phenotype$data_cohort", "", colnames(model_matrix), fixed = TRUE)

contr_matrix <- makeContrasts(contrasts = "PREOPE - MET",
                              levels = colnames(model_matrix))
fit <- lmFit(data, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("protein")
result <- top.table %>%
  dplyr::select(protein, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = protein, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs MET",
                         output_dir_path = "DE_results_2024/proteomics/2_batch_in_limma/p/",
                         plot_file_name = "PREOPEVsMET.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         molecule_names_file_path = "Data/Protein/formatted_data/all_protein_names.csv",
                         molecule_names_file_columns = c(1, 3, 4),
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs MET",
                         output_dir_path = "DE_results_2024/proteomics/2_batch_in_limma/padj/",
                         plot_file_name = "PREOPEVsMET.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         molecule_names_file_path = "Data/Protein/formatted_data/all_protein_names.csv",
                         molecule_names_file_columns = c(1, 3, 4),
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "PREOPE - HC",
                              levels = colnames(model_matrix))
fit <- lmFit(data, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("protein")
result <- top.table %>%
  dplyr::select(protein, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = protein, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs HC",
                         output_dir_path = "DE_results_2024/proteomics/2_batch_in_limma/p/",
                         plot_file_name = "PREOPEVsHC.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         molecule_names_file_path = "Data/Protein/formatted_data/all_protein_names.csv",
                         molecule_names_file_columns = c(1, 3, 4),
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs HC",
                         output_dir_path = "DE_results_2024/proteomics/2_batch_in_limma/padj/",
                         plot_file_name = "PREOPEVsHC.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         molecule_names_file_path = "Data/Protein/formatted_data/all_protein_names.csv",
                         molecule_names_file_columns = c(1, 3, 4),
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "MET - HC",
                              levels = colnames(model_matrix))
fit <- lmFit(data, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("protein")
result <- top.table %>%
  dplyr::select(protein, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = protein, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "MET Vs HC",
                         output_dir_path = "DE_results_2024/proteomics/2_batch_in_limma/p/",
                         plot_file_name = "METVsHC.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         molecule_names_file_path = "Data/Protein/formatted_data/all_protein_names.csv",
                         molecule_names_file_columns = c(1, 3, 4),
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "MET Vs HC",
                         output_dir_path = "DE_results_2024/proteomics/2_batch_in_limma/padj/",
                         plot_file_name = "METVsHC.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         molecule_names_file_path = "Data/Protein/formatted_data/all_protein_names.csv",
                         molecule_names_file_columns = c(1, 3, 4),
                         plot_width_cm = 25)


##############################

#now with combat batch corrected and no batch specified in limma

data.combat <- ComBat(dat = data, batch = phenotype$data_cohort)
data.combat <- as.data.frame(data.combat)

model_matrix <- model.matrix(~ 0 + phenotype$PREOPE_MET_HC)
colnames(model_matrix) <- sub("phenotype$PREOPE_MET_HC", "", colnames(model_matrix), fixed = TRUE)

contr_matrix <- makeContrasts(contrasts = "PREOPE - MET",
                              levels = colnames(model_matrix))
fit <- lmFit(data.combat, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("protein")
result <- top.table %>%
  dplyr::select(protein, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = protein, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "GBM Vs MET",
                         output_dir_path = "DE_results_2024/proteomics/3_combat_corrected/p/",
                         plot_file_name = "GBMVsMET.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         molecule_names_file_path = "Data/Protein/formatted_data/all_protein_names.csv",
                         molecule_names_file_columns = c(1, 3, 4),
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "GBM Vs MET",
                         output_dir_path = "DE_results_2024/proteomics/3_combat_corrected/padj/",
                         plot_file_name = "GBMVsMET.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         molecule_names_file_path = "Data/Protein/formatted_data/all_protein_names.csv",
                         molecule_names_file_columns = c(1, 3, 4),
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "PREOPE - HC",
                              levels = colnames(model_matrix))
fit <- lmFit(data.combat, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("protein")
result <- top.table %>%
  dplyr::select(protein, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = protein, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)
plot_volcano_and_save_DE(result, plot_title = "GBM Vs HC",
                         output_dir_path = "DE_results_2024/proteomics/3_combat_corrected/p/",
                         plot_file_name = "GBMVsHC.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         molecule_names_file_path = "Data/Protein/formatted_data/all_protein_names.csv",
                         molecule_names_file_columns = c(1, 3, 4),
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "GBM Vs HC",
                         output_dir_path = "DE_results_2024/proteomics/3_combat_corrected/padj/",
                         plot_file_name = "GBMVsHC.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         molecule_names_file_path = "Data/Protein/formatted_data/all_protein_names.csv",
                         molecule_names_file_columns = c(1, 3, 4),
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "MET - HC",
                              levels = colnames(model_matrix))
fit <- lmFit(data.combat, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("protein")
result <- top.table %>%
  dplyr::select(protein, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = protein, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)
plot_volcano_and_save_DE(result, plot_title = "MET Vs HC",
                         output_dir_path = "DE_results_2024/proteomics/3_combat_corrected/p/",
                         plot_file_name = "METVsHC.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         molecule_names_file_path = "Data/Protein/formatted_data/all_protein_names.csv",
                         molecule_names_file_columns = c(1, 3, 4),
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "MET Vs HC",
                         output_dir_path = "DE_results_2024/proteomics/3_combat_corrected/padj/",
                         plot_file_name = "METVsHC.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         molecule_names_file_path = "Data/Protein/formatted_data/all_protein_names.csv",
                         molecule_names_file_columns = c(1, 3, 4),
                         plot_width_cm = 25)