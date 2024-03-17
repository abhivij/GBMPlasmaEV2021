#DE analysis for all transcriptomic data - PREOPE, MET, HC - 2024 Jan 

library(edgeR)
library(tidyverse)
library(sva)

data <- read.csv("Data/RNA_all/newquant_Nov2023_umi_counts_PREOPE_MET_HC_filter90.csv", row.names = 1)
phenotype <- read.table("Data/transcriptomic_phenotype_PREOPE_MET_HC_withaddicolumn.txt", header=TRUE, sep="\t")

data <- data[, phenotype$Sample]
all.equal(phenotype$Sample, colnames(data))
dim(data)
# [1] 907  67

keep <- edgeR::filterByExpr(data, group = phenotype$PREOPE_MET_HC)
data <- data[keep, ]

dim(data)
# [1] 207  67

data <- edgeR::cpm(data, log=TRUE)
data <- as.data.frame(data)

summary(factor(paste(phenotype$data_cohort, phenotype$PREOPE_MET_HC)))
# initial HC       initial MET    initial PREOPE validation PREOPE 
# 21                21                10                15 

model_matrix <- model.matrix(~ 0 + phenotype$PREOPE_MET_HC)
colnames(model_matrix) <- sub("phenotype$PREOPE_MET_HC", "", colnames(model_matrix), fixed = TRUE)

# v <- voom(data, model_matrix, plot = TRUE)
v <- data

contr_matrix <- makeContrasts(contrasts = "PREOPE - MET",
                              levels = colnames(model_matrix))
fit <- lmFit(v, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("tra")
result <- top.table %>%
  dplyr::select(tra, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = tra, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs MET",
                         output_dir_path = "DE_results_2024/transcriptomics/1_only_condition/p/",
                         plot_file_name = "PREOPEVsMET.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs MET",
                         output_dir_path = "DE_results_2024/transcriptomics/1_only_condition/padj/",
                         plot_file_name = "PREOPEVsMET.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "PREOPE - HC",
                              levels = colnames(model_matrix))
fit <- lmFit(v, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("tra")
result <- top.table %>%
  dplyr::select(tra, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = tra, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)
plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs HC",
                         output_dir_path = "DE_results_2024/transcriptomics/1_only_condition/p/",
                         plot_file_name = "PREOPEVsHC.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs HC",
                         output_dir_path = "DE_results_2024/transcriptomics/1_only_condition/padj/",
                         plot_file_name = "PREOPEVsHC.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "MET - HC",
                              levels = colnames(model_matrix))
fit <- lmFit(v, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("tra")
result <- top.table %>%
  dplyr::select(tra, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = tra, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)
plot_volcano_and_save_DE(result, plot_title = "MET Vs HC",
                         output_dir_path = "DE_results_2024/transcriptomics/1_only_condition/p/",
                         plot_file_name = "METVsHC.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "MET Vs HC",
                         output_dir_path = "DE_results_2024/transcriptomics/1_only_condition/padj/",
                         plot_file_name = "METVsHC.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         plot_width_cm = 25)

##############################

#now with batch in limma

phenotype <- phenotype %>%
  mutate(data_cohort = case_when(data_cohort == "initial" ~ "Cohort1",
                                 TRUE ~ 'Cohort2'))

model_matrix <- model.matrix(~ 0 + phenotype$PREOPE_MET_HC + phenotype$data_cohort)
colnames(model_matrix) <- sub("phenotype$PREOPE_MET_HC", "", colnames(model_matrix), fixed = TRUE)
colnames(model_matrix) <- sub("phenotype$data_cohort", "", colnames(model_matrix), fixed = TRUE)

# v <- voom(data, model_matrix, plot = TRUE)
v <- data

contr_matrix <- makeContrasts(contrasts = "PREOPE - MET",
                              levels = colnames(model_matrix))
fit <- lmFit(v, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("tra")
result <- top.table %>%
  dplyr::select(tra, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = tra, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs MET",
                         output_dir_path = "DE_results_2024/transcriptomics/2_batch_in_limma/p/",
                         plot_file_name = "PREOPEVsMET.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs MET",
                         output_dir_path = "DE_results_2024/transcriptomics/2_batch_in_limma/padj/",
                         plot_file_name = "PREOPEVsMET.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "PREOPE - HC",
                              levels = colnames(model_matrix))
fit <- lmFit(v, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("tra")
result <- top.table %>%
  dplyr::select(tra, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = tra, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)
plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs HC",
                         output_dir_path = "DE_results_2024/transcriptomics/2_batch_in_limma/p/",
                         plot_file_name = "PREOPEVsHC.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs HC",
                         output_dir_path = "DE_results_2024/transcriptomics/2_batch_in_limma/padj/",
                         plot_file_name = "PREOPEVsHC.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "MET - HC",
                              levels = colnames(model_matrix))
fit <- lmFit(v, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("tra")
result <- top.table %>%
  dplyr::select(tra, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = tra, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)
plot_volcano_and_save_DE(result, plot_title = "MET Vs HC",
                         output_dir_path = "DE_results_2024/transcriptomics/2_batch_in_limma/p/",
                         plot_file_name = "METVsHC.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "MET Vs HC",
                         output_dir_path = "DE_results_2024/transcriptomics/2_batch_in_limma/padj/",
                         plot_file_name = "METVsHC.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         plot_width_cm = 25)

##############################

#now with combat batch corrected and no batch specified in limma

data.combat <- ComBat(dat = data, batch = phenotype$data_cohort)
data.combat <- as.data.frame(data.combat)

model_matrix <- model.matrix(~ 0 + phenotype$PREOPE_MET_HC)
colnames(model_matrix) <- sub("phenotype$PREOPE_MET_HC", "", colnames(model_matrix), fixed = TRUE)

# v <- voom(data.combat, model_matrix, plot = TRUE)
v <- data.combat

contr_matrix <- makeContrasts(contrasts = "PREOPE - MET",
                              levels = colnames(model_matrix))
fit <- lmFit(v, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("tra")
result <- top.table %>%
  dplyr::select(tra, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = tra, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs MET",
                         output_dir_path = "DE_results_2024/transcriptomics/3_combat_corrected/p/",
                         plot_file_name = "PREOPEVsMET.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs MET",
                         output_dir_path = "DE_results_2024/transcriptomics/3_combat_corrected/padj/",
                         plot_file_name = "PREOPEVsMET.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "PREOPE - HC",
                              levels = colnames(model_matrix))
fit <- lmFit(v, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("tra")
result <- top.table %>%
  dplyr::select(tra, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = tra, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)
plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs HC",
                         output_dir_path = "DE_results_2024/transcriptomics/3_combat_corrected/p/",
                         plot_file_name = "PREOPEVsHC.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs HC",
                         output_dir_path = "DE_results_2024/transcriptomics/3_combat_corrected/padj/",
                         plot_file_name = "PREOPEVsHC.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "MET - HC",
                              levels = colnames(model_matrix))
fit <- lmFit(v, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("tra")
result <- top.table %>%
  dplyr::select(tra, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = tra, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)
plot_volcano_and_save_DE(result, plot_title = "MET Vs HC",
                         output_dir_path = "DE_results_2024/transcriptomics/3_combat_corrected/p/",
                         plot_file_name = "METVsHC.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "MET Vs HC",
                         output_dir_path = "DE_results_2024/transcriptomics/3_combat_corrected/padj/",
                         plot_file_name = "METVsHC.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         plot_width_cm = 25)

##########################################################################################
#with edgeR


data <- read.csv("Data/RNA_all/newquant_Nov2023_umi_counts_PREOPE_MET_HC_filter90.csv", row.names = 1)
phenotype <- read.table("Data/transcriptomic_phenotype_PREOPE_MET_HC_withaddicolumn.txt", header=TRUE, sep="\t")

data <- data[, phenotype$Sample]
all.equal(phenotype$Sample, colnames(data))
dim(data)
# [1] 907  67

phenotype <- phenotype %>%
  mutate(data_cohort = case_when(data_cohort == "initial" ~ "Cohort1",
                                 TRUE ~ 'Cohort2'))

keep <- edgeR::filterByExpr(data, group = phenotype$PREOPE_MET_HC)
data <- data[keep, ]

dim(data)
# [1] 207  67

data <- edgeR::cpm(data, log=TRUE)
data <- as.data.frame(data)

summary(factor(paste(phenotype$data_cohort, phenotype$PREOPE_MET_HC)))

data.combat <- ComBat(dat = data, batch = phenotype$data_cohort)
data.combat <- as.data.frame(data.combat)

y <- DGEList(counts = data.combat, group = phenotype$PREOPE_MET_HC)
model_matrix <- model.matrix(~ 0 + phenotype$PREOPE_MET_HC)
colnames(model_matrix) <- sub("phenotype$PREOPE_MET_HC", "", colnames(model_matrix), fixed = TRUE)
y <- estimateDisp(y, model_matrix)
fit <- glmQLFit(y, model_matrix)


contr_matrix <- makeContrasts(contrasts = "PREOPE - MET",
                              levels = colnames(model_matrix))
qlf <- glmQLFTest(fit, contrast = contr_matrix)
topTags(qlf)

result <- topTags(qlf, n = Inf)$table %>%
  rownames_to_column("tra") %>% 
  dplyr::select(tra, logFC, PValue, FDR) %>%
  dplyr::rename(Molecule = tra, adjPVal = FDR, PVal = PValue) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs MET",
                         output_dir_path = "DE_results_2024/transcriptomics/4_edgeR/p/",
                         plot_file_name = "PREOPEVsMET.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PREOPE Vs MET",
                         output_dir_path = "DE_results_2024/transcriptomics/4_edgeR/padj/",
                         plot_file_name = "PREOPEVsMET.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         plot_width_cm = 25)

#nothing significant or close to significant
#this isn't the usual recommended flow for edgeR since filter+logCPM+combat is done prior to creating DGEList and no further norm
#not using this

###############################
