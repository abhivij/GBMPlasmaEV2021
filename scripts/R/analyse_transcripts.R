library(tidyverse)
library(readxl)
library(edgeR)
library(caret)


base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

source("scripts/R/utils.R")
source("scripts/R/plot_data.R")


generated_mirna_counts <- read.table("mirna_counts.csv", header=TRUE, sep = ",", row.names = 1)

data <- read_excel("Data/RNA/158629.all_samples.summary.xlsx", sheet = "miRNA_piRNA")
mirna_data <- data[1:2505,]
pirna_data <- data[2507:2642,] %>%
  separate(miRNA, c("miRNA", NA, NA, NA), sep = "/")

dim(data)[1]
dim(mirna_data)[1] + dim(pirna_data)[1]

data <- rbind(mirna_data, pirna_data) %>%
  column_to_rownames("miRNA")

umi_counts <- data %>%
  select(ends_with("UMIs"))
all_group_names <- sapply(colnames(umi_counts),
                          function(x){
                            prefix <- substr(x, 1, 2)
                            if (prefix == "ME") {
                              label <- "MET"
                            } else if (prefix == "HB") {
                              label <- "GBM"
                            } else {
                              label <- "HC"
                            }
                            return (label)
                          }
)

plot_transcriptomic_data <- function(data, file_name, title, groups = NA, dim_red = "pca"){
  if(length(groups) <= 1 && is.na(groups)){
    groups <- sapply(colnames(data),
                              function(x){
                                prefix <- substr(x, 1, 2)
                                if (prefix == "ME") {
                                  label <- "MET"
                                } else if (prefix == "HB") {
                                  label <- "GBM"
                                } else {
                                  label <- "HC"
                                }
                                return (label)
                              }
    )  
  }
  plot_data(t(data), file_name, title, groups = groups, dim_red = dim_red,
            colour_label = "Labels")
}

plot_transcriptomic_data(umi_counts, "umi_counts.jpg", "UMIs")
plot_transcriptomic_data(umi_counts, "umi_counts.jpg", "UMIs", dim_red = "tsne")
plot_transcriptomic_data(umi_counts, "umi_counts.jpg", "UMIs", dim_red = "umap")




keep <- filterByExpr(umi_counts, group = all_group_names)
umi_counts <- umi_counts[keep, ]

plot_transcriptomic_data(umi_counts, "filtered_umi_counts.jpg", "Filtered UMIs")
plot_transcriptomic_data(umi_counts, "filtered_umi_counts.jpg", "Filtered UMIs", dim_red = "tsne")
plot_transcriptomic_data(umi_counts, "filtered_umi_counts.jpg", "Filtered UMIs", dim_red = "umap")


umi_norm_data <- cpm(umi_counts, log = TRUE)
umi_norm_data <- scale(umi_norm_data)

umi_norm_data <- as.data.frame(t(as.matrix(umi_norm_data)))
umi_norm_data <- predict(preProcess(umi_norm_data), umi_norm_data)
umi_norm_data <- as.data.frame(t(as.matrix(umi_norm_data)))

plot_transcriptomic_data(umi_norm_data, "normlogcpm_umi_counts.jpg", "Norm Log CPM UMIs")
plot_transcriptomic_data(umi_norm_data, "normlogcpm_umi_counts.jpg", "Norm Log CPM UMIs", dim_red = "tsne")
plot_transcriptomic_data(umi_norm_data, "normlogcpm_umi_counts.jpg", "Norm Log CPM UMIs", dim_red = "umap")


all_group_names <- sapply(colnames(umi_norm_data),
                          function(x){
                            prefix <- substr(x, 1, 2)
                            if (prefix == "ME") {
                              label <- "MET"
                            } else if (prefix == "HB") {
                              label <- "GBM"
                            } else {
                              label <- "HC"
                            }
                            return (label)
                          }
)

ggplot() +
  geom_histogram(aes(x = umi_norm_data[,3]), stat = "bin", bins = 200) +
  xlab("UMIs") +
  ylab("Normalized counts")
ggsave("umi_hist.jpg", width = 10, height = 10)


ggplot() +
  geom_histogram(aes(x = umi_counts[,15]), stat = "bin", bins = 200) +
  xlab("UMIs") +
  ylab("Normalized counts")
ggsave("umi_hist.jpg", width = 10, height = 10)



###################

plot_data(generated_mirna_counts, "read_counts.jpg", "Read counts")


all_group_names <- sapply(colnames(generated_mirna_counts),
                          function(x){
                            prefix <- substr(x, 1, 2)
                            if (prefix == "ME") {
                              label <- "MET"
                            } else if (prefix == "HB") {
                              label <- "GBM"
                            } else {
                              label <- "HC"
                            }
                            return (label)
                          }
)
keep <- filterByExpr(generated_mirna_counts, group = all_group_names)
generated_mirna_counts <- generated_mirna_counts[keep, ]
plot_data(generated_mirna_counts, "filtered_read_counts.jpg", "Filtered Reads")


norm_data <- cpm(generated_mirna_counts, log = TRUE)
norm_data <- scale(norm_data)

norm_data <- as.data.frame(t(as.matrix(norm_data)))
norm_data <- predict(preProcess(norm_data), norm_data)
norm_data <- as.data.frame(t(as.matrix(norm_data)))

plot_data(norm_data, "normlogcpm_read_counts.jpg", "Norm Log CPM Reads")


###############

#DE analysis

data <- read_excel("Data/RNA/158629.all_samples.summary.xlsx", sheet = "miRNA_piRNA")
mirna_data <- data[1:2505,]
pirna_data <- data[2507:2642,] %>%
  separate(miRNA, c("miRNA", NA, NA, NA), sep = "/")

dim(data)[1]
dim(mirna_data)[1] + dim(pirna_data)[1]

data <- rbind(mirna_data, pirna_data) %>%
  column_to_rownames("miRNA")

umi_counts <- data %>%
  select(ends_with("UMIs"))
colnames(umi_counts) <- gsub("-UMIs", "", colnames(umi_counts))
colnames(umi_counts) <- sapply(colnames(umi_counts), FUN = 
         function(x){
           strsplit(x, split = "_", fixed = TRUE)[[1]][1]
         }
)
#to get meta data start

protein_data <- read.csv(file = "Data/Protein/output/norm_annotatedQ1-6_NA_equalizeMedians.csv")

group_mapping <- protein_data %>%
  select(SUBJECT_ORIGINAL, GROUP_ORIGINAL) %>%
  arrange(SUBJECT_ORIGINAL)
rm(protein_data)

metadata <- data.frame(SUBJECT_ORIGINAL = colnames(umi_counts)) %>%
  arrange(SUBJECT_ORIGINAL)

#in proteomics data, unlike transcriptomics, HB numbered 01, 02, 03, ...
# so changing that to match transcriptomics sample names
group_mapping <- group_mapping %>%
  mutate(SUBJECT_ORIGINAL = gsub("HB0", "HB", SUBJECT_ORIGINAL))

metadata <- metadata %>%
  left_join(group_mapping)

group <- gsub("-", "_", metadata$GROUP_ORIGINAL)

#meta data end



data_for_de <- umi_counts
keep <- filterByExpr(data_for_de, group=group)
data_for_de <- data_for_de[keep, ]

dge_data <- DGEList(data_for_de)
dge_data <- calcNormFactors(dge_data)

model_matrix <- model.matrix(~0 + group)
colnames(model_matrix) <- gsub("group", "", colnames(model_matrix))

y <- voom(dge_data, model_matrix, plot = TRUE)
# y
vfit <- lmFit(y, model_matrix)
# head(coef(vfit))

contr_matrix <- makeContrasts(contrasts = "MET - HC",
                              levels = colnames(model_matrix))
vfit <- contrasts.fit(vfit, contr_matrix)
efit <- eBayes(vfit)
plotSA(efit)

dt <- decideTests(efit)
summary(dt)

top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("rna")

result <- top.table %>%
  select(rna, logFC, adj.P.Val) %>%
  rename(Molecule = rna, adjPVal = adj.P.Val)

create_volcano_plot(result, 
                    title = "MET Vs HC", 
                    file_name = "de_rna.png", 
                    dir_path = "plots", logFC_cutoff = 2)


###############

data_for_de <- umi_counts
keep <- filterByExpr(data_for_de, group=group)
data_for_de <- data_for_de[keep, ]

dge_data <- DGEList(data_for_de)
dge_data <- calcNormFactors(dge_data)

model_matrix <- model.matrix(~0 + group)
colnames(model_matrix) <- gsub("group", "", colnames(model_matrix))

y <- voom(dge_data, model_matrix, plot = TRUE)
# y
vfit <- lmFit(y, model_matrix)
# head(coef(vfit))

contr_matrix <- makeContrasts(contrasts = "PREOPE - HC",
                              levels = colnames(model_matrix))
vfit <- contrasts.fit(vfit, contr_matrix)
efit <- eBayes(vfit)
plotSA(efit)

dt <- decideTests(efit)
summary(dt)

top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("rna")

result <- top.table %>%
  select(rna, logFC, adj.P.Val) %>%
  rename(Molecule = rna, adjPVal = adj.P.Val)

create_volcano_plot(result, 
                    title = "PREOPE Vs HC", 
                    file_name = "de_rna_ph.png", 
                    dir_path = "plots", logFC_cutoff = 2)



################


data_for_de <- umi_counts
keep <- filterByExpr(data_for_de, group=group)
data_for_de <- data_for_de[keep, ]

dge_data <- DGEList(data_for_de)
dge_data <- calcNormFactors(dge_data)

model_matrix <- model.matrix(~0 + group)
colnames(model_matrix) <- gsub("group", "", colnames(model_matrix))

y <- voom(dge_data, model_matrix, plot = TRUE)
# y
vfit <- lmFit(y, model_matrix)
# head(coef(vfit))

contr_matrix <- makeContrasts(contrasts = "PREOPE - MET",
                              levels = colnames(model_matrix))
vfit <- contrasts.fit(vfit, contr_matrix)
efit <- eBayes(vfit)
plotSA(efit)

dt <- decideTests(efit)
summary(dt)

top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("rna")

result <- top.table %>%
  select(rna, logFC, adj.P.Val) %>%
  rename(Molecule = rna, adjPVal = adj.P.Val)

create_volcano_plot(result, 
                    title = "PREOPE Vs MET", 
                    file_name = "de_rna_pm.png", 
                    dir_path = "plots", logFC_cutoff = 2)

