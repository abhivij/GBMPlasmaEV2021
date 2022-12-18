library(tidyverse)

validation_metadata <- read.csv("Data/metadata_glionet.csv")
data.new <- read.csv("Data/Protein/manual_check_for_quality/MSstatsInput-B1-3.13.7.csv")

data.new <- data.new %>%
  mutate(BioReplicate = sub("20220714_Eclipse-OT_SH_|20220719_Eclipse-OT_SH_", "", File.Name)) %>%
  mutate(BioReplicate = sub(".raw", "", BioReplicate, fixed = ""))

data.new <- data.new %>%
  left_join(validation_metadata %>% select(sample_id, category_old_name), by = c("BioReplicate" = "sample_id")) %>%
  mutate(Condition = case_when(is.na(category_old_name) ~ 'QC',
                               TRUE ~ category_old_name)) %>%
  select(-c(category_old_name))


data.grouped <- data.new %>%
  group_by(BioReplicate, Protein.Name, 
             Peptide.Modified.Sequence, Precursor.Charge,
             Fragment.Ion, Product.Charge) %>%
  summarize(count = n()) %>%
  filter(count > 1) %>%
  ungroup()

data.duplicates <- data.new %>%
  inner_join(data.grouped) %>%
  filter(!is.na(Area) & Area != 0) %>%
  arrange(desc(count),
          BioReplicate, Protein.Name, 
          Peptide.Modified.Sequence, Precursor.Charge, 
          Fragment.Ion, Product.Charge) %>%
  select(c(BioReplicate, Protein.Name, 
           Peptide.Modified.Sequence, Precursor.Charge, 
           Fragment.Ion, Product.Charge,
           count, Condition,
           File.Name, Detection.Q.Value, Area,
           Standard.Type, Isotope.Label.Type))

write.csv(data.new, "Data/Protein/manual_check_for_quality/MSstatsInput-B1-3.13.7_modified.csv",
          row.names = FALSE)  


#################


create_dim_red_plots(comparison = "PREOPEVsPOSTOPE_TP", classes = c("POSTOPE_TP", "PREOPE"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "quantile_train_param",
                     data_to_show = "initial",
                     show_only_common = FALSE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)
create_dim_red_plots(comparison = "PREOPEVsPOSTOPE_TP", classes = c("POSTOPE_TP", "PREOPE"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "quantile_train_param",
                     data_to_show = "initial",
                     show_only_common = TRUE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)
create_dim_red_plots(comparison = "PREOPEVsPOSTOPE_TP", classes = c("POSTOPE_TP", "PREOPE"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "quantile_train_param",
                     data_to_show = "validation",
                     show_only_common = FALSE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)
create_dim_red_plots(comparison = "PREOPEVsPOSTOPE_TP", classes = c("POSTOPE_TP", "PREOPE"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "quantile_train_param",
                     data_to_show = "validation",
                     show_only_common = TRUE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)
create_dim_red_plots(comparison = "PREOPEVsPOSTOPE_TP", classes = c("POSTOPE_TP", "PREOPE"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "quantile_train_param",
                     data_to_show = "validation",
                     show_only_common = FALSE,
                     show_imputed_based_on_initial = TRUE,
                     shownames = FALSE,
                     perform_filter = FALSE)
create_dim_red_plots(comparison = "PREOPEVsPOSTOPE_TP", classes = c("POSTOPE_TP", "PREOPE"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "quantile_train_param",
                     data_to_show = "both",
                     show_only_common = TRUE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)
create_dim_red_plots(comparison = "PREOPEVsPOSTOPE_TP", classes = c("POSTOPE_TP", "PREOPE"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "quantile_train_param",
                     data_to_show = "both",
                     show_only_common = FALSE,
                     show_imputed_based_on_initial = TRUE,
                     shownames = FALSE,
                     perform_filter = FALSE)

create_dim_red_plots(comparison = "PREOPEVsPOSTOPE_TP", classes = c("POSTOPE_TP", "PREOPE"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "none",
                     data_to_show = "initial",
                     show_only_common = FALSE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)
create_dim_red_plots(comparison = "PREOPEVsPOSTOPE_TP", classes = c("POSTOPE_TP", "PREOPE"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "none",
                     data_to_show = "initial",
                     show_only_common = TRUE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)
create_dim_red_plots(comparison = "PREOPEVsPOSTOPE_TP", classes = c("POSTOPE_TP", "PREOPE"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "none",
                     data_to_show = "validation",
                     show_only_common = FALSE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)
create_dim_red_plots(comparison = "PREOPEVsPOSTOPE_TP", classes = c("POSTOPE_TP", "PREOPE"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "none",
                     data_to_show = "validation",
                     show_only_common = TRUE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)
create_dim_red_plots(comparison = "PREOPEVsPOSTOPE_TP", classes = c("POSTOPE_TP", "PREOPE"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "none",
                     data_to_show = "both",
                     show_only_common = TRUE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)

create_dim_red_plots(comparison = "PREOPEVsREC_TP", classes = c("REC_TP", "PREOPE"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "quantile_train_param",
                     data_to_show = "initial",
                     show_only_common = FALSE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)
create_dim_red_plots(comparison = "PREOPEVsREC_TP", classes = c("REC_TP", "PREOPE"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "quantile_train_param",
                     data_to_show = "initial",
                     show_only_common = TRUE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)
create_dim_red_plots(comparison = "PREOPEVsREC_TP", classes = c("REC_TP", "PREOPE"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "quantile_train_param",
                     data_to_show = "validation",
                     show_only_common = FALSE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)
create_dim_red_plots(comparison = "PREOPEVsREC_TP", classes = c("REC_TP", "PREOPE"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "quantile_train_param",
                     data_to_show = "validation",
                     show_only_common = TRUE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)
create_dim_red_plots(comparison = "PREOPEVsREC_TP", classes = c("REC_TP", "PREOPE"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "quantile_train_param",
                     data_to_show = "both",
                     show_only_common = TRUE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)

create_dim_red_plots(comparison = "POSTOPE_TPVsREC_TP", classes = c("REC_TP", "POSTOPE_TP"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "quantile_train_param",
                     data_to_show = "initial",
                     show_only_common = FALSE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)
create_dim_red_plots(comparison = "POSTOPE_TPVsREC_TP", classes = c("REC_TP", "POSTOPE_TP"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "quantile_train_param",
                     data_to_show = "initial",
                     show_only_common = TRUE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)
create_dim_red_plots(comparison = "POSTOPE_TPVsREC_TP", classes = c("REC_TP", "POSTOPE_TP"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "quantile_train_param",
                     data_to_show = "validation",
                     show_only_common = FALSE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)
create_dim_red_plots(comparison = "POSTOPE_TPVsREC_TP", classes = c("REC_TP", "POSTOPE_TP"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "quantile_train_param",
                     data_to_show = "validation",
                     show_only_common = TRUE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)
create_dim_red_plots(comparison = "POSTOPE_TPVsREC_TP", classes = c("REC_TP", "POSTOPE_TP"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "quantile_train_param",
                     data_to_show = "both",
                     show_only_common = TRUE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)

create_dim_red_plots(comparison = "PREOPEVsPOSTOPE_TP", classes = c("POSTOPE_TP", "PREOPE"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "none",
                     data_to_show = "both",
                     show_only_common = TRUE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)
create_dim_red_plots(comparison = "PREOPEVsREC_TP", classes = c("REC_TP", "PREOPE"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "none",
                     data_to_show = "both",
                     show_only_common = TRUE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)
create_dim_red_plots(comparison = "POSTOPE_TPVsREC_TP", classes = c("REC_TP", "POSTOPE_TP"),
                     omics_type = "proteomics",
                     dim_red = "UMAP", norm = "none",
                     data_to_show = "both",
                     show_only_common = TRUE,
                     show_imputed_based_on_initial = FALSE,
                     shownames = FALSE,
                     perform_filter = FALSE)
