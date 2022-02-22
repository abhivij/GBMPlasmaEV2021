library(tidyverse)

# base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
# setwd(base_dir)

source("scripts/R/results_analysis/results_analysis_utils.R")


# comparison <- "METVsHC"
# omics_type = "transcriptomic"
# phenotype_column = "HC_MET_PREOPE_REC_TP"
# conditions <- c("HC", "MET", "PREOPE", "REC_TP")
# best_features_file_path <- "Data/selected_features/best_features_with_add_col.csv"
# plot_dir_path <- "plots/FEMPipeline/biomarker_linegraph/HC_MET_PREOPE_REC_TP"

#plot linegraph with median value of biomarkers as different lines against different conditions
#also plots the associated boxplots to show the variation os biomarkers within the samples
plot_biomarker_linegraph <- function(conditions,
                                     comparison,
                                     omics_type = "transcriptomic",
                                     phenotype_column = "GROUP_Q1to6",
                                     best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                                     plot_dir_path = "plots/FEMPipeline/biomarker_linegraph"){
  
  if(!dir.exists(plot_dir_path)){
    dir.create(plot_dir_path, recursive = TRUE)
  }
  
  best_features <- read.csv(best_features_file_path)  
  
  categories <- strsplit(comparison, split = "Vs", fixed = TRUE)[[1]]
  if(omics_type == "transcriptomic"){
    dataset_id <- paste0("GBMPlasmaEV_transcriptomic_simple_norm_",
                         comparison)
  }else{
    dataset_id <- paste0("GBMPlasmaEV_proteomic_impute50fil_quantile_",
                         comparison)
  }
  best_features_sub <- best_features %>%
    filter(dataset_id == !!dataset_id,
           is_best == 1) 
  biomarkers <- strsplit(best_features_sub$biomarkers, split = "|", fixed = TRUE)[[1]]  
  
  if(omics_type == "transcriptomic"){
    data <- read.table("Data/RNA/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                       nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")  
    norm <- "norm_log_cpm_simple"
    split_str <- "simple_norm_"
    phenotype <- read.table("Data/transcriptomic_phenotype.txt", header=TRUE, sep="\t")
  } else {
    norm <- "quantile"
    split_str <- "quantile_"
    phenotype <- read.table("Data/proteomic_phenotype.txt", header=TRUE, sep="\t")
    if(grepl(pattern = "REC-TP", x = dataset_id, fixed = TRUE)){
      #currently this case will cause issues while reading phenotype 
      # and requiring multiple columns in phenotype
      data <- read.table("Data/Protein/formatted_data/Q7_nonorm_formatted_impute50fil.csv", header=TRUE, sep=",", row.names=1, skip=0,
                         nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
    } else{
      data <- read.table("Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv", header=TRUE, sep=",", row.names=1, skip=0,
                         nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")  
    }
  } 
  
  extracted_samples <- phenotype %>%
    rename("category" = phenotype_column) %>%
    filter(category %in% conditions) %>%
    select(Sample, category)
  
  data <- data %>% select(extracted_samples$Sample) 
  
  norm_data <- normalize_data(data, norm)
  
  data_to_plot <- data.frame(t(norm_data))  #"-" in column names get converted to "."
                                            #while converting to data.frame()
  
  data_to_plot <- data_to_plot %>%
    select(biomarkers) %>%
    rownames_to_column("Sample") %>%
    inner_join(extracted_samples) %>%
    select(-c(Sample))
  
  data_to_plot_boxplot <- data_to_plot %>%
    pivot_longer(cols = !category, names_to = "biomarker", values_to = "norm_expr")
  ggplot(data_to_plot_boxplot, aes(x = biomarker, 
                                   y = norm_expr,
                                   fill = category)) +
    geom_boxplot(size = 0.2, alpha = 0.5) +
    stat_summary(fun = "mean", geom="point", shape=20, size=2, color="blue",
                 position = position_dodge2(width = 0.75,   
                                            preserve = "single")) +
    ggtitle(paste("Biomarkers from", comparison, "comparison", 
                  "from", omics_type, "data")) +
    ylab("normalized expression") +
    theme(axis.text.x = element_text(size=rel(1.2), angle = 45, hjust = 1),
          axis.text.y = element_text(size=rel(1.2)),
          axis.title.x = element_text(size=rel(1.5)),
          axis.title.y = element_text(size=rel(1.5)),
          plot.title  = element_text(size=rel(1.5)))
  
  plot_file_name <- paste0(plot_dir_path, "/",
                           comparison, 
                           "_", omics_type,
                           "_", "boxplot",
                           ".png")
  ggsave(plot_file_name)
  
  
  data_to_plot <- data_to_plot %>%
    group_by(category) %>%
    summarise_all(median)
  
  data_to_plot <- data_to_plot %>%
    pivot_longer(cols = !category, names_to = "biomarker", values_to = "median_norm_expr")
  
  ggplot(data_to_plot, aes(x = category, y = median_norm_expr, 
                           group = biomarker, color = biomarker)) +
    geom_line(alpha = 0.5, linetype = 2) +
    geom_point(size = 2) +
    xlab("Condition") +
    ylab("Normalized expression median across samples") +
    ggtitle(paste("Biomarkers from", comparison, "comparison", 
                  "from", omics_type, "data")) +
    theme(axis.text.x = element_text(size=rel(1.2)),
          axis.text.y = element_text(size=rel(1.2)),
          axis.title.x = element_text(size=rel(1.3)),
          axis.title.y = element_text(size=rel(1.3)),
          plot.title  = element_text(size=rel(1.5)))
  
  plot_file_name <- paste0(plot_dir_path, "/",
                           comparison, 
                           "_", omics_type,
                           "_", "connectedscatter",
                           ".png")
  ggsave(plot_file_name)  
}  