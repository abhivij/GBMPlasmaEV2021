library(tidyverse)
library("ggvenn")

base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
setwd(base_dir)

x <- data
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
dataset_id <- "GBMPlasmaEV_transcriptomic_simple_norm_PREOPEVsMET"
phenotype_file_path <- "Data/transcriptomic_phenotype.txt"
plot_dir_path = "plots/FEMPipeline/selected_features_boxplot"

create_biomarker_boxplots <- function(dataset_id,
                                      phenotype_file_path,
                                      best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                                      plot_dir_path = "plots/FEMPipeline/selected_features_boxplot"){
  best_features <- read.csv(best_features_file_path)  
  best_features <- best_features %>%
    filter(dataset_id == !!dataset_id,
           is_best == 1) 
  biomarkers <- strsplit(best_features$biomarkers, split = "|", fixed = TRUE)[[1]]
  
  if(grepl(pattern = "transcriptomic", x = dataset_id, fixed = TRUE)){
    data <- read.table("Data/RNA/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                       nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")  
    norm <- "norm_log_cpm_simple"
    split_str <- "simple_norm_"
  } else if(grepl(pattern = "proteomic", x = dataset_id, fixed = TRUE)){
    norm <- "quantile"
    split_str <- "quantile_"
    if(grepl(pattern = "REC-TP", x = dataset_id, fixed = TRUE)){
      data <- read.table("Data/Protein/formatted_data/Q7_nonorm_formatted_impute50fil.csv", header=TRUE, sep=",", row.names=1, skip=0,
                         nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
    } else{
      data <- read.table("Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv", header=TRUE, sep=",", row.names=1, skip=0,
                         nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")  
    }
  } else{
    print("Unknown dataset_id type !")
    return
  }
  
  phenotype <- read.table(phenotype_file_path, header=TRUE, sep="\t")
  
  classification_criteria <- strsplit(dataset_id, split = split_str)[[1]]
  classification_criteria <- classification_criteria[length(classification_criteria)]
  
  extracted_samples <- phenotype[!is.na(phenotype[classification_criteria]), ]
  
  data <- data %>% select(extracted_samples$Sample) 
  
  norm_data <- normalize_data(data, norm)
  
  data_sub <- norm_data[gsub(".", "-", biomarkers, fixed = TRUE),]
  
  data_to_plot <- data_sub %>%
    rownames_to_column(var = "biomarker") %>%
    pivot_longer(cols = !biomarker, names_to = "Sample", values_to = "norm_expr") %>%
    inner_join(extracted_samples %>%
                 select(Sample, !!classification_criteria))
  
  data_to_plot <- data_to_plot %>%
    select(-c(Sample)) %>%
    rename("group" = !!classification_criteria)
  
  classes <- strsplit(classification_criteria, split = "Vs", fixed = TRUE)[[1]]
  data_to_plot <- data_to_plot %>%
    mutate(group = factor(group, levels = classes))
  
  ggplot(data_to_plot, aes(x = biomarker, 
                                   y = norm_expr,
                                   fill = group)) +
    geom_boxplot(size = 0.2, alpha = 0.5) +
    ggtitle(dataset_id) +
    ylab("normalized expression") +
    theme(axis.text.x = element_text(size=rel(1.2), angle = 45, hjust = 1),
          axis.text.y = element_text(size=rel(1.2)),
          axis.title.x = element_text(size=rel(1.5)),
          axis.title.y = element_text(size=rel(1.5)),
          plot.title  = element_text(size=rel(1.5)))    
  
  if(!dir.exists(plot_dir_path)){
    dir.create(plot_dir_path, recursive = TRUE)
  }
  
  plot_file_name <- paste0(plot_dir_path, "/",
                           dataset_id,
                           ".png")
  ggsave(plot_file_name)
  
}

best_features <- read.csv("Data/selected_features/best_features_with_add_col.csv")
for(dataset_id in unique(best_features$dataset_id)){
  print(dataset_id)
  if(grepl(pattern = "transcriptomic", x = dataset_id, fixed = TRUE)){
    phenotype_file_path <- "Data/transcriptomic_phenotype.txt"
  } else if(grepl(pattern = "proteomic", x = dataset_id, fixed = TRUE)){
    phenotype_file_path <- "Data/proteomic_phenotype.txt"
  }
  create_biomarker_boxplots(dataset_id = dataset_id,
                            phenotype_file_path = phenotype_file_path)
}




comparisons <- c("PREOPEVsMET", "PREOPEVsHC", "METVsHC")
best_features_file_path = "Data/selected_features/best_features_with_add_col.csv"
plot_dir_path = "plots/FEMPipeline/selected_features_venn"
omics_type = "proteomic"

#plot venn diagram of overlap of biomarkers between groups of interest
plot_biomarker_overlap_venn <- function(comparisons,
                                        omics_type = "transcriptomic",
                                        best_features_file_path = "Data/selected_features/best_features_with_add_col.csv",
                                        plot_dir_path = "plots/FEMPipeline/selected_features_venn"){
  best_features <- read.csv(best_features_file_path)  
  
  venn_diagram_list <- list()
  for(i in c(1 : length(comparisons))){
    if(omics_type == "transcriptomic"){
      dataset_id <- paste0("GBMPlasmaEV_transcriptomic_simple_norm_",
                           comparisons[i])
    }else{
      dataset_id <- paste0("GBMPlasmaEV_proteomic_impute50fil_quantile_",
                           comparisons[i])
    }
    best_features_sub <- best_features %>%
      filter(dataset_id == !!dataset_id,
             is_best == 1) 
    biomarkers <- strsplit(best_features_sub$biomarkers, split = "|", fixed = TRUE)[[1]]
    venn_diagram_list[[comparisons[i]]] <- biomarkers
  }
  
  if(length(comparisons) == 2){
    color_list <- c("lightskyblue", "mistyrose")
  } else if(length(comparisons) == 3){
    color_list <- c("lightskyblue", "mistyrose", "palegreen")
  } else{
    color_list <- c()
    print("invalid case !")
  }
  
  ggvenn(venn_diagram_list,
         fill_color = color_list,
         stroke_size = 0.1,
         set_name_size = 5,
         text_size = 3)
  if(!dir.exists(plot_dir_path)){
    dir.create(plot_dir_path, recursive = TRUE)
  }
  file_name <- paste0(omics_type, "_", paste(comparisons, collapse = "_"), ".png")
  ggsave(paste(plot_dir_path, file_name, sep = "/"))
}


plot_biomarker_overlap_venn(comparisons = c("PREOPEVsMET", 
                                            "PREOPEVsHC", 
                                            "METVsHC"),
                            omics_type = "transcriptomic")
plot_biomarker_overlap_venn(comparisons = c("PREOPEVsPOSTOPE_T", 
                                            "PREOPEVsPOSTOPE_P", 
                                            "POSTOPE_TVsPOSTOPE_P"),
                            omics_type = "transcriptomic")
plot_biomarker_overlap_venn(comparisons = c("POSTOPE_TVsREC_T", 
                                            "POSTOPE_PVsREC_P"),
                            omics_type = "transcriptomic")
plot_biomarker_overlap_venn(comparisons = c("POSTOPE_TVsREC_T", 
                                            "POSTOPE_TVsPREREC"),
                            omics_type = "transcriptomic")
plot_biomarker_overlap_venn(comparisons = c("POSTOPE_TVsREC_T", 
                                            "PREOPEVsREC_TP"),
                            omics_type = "transcriptomic")


plot_biomarker_overlap_venn(comparisons = c("PREOPEVsMET", 
                                            "PREOPEVsHC", 
                                            "METVsHC"),
                            omics_type = "proteomic")
plot_biomarker_overlap_venn(comparisons = c("PREOPEVsPOSTOPE_T", 
                                            "PREOPEVsPOSTOPE_P", 
                                            "POSTOPE_TVsPOSTOPE_P"),
                            omics_type = "proteomic")
plot_biomarker_overlap_venn(comparisons = c("POSTOPE_TVsREC_T", 
                                            "POSTOPE_PVsREC_P"),
                            omics_type = "proteomic")
plot_biomarker_overlap_venn(comparisons = c("POSTOPE_TVsREC_T", 
                                            "POSTOPE_TVsPREREC"),
                            omics_type = "proteomic")
plot_biomarker_overlap_venn(comparisons = c("POSTOPE_TVsREC_T", 
                                            "PREOPEVsREC_TP"),
                            omics_type = "proteomic")

