#to check the expression levels of proteins and transcripts in PREOPE, MET, HC

library(tidyverse)
library(ggrepel)
library(ggvenn)
library(sva)
library(ggplot2)


comparison = "PREOPEVsHC"
classes = c("HC", "PREOPE")
omics_type = "transcriptomics"
norm = "log_cpm"
perform_filter = TRUE
batch_effect_correction = "combat"
plot_dir_path = "plots_RNA_all/PREOPE_MET_HC/qc/expression_sample_level"
data_file_path = "Data/RNA_all/newquant_Nov2023_umi_counts_PREOPE_MET_HC_filter90.csv"
validation_data_file_path = "Data/RNA_all/newquant_Nov2023_umi_counts_PREOPE_MET_HC_filter90.csv"
phenotype_file_path = "Data/transcriptomic_phenotype_PREOPE_MET_HC.txt"
plot_type = "per_sample"

create_per_feature_or_sample_plot <- function(comparison, classes,
                                              omics_type, norm,
                                              perform_filter = FALSE,
                                              batch_effect_correction = "none",
                                              plot_type = "per_feature",
                                              plot_dir_path,
                                              data_file_path = NA,
                                              validation_data_file_path = NA,
                                              phenotype_file_path = NA){
  if(omics_type == "proteomics"){
    if(is.na(data_file_path)){
      data_file_path <- "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv"
    }
    if(is.na(validation_data_file_path)){
      validation_data_file_path <- "Data/Protein/formatted_data/newcohort_nonorm_formatted_impute50fil.csv"
    }
    if(is.na(phenotype_file_path)){
      phenotype_file_path <- "Data/proteomic_phenotype_PREOPE_MET_HC.txt"  
    }
  } else if(omics_type == "transcriptomics"){
    
    if(is.na(data_file_path)){
      data_file_path <- "Data/RNA/umi_counts_initial_cohort.csv"
    }
    if(is.na(validation_data_file_path)){
      validation_data_file_path <- "Data/RNA/umi_counts_validation_cohort.csv"
    }
    if(is.na(phenotype_file_path)){
      phenotype_file_path <- "Data/transcriptomic_phenotype_PREOPE_MET_HC.txt"  
    }
    # validation_metadata <- read.csv("Data/RNA_validation/metadata_glionet.csv") %>%
    #   mutate(sample_category = factor(sample_category)) %>%
    #   mutate(sample_category = recode_factor(sample_category, "PRE-OP" = "PREOPE",
    #                                          "POST-OP" = "POSTOPE_TP",
    #                                          "RECURRENCE" = "REC_TP"))
  }
  
  data <- read.csv(data_file_path, row.names = 1)
  validation_data <- read.csv(validation_data_file_path, row.names = 1)
  
  phenotype <- read.table(phenotype_file_path, header=TRUE, sep="\t")
  
  if(omics_type == "proteomics"){
    colnames(validation_data)[colnames(validation_data) == "SB12_01"] = "SB12"
    #use SB22.02
    colnames(validation_data)[colnames(validation_data) == "SB22.02"] = "SBtobeused22"
    colnames(validation_data)[colnames(validation_data) == "SB22"] = "SB22_dont_include"
    colnames(validation_data)[colnames(validation_data) == "SBtobeused22"] = "SB22"
    
    #SB7 sample not required for this analysis
    # so no need to filter out
    # validation_metadata <- validation_metadata %>%
    #   filter(Sample != "SB7")    
    
  } else if(omics_type == "transcriptomics" & data_file_path != validation_data_file_path){
    colnames(validation_data) <- paste0("S", colnames(validation_data))
  }
  
  title <- paste0(norm, " BEC-", batch_effect_correction)
  
  output_labels.cohort1 <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes) %>%
    dplyr::select(Sample, Label, data_cohort, Subgroup, Sex, Age) %>%
    filter(data_cohort == "initial")
  output_labels.cohort2 <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes) %>%
    dplyr::select(Sample, Label, data_cohort, Subgroup, Sex, Age) %>%
    filter(data_cohort == "validation")
  
  #currently data format : (transcripts x samples)
  
  data.cohort1 <- data %>% dplyr::select(output_labels.cohort1$Sample)
  data.cohort2 <- validation_data %>% dplyr::select(output_labels.cohort2$Sample)
  
  #take common only if both cohorts contain samples
  if(ncol(data.cohort1) > 0 & ncol(data.cohort2) > 0){
    common <- intersect(rownames(data.cohort1), rownames(data.cohort2))  
    data.cohort1 <- data.cohort1[common, ]
    data.cohort2 <- data.cohort2[common, ]
    data <- cbind(data.cohort1, data.cohort2)    
  } else if(ncol(data.cohort1) > 0){
    data <- data.cohort1
  } else if(ncol(data.cohort2) > 0){
    data <- data.cohort2
  } else{
    print("no samples")
  }
  
  output_labels <- rbind(output_labels.cohort1, output_labels.cohort2)
  
  if(perform_filter){
    keep <- edgeR::filterByExpr(data, group = output_labels$Label)
    data_left_out <- data[!keep, ]
    data <- data[keep, ]
  }
  
  if(norm == "quantile"){
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
    data <- data.norm
    
  } else if(norm == "log_cpm"){
    #calculating norm log cpm
    data <- edgeR::cpm(data, log=TRUE)
  } else if(norm == "none"){
    print("just log")
    data <- log2(data + 2^-10)
  }
  
  data <- as.data.frame(t(as.matrix(data)))
  
  output_labels <- output_labels %>%
    mutate(Label = factor(Label), data_cohort = factor(data_cohort))
  
  group_counts <- output_labels %>%
    dplyr::mutate(Label = paste(data_cohort, Label, sep = "_")) %>%
    group_by(Label) %>%
    summarise(n = n())
  
  group_counts_text <- paste(apply(group_counts, MARGIN = 1, FUN = function(x){paste(x[1], x[2], sep = ":")}),
                             collapse = "  ")
  
  all.equal(rownames(data), output_labels$Sample)
    
  if(batch_effect_correction == "combat"){
    data <- as.data.frame(t(as.matrix(data)))
    data.combat = ComBat(dat=data, batch=output_labels$data_cohort)
    data.combat <- as.data.frame(t(as.matrix(data.combat)))
    data <- data.combat
  }
  
  if(plot_type == "per_feature"){
    
    data_to_plot <- data %>%
      rownames_to_column(var = "Sample") %>%
      pivot_longer(cols = !Sample, names_to = "features", values_to = "expr") %>%
      inner_join(output_labels %>%
                   mutate(Label = paste(Label, data_cohort, sep = "_")) %>%
                   dplyr::select(c(Sample, Label))
      )
    
    data_to_plot <- data_to_plot %>%
      group_by(features, Label) %>%
      summarize(med_expr = median(expr), 
                q25up = quantile(expr, probs = 0.75), 
                q25down = quantile(expr, probs = 0.25)) %>%
      ungroup()
    
    median_expr_data <- data_to_plot %>%
      group_by(features) %>%
      summarize(med_med_expr = median(med_expr)) %>%
      arrange(desc(med_med_expr))
    
    data_to_plot <- data_to_plot %>%
      pivot_longer(cols = c(med_expr, q25up, q25down), names_to = "expr_type", values_to = "expr") %>%
      mutate(features = factor(features, levels = median_expr_data$features))
    
    data_to_plot <- data_to_plot %>%
      dplyr::filter(expr_type == "med_expr")
    
    y_lab <- "Expression"
    if(omics_type == "transcriptomics"){
      x_lab <- "transcripts"
    } else{
      x_lab <- "proteins"
      
      # data_to_plot <- data_to_plot %>%
      #   left_join(all_protein_names %>% dplyr::select(c(from_id, primary_gene_id)), 
      #             by = c("biomarker" = "from_id")) %>%
      #   dplyr::select(-c(biomarker)) %>%
      #   dplyr::rename(c("biomarker" = "primary_gene_id"))
    }
    x_lab <- paste0(x_lab, "(", length(median_expr_data$features), ")")
    
    # biomarker_agg <- data_to_plot %>%
    #   group_by(biomarker) %>%
    #   summarize(med_expr = median(norm_expr)) %>%
    #   arrange(desc(med_expr))
    # 
    # # data_to_plot <- data_to_plot %>%
    # #   mutate(Label = factor(Label, levels = rev(classes)),
    # #          biomarker = factor(biomarker, biomarker_agg$biomarker)) 
    # data_to_plot <- data_to_plot %>%
    #   mutate(biomarker = factor(biomarker, biomarker_agg$biomarker))
    # # print("here")
    
    # ggplot(data_to_plot) +
    #   geom_point(aes(x = features, y = expr, fill = Label, shape = expr_type), stroke = 0.1) +
    #   theme(axis.text.x = element_blank()) +
    #   scale_shape_manual(name = "Expression type", values = c(21, 24, 25)) +
    #   guides(fill = guide_legend(override.aes = list(shape = 21, colour = NA))) +
    #   ggtitle(paste(sub("Vs", " Vs ", comparison), x_lab)) +
    #   xlab(x_lab) +
    #   ylab(y_lab)
    
    ggplot(data_to_plot) +
      geom_point(aes(x = features, y = expr, color = Label)) +
      theme(axis.text.x = element_blank()) +
      ggtitle(paste(sub("Vs", " Vs ", comparison), x_lab, title)) +
      xlab(x_lab) +
      ylab(y_lab)    
  }  else if(plot_type == "per_sample"){
    
    data_to_plot <- data %>%
      rownames_to_column(var = "Sample") %>%
      pivot_longer(cols = !Sample, names_to = "features", values_to = "expr") %>%
      inner_join(output_labels %>%
                   mutate(Label = paste(Label, data_cohort, sep = "_")) %>%
                   dplyr::select(c(Sample, Label))
      )
    data_to_plot <- data_to_plot %>%
      mutate(Label = factor(Label))
    
    median_expr_data <- data_to_plot %>%
      group_by(Sample, Label) %>%
      summarize(med_expr = median(expr)) %>%
      arrange(Label, desc(med_expr))
    
    data_to_plot <- data_to_plot %>%
      mutate(Sample = factor(Sample, levels = median_expr_data$Sample))
    
    y_lab <- "Expression"
    x_lab <- "Sample"
    # x_lab <- paste0(x_lab, "(", length(data_to_plot$Sample), ")")
    
    ggplot(data_to_plot) +
      geom_boxplot(aes(x = Sample, y = expr, color = Label)) +
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.99)) +
      ggtitle(paste(sub("Vs", " Vs ", comparison), x_lab, title)) +
      xlab("Sample") +
      ylab(y_lab)    
  }

  
  if(!dir.exists(plot_dir_path)){
    dir.create(plot_dir_path, recursive = TRUE)
  }
  file_name <- paste(comparison, x_lab, title, ".png", sep = "_")
  file_path <- paste(plot_dir_path, file_name, sep = "/")
  ggplot2::ggsave(file_path, units = "cm", width = 30)
  
}



