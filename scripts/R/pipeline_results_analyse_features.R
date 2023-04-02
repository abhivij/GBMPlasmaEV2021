setwd("~/UNSW/VafaeeLab/GBMPlasmaEV/")
source("scripts/R/dataset_pipeline_arguments.R")
source("scripts/R/dataset_pipeline_arguments_transcriptomic.R")
source("scripts/R/dataset_pipeline_arguments_proteomic.R")
source("scripts/R/utils.R")
library(tidyverse)
library(viridis)
library(ComplexHeatmap)
library(UpSetR)


# best_fsm <- "mrmr100"

dparg_id = 41
best_fsm_vec = c("t-test",
                 "wilcoxontest",
                 "mrmr75",
                 "wilcoxontest_pval_0.005"
)
min_iter_feature_presence = 28

dparg_id = 139
best_fsm_vec = c("t-test_pval_0.025",
                 "mrmr100",
                 "wilcoxontest_BH",
                 "ranger_impu_cor"
)
min_iter_feature_presence = 28



explore_common_features <- function(dparg_id, best_fsm_vec, 
                                    min_iter_feature_presence,
                                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                                    results_dir = "fem_pipeline_results",
                                    dir_path = "plots/FEMPipeline/common_features_upset"){
  
  ds <- dataset_pipeline_arguments[[dparg_id]]
  dataset_id <- paste(ds$dataset_id, ds$classification_criteria, sep = "_")
  print(dataset_id)
  
  features_file <- paste(dataset_id, "features.csv", sep = "_")
  features_file <- paste(results_dir, features_file, sep = "/")
  print(features_file)
  
  selected_features <- list()
  max_size <- 0  #for upset plot max bar length
  for(best_fsm in best_fsm_vec){
    features_info <- read.table(features_file, sep = ',', header = TRUE)
    
    features_info <- features_info %>%
      filter(FSM == best_fsm) %>%
      select(-c(FSM, Iter))
    print(best_fsm)
    print(dim(features_info))
    
    print(dim(features_info[,colSums(features_info) >= min_iter_feature_presence, drop = FALSE]))
    selected_features[[best_fsm]] <- colnames(features_info[,colSums(features_info) >= min_iter_feature_presence, drop = FALSE])   
    
    size <- length(selected_features[[best_fsm]])
    if(size > max_size){
      max_size <- size
    }
  }
  ###########################################
  #upsetplot
  file_name <- paste(dataset_id, min_iter_feature_presence, 
                     "upset.png", sep = "_")
  
  if(!dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE)
  }
  
  print("creating upset plot...")
  print(paste(dir_path, file_name, sep = "/"))
  # dev.copy(png, paste(dir_path, file_name, sep = "/"), 
  #          units = "cm", width = 30, height = 15, res = 1200)
  png(filename = paste(dir_path, file_name, sep = "/"),
      units = "cm", width = 30, height = 15, res = 1200)
  print({
    upset(fromList(selected_features), set_size.show = TRUE,  
          set_size.scale_max = max_size + max_size/10)  
  })
  grid.text(paste(dataset_id, min_iter_feature_presence),
            x = 0.6, y = 0.95)
  dev.off()
  
  
  ####write best features to file
  output_dir <- "Data/selected_features/"
  file_name <- "best_FSM_common_features.csv"

  selected_features_df <- data.frame(matrix(nrow = 0, ncol = 4))
  for(fsm in names(selected_features)){
    print(fsm)
    row <- c(dataset_id,
             fsm,
             min_iter_feature_presence,
             paste(selected_features[[fsm]], collapse = "|")
             )
    selected_features_df <- rbind(selected_features_df, row)
  }
  colnames(selected_features_df) <- c("datasetid",
                                      "fsm",
                                      "miniter_feature_presence",
                                      "biomarkers")
  if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
  }
  
  if(length(best_fsm_vec) == dim(selected_features_df)[1]){
    file_path <- paste0(output_dir, file_name)
    write.table(selected_features_df, 
                file_path,
                sep = ",",
                row.names = FALSE, append = TRUE,
                col.names = !file.exists(file_path))    
  } else{
    print("Not writing selected features to file")
    print("all FSMs dont have common features across iter")
  }
}




get_features_from_df <- function(features_df, fsm){
  features <- features_df[features_df$fsm == fsm, "biomarkers"]
  if(length(features) != 0){
    features <- strsplit(features, split = "|", fixed = TRUE)[[1]]  
  } else{
    features <- c()
  }
  features
}

write_subset_file <- function(data, features, subset_file_path){
  data_sub <- data[gsub(".", "-", features, fixed = TRUE),]
  print(dim(data_sub))
  print(sum(is.na(data_sub)))
  print(subset_file_path)
  write.csv(data_sub, subset_file_path)
}

# dparg_id = 42
# min_iter_feature_presence = 28
# subset_creation_criteria <- list("i"= c("ranger_impu_cor"))
# subset_file_name_substr = "ranger"
# 
# create_all_common = TRUE


dparg_id = 13
dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic
min_iter_feature_presence = 28
subset_creation_criteria <- list("i"= c("ga_rf"))
subset_file_name_substr = "ga_rf"
create_all_common = FALSE



create_data_subsets <- function(dparg_id, 
                                min_iter_feature_presence,
                                subset_creation_criteria,
                                subset_file_name_substr = "common3",
                                create_all_common = TRUE,
                                dataset_pipeline_arguments = dataset_pipeline_arguments,
                                data_file_path = NA){
  
  ds <- dataset_pipeline_arguments[[dparg_id]]
  dataset_id <- paste(ds$dataset_id, ds$classification_criteria, sep = "_")
  print(dataset_id)
  
  output_dir <- "Data/selected_features/"
  file_name <- "best_FSM_common_features.csv"
  features_df <- read.csv(paste0(output_dir, file_name))
  
  features_df <- features_df %>%
    filter(datasetid == dataset_id) %>%
    filter(miniter_feature_presence == min_iter_feature_presence)
  
  
  ######create new data subsets
  
  if(grepl(pattern = "transcript", x = dataset_id, fixed = TRUE)){
    if(is.na(data_file_path)){
      data_file_path <- "Data/RNA/umi_counts_initial_cohort.csv"
    }
    data <- read.table(data_file_path, header=TRUE, sep=",", row.names=1, skip=0,
                       nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")   
    output_dir <- "Data/RNA/subset_initial_cohort/"
  } else if(grepl(pattern = "prot", x = dataset_id, fixed = TRUE)){
    if(is.na(data_file_path)){
      data_file_path <- "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil.csv"
    }
    data <- read.table(data_file_path, header=TRUE, sep=",", row.names=1, skip=0,
                         nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")  
    output_dir <- "Data/Protein/subset_initial_cohort/"
  } else{
    print("Unknown dataset_id type !")
    return
  }
  if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
  }

  best_features_df <- data.frame(matrix(nrow = 0, ncol = 5))
  
  #create datasubset with features common in all FSMs
  if(create_all_common == TRUE){
    intersect_list <- list()
    i <- 1
    for(i_fsm in features_df$fsm){
      intersect_list[[i]] <- get_features_from_df(features_df, i_fsm)
      i <- i + 1
    }
    features <- Reduce(intersect, intersect_list)  
    write_subset_file(data, features, 
                      subset_file_path = paste0(output_dir, dataset_id, "_all_common_",
                                        min_iter_feature_presence, ".csv"))
  
    best_features_df <- 
      rbind(
        best_features_df,
        data.frame(
          dataset_id = dataset_id,
          description = "all_common",
          min_iter_feature_presence = min_iter_feature_presence,
          biomarkers = paste(features, collapse = "|"),
          size = length(features)
        )        
      )
    
    if(sum(grepl("piR", features, fixed = TRUE)) > 0){
      print("all common features contain piRNA !")
      features_without_pir <- features[!grepl("piR", features, fixed = TRUE)]
      write_subset_file(data, features_without_pir, 
                        subset_file_path = paste0(output_dir, dataset_id, 
                                                  "_all_common_nopir_",
                                                  min_iter_feature_presence, ".csv"))
      best_features_df <- 
        rbind(
          best_features_df,
          data.frame(
            dataset_id = dataset_id,
            description = "all_common_nopir",
            min_iter_feature_presence = min_iter_feature_presence,
            biomarkers = paste(features_without_pir, collapse = "|"),
            size = length(features_without_pir)
          )        
        )
    }
  }
  
  if(length(subset_creation_criteria) != 0){
    # example subset_creation_criteria
    # subset_creation_criteria <- list("i"= c("t-test",
    #                                         "wilcoxontest",
    #                                         "wilcoxontest_pval_0.005"),
    #                                  "d"= c("mrmr75"))
    intersect_list <- list()
    i <- 1
    for(i_fsm in subset_creation_criteria[["i"]]){
      print(i_fsm)
      intersect_list[[i]] <- get_features_from_df(features_df, i_fsm)
      i <- i + 1
    }

    features <- Reduce(intersect, intersect_list)
    #currently handles case when "d" has single value only
    if(length(subset_creation_criteria[["d"]]) == 1){
      features <- setdiff(features,
                          get_features_from_df(features_df, subset_creation_criteria[["d"]])
                          )
    }
    
    write_subset_file(data, features, 
                      subset_file_path = paste0(output_dir, dataset_id, 
                                               "_", subset_file_name_substr, "_",
                                               min_iter_feature_presence, ".csv"))
    best_features_df <- 
      rbind(
        best_features_df,
        data.frame(
          dataset_id = dataset_id,
          description = subset_file_name_substr,
          min_iter_feature_presence = min_iter_feature_presence,
          biomarkers = paste(features, collapse = "|"),
          size = length(features)
        )        
      )
    
    if(sum(grepl("piR", features, fixed = TRUE)) > 0){
      print("subset creation criteria features contain piRNA !")
      features_without_pir <- features[!grepl("piR", features, fixed = TRUE)]
      write_subset_file(data, features_without_pir, 
                        subset_file_path = paste0(output_dir, dataset_id, 
                                                  "_", subset_file_name_substr, 
                                                  "_nopir_",
                                                  min_iter_feature_presence, ".csv"))
      best_features_df <- 
        rbind(
          best_features_df,
          data.frame(
            dataset_id = dataset_id,
            description = paste0(subset_file_name_substr, "_nopir"),
            min_iter_feature_presence = min_iter_feature_presence,
            biomarkers = paste(features_without_pir, collapse = "|"),
            size = length(features_without_pir)
          )        
        )
    }
    
  }
  
  output_dir <- "Data/selected_features/"
  file_name <- "best_features.csv"
  if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
  }
  file_path <- paste0(output_dir, file_name)
  write.table(best_features_df, 
              file_path,
              sep = ",",
              row.names = FALSE, append = TRUE,
              col.names = !file.exists(file_path))

}

#p PREOPEVsMET
explore_common_features(dparg_id = 41,
                        best_fsm_vec = c("t-test",
                                         "wilcoxontest",
                                         "mrmr75",
                                         "ranger_impu_cor"),
                        min_iter_feature_presence = 28
)
explore_common_features(dparg_id = 41,
                        best_fsm_vec = c("t-test",
                                         "wilcoxontest",
                                         "mrmr75",
                                         "ranger_impu_cor"),
                        min_iter_feature_presence = 29
)
create_data_subsets(dparg_id = 41,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria = list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75")
create_data_subsets(dparg_id = 41,
                    min_iter_feature_presence = 29,
                    subset_creation_criteria = list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)

#p PREOPEVsHC 
explore_common_features(dparg_id = 42,
                        best_fsm_vec = c("t-test",
                                         "wilcoxontest",
                                         "mrmr50",
                                         "ranger_impu_cor"),
                        min_iter_feature_presence = 28
)
explore_common_features(dparg_id = 42,
                        best_fsm_vec = c("t-test",
                                         "wilcoxontest",
                                         "mrmr50",
                                         "ranger_impu_cor"),
                        min_iter_feature_presence = 29
)

create_data_subsets(dparg_id = 42,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("ranger_impu_cor")),
                    subset_file_name_substr = "ranger")

create_data_subsets(dparg_id = 42,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr50")),
                    subset_file_name_substr = "mrmr50",
                    create_all_common = FALSE) 

create_data_subsets(dparg_id = 42,
                    min_iter_feature_presence = 29,
                    subset_creation_criteria <- list("i"= c("ranger_impu_cor")),
                    subset_file_name_substr = "ranger")

create_data_subsets(dparg_id = 42,
                    min_iter_feature_presence = 29,
                    subset_creation_criteria <- list("i"= c("mrmr50")),
                    subset_file_name_substr = "mrmr50",
                    create_all_common = FALSE) 

#p METVsHC
explore_common_features(dparg_id = 43,
                        best_fsm_vec = c("t-test",
                                         "wilcoxontest",
                                         "mrmr100",
                                         "ranger_impu_cor"),
                        min_iter_feature_presence = 28
)
explore_common_features(dparg_id = 43,
                        best_fsm_vec = c("t-test",
                                         "wilcoxontest",
                                         "mrmr100"),
                        min_iter_feature_presence = 29
)

create_data_subsets(dparg_id = 43,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("ranger_impu_cor")),
                    subset_file_name_substr = "ranger",
                    create_all_common = FALSE) 
create_data_subsets(dparg_id = 43,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE) 
create_data_subsets(dparg_id = 43,
                    min_iter_feature_presence = 29,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE) 


#p PREOPEVsPOSTOPE_T 
explore_common_features(dparg_id = 125,
                        best_fsm_vec = c("wilcoxontest",
                                         "mrmr100"),
                        min_iter_feature_presence = 28
)

create_data_subsets(dparg_id = 125,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria = list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 125,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria = list("i"= c("wilcoxontest",
                                                           "mrmr100")),
                    subset_file_name_substr = "common2",
                    create_all_common = FALSE)



#p PREOPEVsPOSTOPE_P 
explore_common_features(dparg_id = 126,
                        best_fsm_vec = c("t-test",
                                         "wilcoxontest",
                                         "mrmr100"),
                        min_iter_feature_presence = 28
)
explore_common_features(dparg_id = 126,
                        best_fsm_vec = c("t-test",
                                         "wilcoxontest",
                                         "mrmr100"),
                        min_iter_feature_presence = 29
)

create_data_subsets(dparg_id = 126,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("t-test",
                                                            "wilcoxontest"),
                                                     "d"= c("mrmr100")),
                    subset_file_name_substr = "common2only")
create_data_subsets(dparg_id = 126,
                    min_iter_feature_presence = 29,
                    subset_creation_criteria <- list("i"= c("t-test",
                                                            "wilcoxontest")),
                    subset_file_name_substr = "common2")

#p POSTOPE_TVsPOSTOPE_P
explore_common_features(dparg_id = 127,
                        best_fsm_vec = c("t-test",
                                         "wilcoxontest",
                                         "ranger_impu_cor",
                                         "mrmr100"),
                        min_iter_feature_presence = 28
)
explore_common_features(dparg_id = 127,
                        best_fsm_vec = c("t-test",
                                         "wilcoxontest",
                                         "ranger_impu_cor",
                                         "mrmr100"),
                        min_iter_feature_presence = 29
)

create_data_subsets(dparg_id = 127,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("ranger_impu_cor")),
                    subset_file_name_substr = "ranger",
                    create_all_common = FALSE) 
create_data_subsets(dparg_id = 127,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE) 

#p POSTOPE_TVsREC_T
explore_common_features(dparg_id = 131,
                        best_fsm_vec = c("t-test",
                                         "wilcoxontest",
                                         "mrmr20"),
                        min_iter_feature_presence = 28
)

explore_common_features(dparg_id = 131,
                        best_fsm_vec = c("t-test",
                                         "wilcoxontest"),
                        min_iter_feature_presence = 29
)

create_data_subsets(dparg_id = 131,
                    min_iter_feature_presence = 29,
                    subset_creation_criteria <- list("i"= c("t-test")),
                    subset_file_name_substr = "t-test")
create_data_subsets(dparg_id = 131,
                    min_iter_feature_presence = 29,
                    subset_creation_criteria <- list("i"= c("wilcoxontest")),
                    subset_file_name_substr = "wilcoxontest",
                    create_all_common = FALSE) 

#p POSTOPE_PVsREC_P
explore_common_features(dparg_id = 133,
                        best_fsm_vec = c("t-test",
                                         "wilcoxontest",
                                         "mrmr50",
                                         "ranger_impu_cor"),
                        min_iter_feature_presence = 28
)
explore_common_features(dparg_id = 133,
                        best_fsm_vec = c("t-test",
                                         "wilcoxontest",
                                         "mrmr50",
                                         "ranger_impu_cor"),
                        min_iter_feature_presence = 29
)

create_data_subsets(dparg_id = 133,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr50")),
                    subset_file_name_substr = "mrmr50",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 133,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("ranger_impu_cor")),
                    subset_file_name_substr = "ranger",
                    create_all_common = FALSE)

create_data_subsets(dparg_id = 133,
                    min_iter_feature_presence = 29,
                    subset_creation_criteria <- list("i"= c("mrmr50")),
                    subset_file_name_substr = "mrmr50",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 133,
                    min_iter_feature_presence = 29,
                    subset_creation_criteria <- list("i"= c("ranger_impu_cor")),
                    subset_file_name_substr = "ranger",
                    create_all_common = FALSE)

#p POSTOPE_TVsPREREC
explore_common_features(dparg_id = 135,
                        best_fsm_vec = c("ga_rf",
                                         "wilcoxontest",
                                         "mrmr100"),
                        min_iter_feature_presence = 28
)
create_data_subsets(dparg_id = 135,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("ga_rf")),
                    subset_file_name_substr = "ga_rf",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 135,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("wilcoxontest",
                                                            "mrmr100")),
                    subset_file_name_substr = "common2",
                    create_all_common = FALSE)


#p PREOPEVsREC_TP
explore_common_features(dparg_id = 137,
                        best_fsm_vec = c("t-test",
                                         "wilcoxontest",
                                         "mrmr100"),
                        min_iter_feature_presence = 28
)
explore_common_features(dparg_id = 137,
                        best_fsm_vec = c("t-test",
                                         "wilcoxontest",
                                         "mrmr100"),
                        min_iter_feature_presence = 29
)

create_data_subsets(dparg_id = 137,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("t-test",
                                                            "wilcoxontest")),
                    subset_file_name_substr = "common2")
create_data_subsets(dparg_id = 137,
                    min_iter_feature_presence = 29,
                    subset_creation_criteria <- list("i"= c("t-test",
                                                            "wilcoxontest")),
                    subset_file_name_substr = "common2")




#t PREOPEVsMET
explore_common_features(dparg_id = 139,
                        best_fsm_vec = c("t-test_pval_0.025",
                                         "mrmr100",
                                         "wilcoxontest_BH",
                                         "ranger_impu_cor"),
                        min_iter_feature_presence = 28
)
explore_common_features(dparg_id = 139,
                        best_fsm_vec = c("t-test_pval_0.025",
                                         "mrmr100",
                                         "wilcoxontest_BH",
                                         "ranger_impu_cor"),
                        min_iter_feature_presence = 29
)

create_data_subsets(dparg_id = 139,
                    min_iter_feature_presence = 29,
                    subset_creation_criteria <- list("i"= c("t-test_pval_0.025",
                                                            "mrmr100",
                                                            "ranger_impu_cor")),
                    subset_file_name_substr = "common3")

#t PREOPEVsHC 
explore_common_features(dparg_id = 140,
                        best_fsm_vec = c("t-test_pval_0.025",
                                         "wilcoxontest",
                                         "mrmr30",
                                         "ranger_impu_cor"),
                        min_iter_feature_presence = 28
)
explore_common_features(dparg_id = 140,
                        best_fsm_vec = c("t-test_pval_0.025",
                                         "wilcoxontest",
                                         "mrmr30",
                                         "ranger_impu_cor"),
                        min_iter_feature_presence = 29
)
create_data_subsets(dparg_id = 140,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria = list("i"= c("mrmr30")),
                    subset_file_name_substr = "mrmr30")
create_data_subsets(dparg_id = 140,
                    min_iter_feature_presence = 29,
                    subset_creation_criteria = list("i"= c("mrmr30")),
                    subset_file_name_substr = "mrmr30")

#t METVsHC
explore_common_features(dparg_id = 141,
                        best_fsm_vec = c("ranger_impu_cor",
                                         "wilcoxontest",
                                         "RF_RFE",
                                         "mrmr100"),
                        min_iter_feature_presence = 28
)
explore_common_features(dparg_id = 141,
                        best_fsm_vec = c("ranger_impu_cor",
                                         "wilcoxontest",
                                         "RF_RFE",
                                         "mrmr100"),
                        min_iter_feature_presence = 29
)
create_data_subsets(dparg_id = 141,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("ranger_impu_cor",
                                                            "wilcoxontest",
                                                            "mrmr100")),
                    subset_file_name_substr = "common3")
create_data_subsets(dparg_id = 141,
                    min_iter_feature_presence = 29,
                    subset_creation_criteria <- list("i"= c("ranger_impu_cor",
                                                            "wilcoxontest",
                                                            "mrmr100")),
                    subset_file_name_substr = "common3")


#t PREOPEVsPOSTOPE_T 
explore_common_features(dparg_id = 145,
                        best_fsm_vec = c("t-test",
                                         "mrmr30"),
                        min_iter_feature_presence = 28
)
create_data_subsets(dparg_id = 145,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria = list("i"= c("t-test")),
                    subset_file_name_substr = "t-test")

#t PREOPEVsPOSTOPE_P 
explore_common_features(dparg_id = 146,
                        best_fsm_vec = c("ranger_impu_cor",
                                         "mrmr100"),
                        min_iter_feature_presence = 27
)
create_data_subsets(dparg_id = 146,
                    min_iter_feature_presence = 27,
                    subset_creation_criteria <- list("i"= c("ranger_impu_cor")),
                    subset_file_name_substr = "ranger",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 146,
                    min_iter_feature_presence = 27,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE)


#t POSTOPE_TVsPOSTOPE_P
explore_common_features(dparg_id = 147,
                        best_fsm_vec = c("ranger_impu_cor",
                                         "mrmr100"),
                        min_iter_feature_presence = 28
)
explore_common_features(dparg_id = 147,
                        best_fsm_vec = c("ranger_impu_cor",
                                         "mrmr100"),
                        min_iter_feature_presence = 29
)
create_data_subsets(dparg_id = 147,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("ranger_impu_cor")),
                    subset_file_name_substr = "ranger",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 147,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE)


#t POSTOPE_TVsREC_T
explore_common_features(dparg_id = 151,
                        best_fsm_vec = c("ranger_impu_cor",
                                         "mrmr100",
                                         "t-test_pval_0.025",
                                         "wilcoxontest"),
                        min_iter_feature_presence = 28
)
explore_common_features(dparg_id = 151,
                        best_fsm_vec = c("ranger_impu_cor",
                                         "mrmr100",
                                         "t-test_pval_0.025",
                                         "wilcoxontest"),
                        min_iter_feature_presence = 29
)
create_data_subsets(dparg_id = 151,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("ranger_impu_cor")),
                    subset_file_name_substr = "ranger",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 151,
                    min_iter_feature_presence = 29,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE)

#t POSTOPE_PVsREC_P
explore_common_features(dparg_id = 153,
                        best_fsm_vec = c("ranger_impu_cor", "mrmr100"),
                        min_iter_feature_presence = 27
)
create_data_subsets(dparg_id = 153,
                    min_iter_feature_presence = 27,
                    subset_creation_criteria <- list("i"= c("ranger_impu_cor")),
                    subset_file_name_substr = "ranger",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 153,
                    min_iter_feature_presence = 27,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE)

#t POSTOPE_TVsPREREC
explore_common_features(dparg_id = 155,
                        best_fsm_vec = c("t-test",
                                         "mrmr20",
                                         "wilcoxontest"),
                        min_iter_feature_presence = 28
)
create_data_subsets(dparg_id = 155,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("wilcoxontest")),
                    subset_file_name_substr = "wilcoxontest",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 155,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("t-test")),
                    subset_file_name_substr = "t-test",
                    create_all_common = FALSE)

#t PREOPEVsREC_TP
explore_common_features(dparg_id = 157,
                        best_fsm_vec = c("ranger_impu_cor",
                                         "mrmr100",
                                         "t-test"),
                        min_iter_feature_presence = 28
)
create_data_subsets(dparg_id = 157,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("ranger_impu_cor")),
                    subset_file_name_substr = "ranger",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 157,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE)


####new comparisons

#t PREOPEVsPOSTOPE_TP
explore_common_features(dparg_id = 229,
                        best_fsm_vec = c("t-test",
                                         "ranger_impu_cor",
                                         "mrmr100"),
                        min_iter_feature_presence = 28
)

create_data_subsets(dparg_id = 229,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("t-test")),
                    subset_file_name_substr = "t-test",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 229,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("ranger_impu_cor")),
                    subset_file_name_substr = "ranger",
                    create_all_common = FALSE)


#t POSTOPE_TPVsREC_TP
explore_common_features(dparg_id = 231,
                        best_fsm_vec = c("ranger_impu_cor",
                                         "mrmr100",
                                         "t-test_pval_0.01"),
                        min_iter_feature_presence = 28
)
create_data_subsets(dparg_id = 231,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 231,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("ranger_impu_cor",
                                                            "mrmr100")),
                    subset_file_name_substr = "r_m100",
                    create_all_common = FALSE)

#p PREOPEVsPOSTOPE_TP
explore_common_features(dparg_id = 239,
                        best_fsm_vec = c("ga_rf",
                                         "t-test",
                                         "wilcoxontest_pval_0.025",
                                         "mrmr100"),
                        min_iter_feature_presence = 28
)
create_data_subsets(dparg_id = 239,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 239,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("wilcoxontest_pval_0.025")),
                    subset_file_name_substr = "w_025",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 239,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("ga_rf")),
                    subset_file_name_substr = "ga_rf",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 239,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100",
                                                            "wilcoxontest_pval_0.025",
                                                            "t-test")),
                    subset_file_name_substr = "common3",
                    create_all_common = FALSE)

#p POSTOPE_TPVsREC_TP

explore_common_features(dparg_id = 241,
                        best_fsm_vec = c("t-test",
                                         "wilcoxontest"),
                        min_iter_feature_presence = 29
)

create_data_subsets(dparg_id = 241,
                    min_iter_feature_presence = 29,
                    subset_creation_criteria <- list("i"= c("t-test")),
                    subset_file_name_substr = "t-test",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 241,
                    min_iter_feature_presence = 29,
                    subset_creation_criteria <- list("i"= c("wilcoxontest")),
                    subset_file_name_substr = "wilcoxontest",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 241,
                    min_iter_feature_presence = 29,
                    subset_creation_criteria <- list("i"= c("t-test", "wilcoxontest")),
                    subset_file_name_substr = "t_w",
                    create_all_common = FALSE)


##############common features across all

best_features <- read.csv("Data/selected_features/best_features_with_add_col.csv")
best_features <- best_features %>%
  filter(is_best == 1)

best_features_pr <- best_features %>%
  filter(grepl("proteomic", dataset_id, fixed = TRUE))
i <- 1
intersect_list <- list()
for(ds in best_features_pr$dataset_id){
  print(ds)
  features <- best_features_pr[best_features_pr$dataset_id == ds, "biomarkers"]
  if(length(features) != 0){
    features <- strsplit(features, split = "|", fixed = TRUE)[[1]]  
  } else{
    features <- c()
  }
  intersect_list[[i]] <- features
  i <- i + 1
}
features <- Reduce(intersect, intersect_list)
names(intersect_list) <- best_features_pr$dataset_id
upset(fromList(intersect_list), set_size.show = TRUE)  

best_features_tr <- best_features %>%
  filter(grepl("transcriptomic", dataset_id, fixed = TRUE))
i <- 1
intersect_list <- list()
for(ds in best_features_tr$dataset_id){
  print(ds)
  features <- best_features_tr[best_features_tr$dataset_id == ds, "biomarkers"]
  if(length(features) != 0){
    features <- strsplit(features, split = "|", fixed = TRUE)[[1]]  
  } else{
    features <- c()
  }
  intersect_list[[i]] <- features
  i <- i + 1
}
features <- Reduce(intersect, intersect_list)
names(intersect_list) <- best_features_tr$dataset_id
upset(fromList(intersect_list), set_size.show = TRUE)  
dev.off()


#####################################

#initial cohort with new quantified results


dparg_id = 1
dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic
best_fsm_vec = c("t-test")
min_iter_feature_presence = 28
results_dir = "fem_pipeline_results_tr"
dir_path = "plots/FEMPipeline_new_quant/common_features_upset"


#t PREOPEVsPOSTOPE_TP
explore_common_features(dparg_id = 1,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                        best_fsm_vec = c("t-test", "mrmr10", 
                                         "mrmr75", "mrmr_perc50"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_tr",
                        dir_path = "plots/FEMPipeline_new_quant/common_features_upset")



create_data_subsets(dparg_id = 1,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)



#t POSTOPE_TPVSREC_TP
explore_common_features(dparg_id = 5,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                        best_fsm_vec = c("t-test", "mrmr_perc50", 
                                         "mrmr75"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_tr",
                        dir_path = "plots/FEMPipeline_new_quant/common_features_upset")

create_data_subsets(dparg_id = 5,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)


#t PREOPEVSREC_TP
explore_common_features(dparg_id = 9,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                        best_fsm_vec = c("mrmr30", "mrmr_perc50", 
                                         "mrmr75"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_tr",
                        dir_path = "plots/FEMPipeline_new_quant/common_features_upset")

create_data_subsets(dparg_id = 9,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)


#####################################
#proteomic with no norm

explore_common_features(dparg_id = 1,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                        best_fsm_vec = c("ga_rf", "mrmr_perc50", "mrmr75", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_pr",
                        dir_path = "plots/FEMPipeline_prot_no_norm/common_features_upset")

explore_common_features(dparg_id = 5,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                        best_fsm_vec = c("RF_RFE", "ranger_pos_impu_cor", "ga_rf", 
                                         "mrmr_perc50", "mrmr50", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_pr",
                        dir_path = "plots/FEMPipeline_prot_no_norm/common_features_upset")

explore_common_features(dparg_id = 9,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                        best_fsm_vec = c("mrmr30", "mrmr50", "mrmr75", 
                                         "wilcoxontest", "t-test"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_pr",
                        dir_path = "plots/FEMPipeline_prot_no_norm/common_features_upset")


create_data_subsets(dparg_id = 1,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 1,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)

create_data_subsets(dparg_id = 5,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 5,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE)

create_data_subsets(dparg_id = 9,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 9,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("t-test")),
                    subset_file_name_substr = "t-test",
                    create_all_common = FALSE)

#####################################
#proteomic with quantile train param norm

explore_common_features(dparg_id = 13,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                        best_fsm_vec = c("ga_rf", "mrmr30", "mrmr75", "mrmr100", "mrmr50"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_pr",
                        dir_path = "plots/FEMPipeline_prot_quantile_norm_with_train_param/common_features_upset")

explore_common_features(dparg_id = 17,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                        best_fsm_vec = c("mrmr_perc50", "ranger_pos_impu_cor", "ga_rf", 
                                         "RF_RFE", "t-test"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_pr",
                        dir_path = "plots/FEMPipeline_prot_quantile_norm_with_train_param/common_features_upset")

explore_common_features(dparg_id = 21,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                        best_fsm_vec = c("mrmr30", "t-test", "ga_rf", "RF_RFE", 
                                         "wilcoxontest"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_pr",
                        dir_path = "plots/FEMPipeline_prot_quantile_norm_with_train_param/common_features_upset")

create_data_subsets(dparg_id = 13,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("ga_rf")),
                    subset_file_name_substr = "ga_rf",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 13,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE)

create_data_subsets(dparg_id = 17,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 17,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("t-test")),
                    subset_file_name_substr = "t-test",
                    create_all_common = FALSE)

create_data_subsets(dparg_id = 21,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("t-test")),
                    subset_file_name_substr = "t-test",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 21,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("wilcoxontest")),
                    subset_file_name_substr = "wilcoxontest",
                    create_all_common = FALSE)


#####################################
#proteomic common with no norm

explore_common_features(dparg_id = 41,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                        best_fsm_vec = c("mrmr75", "mrmr100", "mrmr50", "RF_RFE"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_pr_common",
                        dir_path = "plots/FEMPipeline_prot_common_no_norm/common_features_upset")

explore_common_features(dparg_id = 37,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                        best_fsm_vec = c("mrmr_perc50", "ga_rf", "t-test", "mrmr30", "wilcoxontest"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_pr_common",
                        dir_path = "plots/FEMPipeline_prot_common_no_norm/common_features_upset")

explore_common_features(dparg_id = 45,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                        best_fsm_vec = c("mrmr30", "mrmr50", "wilcoxontest", "t-test", "mrmr75"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_pr_common",
                        dir_path = "plots/FEMPipeline_prot_common_no_norm/common_features_upset")


create_data_subsets(dparg_id = 41,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE, 
                    data_file_path = "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil_common.csv")
create_data_subsets(dparg_id = 37,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE, 
                    data_file_path = "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil_common.csv")
create_data_subsets(dparg_id = 37,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("t-test")),
                    subset_file_name_substr = "t-test",
                    create_all_common = FALSE, 
                    data_file_path = "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil_common.csv")
create_data_subsets(dparg_id = 45,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr50")),
                    subset_file_name_substr = "mrmr50",
                    create_all_common = FALSE, 
                    data_file_path = "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil_common.csv")

#####################################
#proteomic common with quantile train param norm

explore_common_features(dparg_id = 53,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                        best_fsm_vec = c("ranger_pos_impu_cor", "mrmr_perc50", "RF_RFE", "mrmr50",
                                         "t-test", "wilcoxontest", "mrmr75", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_pr_common",
                        dir_path = "plots/FEMPipeline_prot_common_quantile_norm_with_train_param/common_features_upset")

explore_common_features(dparg_id = 49,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                        best_fsm_vec = c("ga_rf", "mrmr75", "mrmr50", "mrmr30"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_pr_common",
                        dir_path = "plots/FEMPipeline_prot_common_quantile_norm_with_train_param/common_features_upset")

explore_common_features(dparg_id = 57,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                        best_fsm_vec = c("mrmr100", "mrmr75", "mrmr50", "RF_RFE"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_pr_common",
                        dir_path = "plots/FEMPipeline_prot_common_quantile_norm_with_train_param/common_features_upset")


create_data_subsets(dparg_id = 53,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE, 
                    data_file_path = "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil_common.csv")
create_data_subsets(dparg_id = 53,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("t-test")),
                    subset_file_name_substr = "t-test",
                    create_all_common = FALSE, 
                    data_file_path = "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil_common.csv")
create_data_subsets(dparg_id = 49,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE, 
                    data_file_path = "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil_common.csv")
create_data_subsets(dparg_id = 57,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE, 
                    data_file_path = "Data/Protein/formatted_data/Q1-6_nonorm_formatted_impute50fil_common.csv")



#####################################
#proteomic common combat


explore_common_features(dparg_id = 69,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                        best_fsm_vec = c("ranger_pos_impu_cor", "t-test", "wilcoxontest", 
                                         "mrmr50", "mrmr75"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_pr_common_combat",
                        dir_path = "plots/FEMPipeline_prot_common_combat/common_features_upset")
explore_common_features(dparg_id = 73,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                        best_fsm_vec = c("ga_rf", "mrmr100", "wilcoxontest",
                                         "RF_RFE", "mrmr75"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_pr_common_combat",
                        dir_path = "plots/FEMPipeline_prot_common_combat/common_features_upset")
explore_common_features(dparg_id = 77,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                        best_fsm_vec = c("mrmr75", "mrmr100", "t-test"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_pr_common_combat",
                        dir_path = "plots/FEMPipeline_prot_common_combat/common_features_upset")

create_data_subsets(dparg_id = 69,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("t-test")),
                    subset_file_name_substr = "t-test",
                    create_all_common = FALSE, 
                    data_file_path = "Data/Protein/initial_data.combat.POSTOPE_TPVsREC_TP.csv")
create_data_subsets(dparg_id = 73,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("ga_rf")),
                    subset_file_name_substr = "ga_rf",
                    create_all_common = FALSE, 
                    data_file_path = "Data/Protein/initial_data.combat.PREOPEVsPOSTOPE_TP.csv")
create_data_subsets(dparg_id = 73,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE, 
                    data_file_path = "Data/Protein/initial_data.combat.PREOPEVsPOSTOPE_TP.csv")
create_data_subsets(dparg_id = 77,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE, 
                    data_file_path = "Data/Protein/initial_data.combat.PREOPEVsREC_TP.csv")


#####################################
#transcriptomic common combat


explore_common_features(dparg_id = 19,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                        best_fsm_vec = c("t-test", "mrmr75", "mrmr100", 
                                         "mrmr_perc50", "mrmr50"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_tr_common_combat",
                        dir_path = "plots/FEMPipeline_tr_common_combat/common_features_upset")
explore_common_features(dparg_id = 19,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                        best_fsm_vec = c("t-test", "mrmr75", "mrmr100", 
                                         "mrmr_perc50", "mrmr50"),
                        min_iter_feature_presence = 29,
                        results_dir = "fem_pipeline_results_tr_common_combat",
                        dir_path = "plots/FEMPipeline_tr_common_combat/common_features_upset")
explore_common_features(dparg_id = 23,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                        best_fsm_vec = c("ga_rf", "RF_RFE", "mrmr75",
                                         "mrmr30", "mrmr50", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_tr_common_combat",
                        dir_path = "plots/FEMPipeline_tr_common_combat/common_features_upset")
explore_common_features(dparg_id = 23,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                        best_fsm_vec = c("ga_rf", "RF_RFE", "mrmr75",
                                         "mrmr30", "mrmr50", "mrmr100"),
                        min_iter_feature_presence = 29,
                        results_dir = "fem_pipeline_results_tr_common_combat",
                        dir_path = "plots/FEMPipeline_tr_common_combat/common_features_upset")
explore_common_features(dparg_id = 27,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                        best_fsm_vec = c("ga_rf", "ranger_pos_impu_cor", "mrmr100",
                                         "mrmr50", "mrmr75", "mrmr_perc50"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_tr_common_combat",
                        dir_path = "plots/FEMPipeline_tr_common_combat/common_features_upset")
explore_common_features(dparg_id = 27,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                        best_fsm_vec = c("ga_rf", "ranger_pos_impu_cor", "mrmr100",
                                         "mrmr50", "mrmr75", "mrmr_perc50"),
                        min_iter_feature_presence = 29,
                        results_dir = "fem_pipeline_results_tr_common_combat",
                        dir_path = "plots/FEMPipeline_tr_common_combat/common_features_upset")
explore_common_features(dparg_id = 27,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                        best_fsm_vec = c("ga_rf", "ranger_pos_impu_cor", "mrmr100",
                                         "mrmr50", "mrmr75", "mrmr_perc50"),
                        min_iter_feature_presence = 30,
                        results_dir = "fem_pipeline_results_tr_common_combat",
                        dir_path = "plots/FEMPipeline_tr_common_combat/common_features_upset")


create_data_subsets(dparg_id = 19,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                    min_iter_feature_presence = 29,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE, 
                    data_file_path = "Data/RNA/initial_data.combat.POSTOPE_TPVsREC_TP.csv")
create_data_subsets(dparg_id = 23,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                    min_iter_feature_presence = 29,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE, 
                    data_file_path = "Data/RNA/initial_data.combat.PREOPEVsPOSTOPE_TP.csv")
create_data_subsets(dparg_id = 27,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("ranger_pos_impu_cor")),
                    subset_file_name_substr = "ranger_pos_impu_cor",
                    create_all_common = FALSE, 
                    data_file_path = "Data/RNA/initial_data.combat.PREOPEVsREC_TP.csv")
create_data_subsets(dparg_id = 27,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                    min_iter_feature_presence = 30,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE, 
                    data_file_path = "Data/RNA/initial_data.combat.PREOPEVsREC_TP.csv")



#####################################
#transcriptomic common 

explore_common_features(dparg_id = 31,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                        best_fsm_vec = c("mrmr_perc50", "mrmr50", "mrmr30",
                                         "ranger_pos_impu_cor", "ga_rf", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_tr_common",
                        dir_path = "plots/FEMPipeline_tr_common/common_features_upset")
explore_common_features(dparg_id = 35,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                        best_fsm_vec = c("mrmr_perc50", "mrmr100", "ga_rf",
                                         "mrmr30", "mrmr50"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_tr_common",
                        dir_path = "plots/FEMPipeline_tr_common/common_features_upset")
explore_common_features(dparg_id = 35,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                        best_fsm_vec = c("mrmr_perc50", "mrmr100", "ga_rf",
                                         "mrmr30", "mrmr50"),
                        min_iter_feature_presence = 29,
                        results_dir = "fem_pipeline_results_tr_common",
                        dir_path = "plots/FEMPipeline_tr_common/common_features_upset")
explore_common_features(dparg_id = 35,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                        best_fsm_vec = c("mrmr_perc50", "mrmr100", "ga_rf",
                                         "mrmr30", "mrmr50"),
                        min_iter_feature_presence = 30,
                        results_dir = "fem_pipeline_results_tr_common",
                        dir_path = "plots/FEMPipeline_tr_common/common_features_upset")
explore_common_features(dparg_id = 39,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                        best_fsm_vec = c("mrmr75", "mrmr100", "ga_rf", "ranger_pos_impu_cor"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_tr_common",
                        dir_path = "plots/FEMPipeline_tr_common/common_features_upset")
explore_common_features(dparg_id = 39,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                        best_fsm_vec = c("mrmr75", "mrmr100", "ga_rf", "ranger_pos_impu_cor"),
                        min_iter_feature_presence = 30,
                        results_dir = "fem_pipeline_results_tr_common",
                        dir_path = "plots/FEMPipeline_tr_common/common_features_upset")



create_data_subsets(dparg_id = 31,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE, 
                    data_file_path = "Data/RNA/umi_counts_initial_cohort_common_tr.csv")
create_data_subsets(dparg_id = 31,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr50")),
                    subset_file_name_substr = "mrmr50",
                    create_all_common = FALSE, 
                    data_file_path = "Data/RNA/umi_counts_initial_cohort_common_tr.csv")
create_data_subsets(dparg_id = 35,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                    min_iter_feature_presence = 30,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE, 
                    data_file_path = "Data/RNA/umi_counts_initial_cohort_common_tr.csv")
create_data_subsets(dparg_id = 35,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                    min_iter_feature_presence = 30,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE, 
                    data_file_path = "Data/RNA/umi_counts_initial_cohort_common_tr.csv")
create_data_subsets(dparg_id = 39,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                    min_iter_feature_presence = 30,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE, 
                    data_file_path = "Data/RNA/umi_counts_initial_cohort_common_tr.csv")


#####################################
#proteomic combined combat
explore_common_features(dparg_id = 85,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                        best_fsm_vec = c("mrmr_perc50", "RF_RFE", "ranger_pos_impu_cor",
                                         "t-test", "wilcoxontest"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_combined_pr_common_combat",
                        dir_path = "plots/fem_pipeline_combined_pr/common_features_upset")

explore_common_features(dparg_id = 89,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                        best_fsm_vec = c("ga_rf", "mrmr30", "mrmr50", "mrmr75", "RF_RFE"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_combined_pr_common_combat",
                        dir_path = "plots/fem_pipeline_combined_pr/common_features_upset")

explore_common_features(dparg_id = 93,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                        best_fsm_vec = c("ga_rf", "mrmr100", "mrmr75", "RF_RFE",
                                         "mrmr50"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_combined_pr_common_combat",
                        dir_path = "plots/fem_pipeline_combined_pr/common_features_upset")


create_data_subsets(dparg_id = 85,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("t-test")),
                    subset_file_name_substr = "t-test",
                    create_all_common = FALSE, 
                    data_file_path = "Data/Protein/combined_data.combat.POSTOPE_TPVsREC_TP.csv")
create_data_subsets(dparg_id = 89,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE, 
                    data_file_path = "Data/Protein/combined_data.combat.POSTOPE_TPVsREC_TP.csv")
create_data_subsets(dparg_id = 93,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_proteomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE, 
                    data_file_path = "Data/Protein/combined_data.combat.POSTOPE_TPVsREC_TP.csv")


#####################################
#transcriptomic combined combat

explore_common_features(dparg_id = 43,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                        best_fsm_vec = c("mrmr30", "ranger_pos_impu_cor",
                                         "mrmr10", "mrmr_perc50", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_combined_tr_common_combat",
                        dir_path = "plots/fem_pipeline_combined_tr/common_features_upset")
explore_common_features(dparg_id = 47,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                        best_fsm_vec = c("t-test", "wilcoxontest", "mrmr30",
                                         "mrmr50", "mrmr_perc50"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_combined_tr_common_combat",
                        dir_path = "plots/fem_pipeline_combined_tr/common_features_upset")
explore_common_features(dparg_id = 51,
                        dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                        best_fsm_vec = c("mrmr_perc50", "mrmr50", "mrmr30",
                                         "RF_RFE", "ga_rf"),
                        min_iter_feature_presence = 28,
                        results_dir = "fem_pipeline_results_combined_tr_common_combat",
                        dir_path = "plots/fem_pipeline_combined_tr/common_features_upset")

create_data_subsets(dparg_id = 43,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr30")),
                    subset_file_name_substr = "mrmr30",
                    create_all_common = FALSE, 
                    data_file_path = "Data/RNA/combined_data.combat.POSTOPE_TPVsREC_TP.csv")
create_data_subsets(dparg_id = 47,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr30")),
                    subset_file_name_substr = "mrmr30",
                    create_all_common = FALSE, 
                    data_file_path = "Data/RNA/combined_data.combat.PREOPEVsPOSTOPE_TP.csv")
create_data_subsets(dparg_id = 51,
                    dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE, 
                    data_file_path = "Data/RNA/combined_data.combat.PREOPEVsREC_TP.csv")
