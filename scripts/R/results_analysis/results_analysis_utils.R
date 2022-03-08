library(tidyverse)
library(UniProt.ws)

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

# best_features_file_path = "Data/selected_features/best_features_with_add_col.csv"
create_protein_biomarker_mapping <- function(best_features_file_path = "Data/selected_features/best_features_with_add_col.csv"){
  best_features <- read.csv(best_features_file_path)
  best_features <- best_features %>%
    filter(grepl("proteomic", dataset_id)) 
  all_biomarkers <- c()
  for(i in c(1:dim(best_features)[1])){
    print(i)
    biomarkers <- strsplit(best_features[i, "biomarkers"], split = "|", fixed = TRUE)[[1]]
    all_biomarkers <- c(all_biomarkers, biomarkers)
  }
  all_biomarkers <- unique(all_biomarkers) 
  up <- UniProt.ws()
  info <- AnnotationDbi::select(x = up, columns = c("ENTRY-NAME", "PROTEIN-NAMES", "GENENAME"),
                                keys = all_biomarkers)
  write.csv(info, "Data/selected_features/protein_biomarker_info.csv", row.names = FALSE)
}