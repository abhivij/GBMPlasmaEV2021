library(missForest)

append_path <- function(dirname, filename){
  if (dirname != "") {
    filename <- paste(dirname, filename, sep = "/")
  }
  return (filename)
}

impute_data <- function(data, filter_na_perc = 75, impute = TRUE){
  data <- data[, colSums(is.na(data)) < (filter_na_perc/100*nrow(data))]
  if (impute) {
    data <- missForest(data)$ximp
  }
  return (data)
}

compare_filter_na <- function(data){
  filter_perc_vec <- c(100, 90, 80, 75, 50, 25)
  info_df <- data.frame()
  row_info <- c(dim(data)[1], "-", dim(data)[2], 100)
  info_df <- rbind(info_df, row_info)
  colnames(info_df) <- c("num_samples", "NA_filter_percent", "num_cols", "col_percent")
  for(perc in filter_perc_vec){
    print(perc)
    filtered_data <- data[, colSums(is.na(data)) < (perc/100*nrow(data))]
    row_info <- c(dim(filtered_data)[1], perc, dim(filtered_data)[2], 
                  (dim(filtered_data)[2]*100/dim(data)[2]))
    info_df <- rbind(info_df, row_info)
  }
  return (info_df)
}
