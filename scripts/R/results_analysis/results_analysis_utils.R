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