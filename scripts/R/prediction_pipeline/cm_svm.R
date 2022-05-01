source("scripts/R/prediction_pipeline/compute_metrics.R")

svm_model <- function(data.train, label.train, data.test, label.test, 
                      data.test2, label.test2,
                      classes, kernel = "sigmoid", 
                      ...){
  
  # kernel = "sigmoid"
  kernel_name <- paste(toupper(substring(kernel, 1, 1)), substring(kernel, 2), sep = "")
  model_name <- paste(kernel_name, "Kernel SVM")
  #setting default value for metrics, to handle case where unable to train / execute classification model
  metrics <- c(0, 0) 
  
  try({
    model <- e1071::svm(data.train, factor(label.train$Label, levels = classes), probability = TRUE, kernel = kernel)
    
    pred <- predict(model, data.test, probability = TRUE)
    pred_prob <- data.frame(attr(pred, 'probabilities'))[classes[2]]
    
    pred2 <- predict(model, data.test2, probability = TRUE)
    pred_prob2 <- data.frame(attr(pred2, 'probabilities'))[classes[2]]
    
    result_df1 <- data.frame("TestDataClassName" = label.test$Label,
                             "TestDataExpectedClassName" = label.test$Label,
                             "Pred_prob" = pred_prob[,1],
                             "Prediction" = pred)
    
    result_df2 <- data.frame("TestDataClassName" = "PREREC",
                             "TestDataExpectedClassName" = "REC_TP", 
                             "Pred_prob" = pred_prob2[,1],
                             "Prediction" = pred2)
    
    result_df <- rbind(result_df1, result_df2) %>%
      mutate(Pred_prob = as.double(Pred_prob)) 
    
    # result_file_name <- "Data/prediction_result/transcriptomics.csv"
    # write.csv(result_df, result_file_name)  
    
    
    metrics <- compute_metrics(pred = pred, pred_prob = pred_prob, true_label = label.test$Label, classes = classes)    
    print(metrics)
    
    
    return (result_df)
    
    })
  
  
}
