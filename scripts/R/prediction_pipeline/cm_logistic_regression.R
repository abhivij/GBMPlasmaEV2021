# library(glmnet)
# source("R/metrics/compute_metrics.R")

logistic_regression <- function(data.train, label.train, data.test, label.test, 
                                classes, regularize = NA, ...){
  
  model_name <- "logistic regression"
  if (is.na(regularize)) {
    model_name <- paste("Simple", model_name)
  }
  else if(regularize == 'l1') {
    model_name <- paste("L1 Regularized", model_name)
  }
  else {
    model_name <- paste("L2 Regularized", model_name)
  }

  #setting default value for metrics, to handle case where unable to train / execute classification model
  metrics <- c(0, 0) 
  
  try({
    label.train$Label <- ifelse(label.train$Label == classes[1], 0, 1)
    label.test$Label <- ifelse(label.test$Label == classes[1], 0, 1)
    
    if (!is.na(regularize)) {
      #alpha = 1 => l1 regularization (lasso)
      #alpha = 0 => l2 regularization (ridge)
      if(regularize == 'l1') {
        alpha <- 1
      }
      else {
        alpha <- 0
      }

      model <- glmnet::cv.glmnet(as.matrix(data.train), label.train$Label, alpha = alpha, family = 'binomial', type.measure = 'mse')
      # plot(model)
      
      lambda_min <- model$lambda.min
      lambda_1se <- model$lambda.1se
      
      pred_prob <- predict(model, newx = as.matrix(data.test), s = lambda_1se, type = 'response')
      pred <- ifelse(pred_prob > 0.5, 1, 0)
    }
    else {
      model <- glm(label.train$Label ~., data = data.train, family = binomial)
      pred_prob <- predict(model, newdata = data.test, type = 'response')
      pred <- ifelse(pred_prob > 0.5, 1, 0)
    }
    metrics <- compute_metrics(pred = pred, pred_prob = pred_prob, true_label = label.test$Label, classes = c(0, 1))  
  
    label.test$Label <- ifelse(label.test$Label == 0, classes[1], classes[2])
    result_df <- data.frame("TestLabel" = label.test$Label,
                            "Pred_prob" = pred_prob[,1],
                            "Prediction" = pred[,1])
      
  })
  
  return (list(model_name, metrics))    
}


