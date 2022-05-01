source("scripts/R/prediction_pipeline/compute_metrics.R")

logistic_regression <- function(data.train, label.train, data.test, label.test, 
                                data.test2, label.test2,
                                classes, regularize = NA, 
                                ...){
  
  model_name <- "logistic regression"
  if (is.na(regularize)) {
    model_name <- paste("Simple", model_name)
  } else if(regularize == 'l1') {
    model_name <- paste("L1 Regularized", model_name)
  } else {
    model_name <- paste("L2 Regularized", model_name)
  }

  #setting default value for metrics, to handle case where unable to train / execute classification model
  metrics <- c(0, 0) 
  
  try({
    label.train$Label <- ifelse(label.train$Label == classes[1], 0, 1)
    label.test$Label <- ifelse(label.test$Label == classes[1], 0, 1)
    
    #PREREC samples expected to be shown as REC_TP
    label.test2$Label <- ifelse(label.test2$Label == "PREREC", 1, 0)
    
    if (!is.na(regularize)) {
      #alpha = 1 => l1 regularization (lasso)
      #alpha = 0 => l2 regularization (ridge)
      if(regularize == 'l1') {
        alpha <- 1
      }
      else {
        alpha <- 0
      }
      set.seed(1000)
      model <- glmnet::cv.glmnet(as.matrix(data.train), label.train$Label, alpha = alpha, family = 'binomial', type.measure = 'mse')
      # plot(model)
      
      lambda_min <- model$lambda.min
      lambda_1se <- model$lambda.1se
      
      best_acc <- -1
      best_cut_off <- 0.3
      for(cut_off in seq(0.3, 0.7, 0.01)){
        print(cut_off)
        pred_prob.train <- predict(model, newx = as.matrix(data.train), s = lambda_1se, type = 'response')
        pred.train <- ifelse(pred_prob.train > cut_off, 1, 0)
        acc <- mean(pred.train == label.train$Label)
        if(acc > best_acc){
          best_acc <- acc
          best_cut_off <- cut_off
        }
      }
      
      
      pred_prob <- predict(model, newx = as.matrix(data.test), s = lambda_1se, type = 'response')
      pred <- ifelse(pred_prob > best_cut_off, 1, 0)
      mean(pred == label.test$Label)
      
      pred_prob2 <- predict(model, newx = as.matrix(data.test2), s = lambda_1se, type = 'response')
      pred2 <- ifelse(pred_prob2 > best_cut_off, 1, 0)
      mean(pred2 == 1)
    }
    else {
      model <- glm(label.train$Label ~., data = data.train, family = binomial)
      pred_prob <- predict(model, newdata = data.test, type = 'response')
      pred <- ifelse(pred_prob > 0.5, 1, 0)
      
      pred_prob2 <- predict(model, newdata = data.test2, type = 'response')
      pred2 <- ifelse(pred_prob2 > 0.5, 1, 0)
    }
    metrics <- compute_metrics(pred = pred, pred_prob = pred_prob, true_label = label.test$Label, classes = c(0, 1))  
    print(metrics)
    result_df1 <- data.frame("TestDataClassName" = ifelse(label.test$Label == 0, classes[1], classes[2]),
                            "TestDataExpectedClassName" = ifelse(label.test$Label == 0, classes[1], classes[2]),
                            "Pred_prob" = pred_prob[,1],
                            "Prediction" = pred[,1])
    
    result_df2 <- data.frame("TestDataClassName" = "PREREC",
                             "TestDataExpectedClassName" = "REC_TP", 
                             "Pred_prob" = pred_prob2[,1],
                             "Prediction" = pred2[,1])
    
    result_df <- rbind(result_df1, result_df2) %>%
      mutate(Pred_prob = as.double(Pred_prob)) %>%
      mutate(Prediction = ifelse(Prediction == 0, classes[1], classes[2]))
    colnames(result_df)[4] <- paste0("prediction_with_cutoff_", best_cut_off) 
    
    # result_file_name <- "Data/prediction_result/transcriptomics.csv"
    # write.csv(result_df, result_file_name) 
    
    return (result_df)
  })
  

}


