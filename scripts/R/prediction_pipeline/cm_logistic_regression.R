log_reg_model <- function(data.train, label.train, data.test, label.test, 
                          classes, regularize = NA, ...){
  
  model_name <- "logistic regression"
  if (is.na(regularize)) {
    model_name <- paste("Simple", model_name)
  } else if(regularize == 'l1') {
    model_name <- paste("L1 Regularized", model_name)
  } else {
    model_name <- paste("L2 Regularized", model_name)
  }
  
  try({
    label.train$Label <- ifelse(label.train$Label == classes[1], 0, 1)
    #label.test contains UNK also. So can't change labels to 1 or 0  
    #label.test$Label <- ifelse(label.test$Label == classes[1], 0, 1)
    
    if (!is.na(regularize)) {
      #alpha = 1 => l1 regularization (lasso)
      #alpha = 0 => l2 regularization (ridge)
      #alpha = 0.5 => elastic net regularization (combination of lasso and ridge)
      if(regularize == "el"){
        alpha <- 0.5
      }
      else if(regularize == 'l1') {
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
      best_cut_off <- 0.5
      for(cut_off in c(0.5, seq(0.2, 0.8, 0.01))){
        print(cut_off)
        pred_prob.train <- predict(model, newx = as.matrix(data.train), s = lambda_1se, type = 'response')
        pred.train <- ifelse(pred_prob.train > cut_off, 1, 0)
        acc <- mean(pred.train == label.train$Label)
        if(acc > best_acc){
          best_acc <- acc
          best_cut_off <- cut_off
        }
      }
      
      pred_prob.train <- predict(model, newx = as.matrix(data.train), s = lambda_1se, type = 'response')
      pred.train <- ifelse(pred_prob.train > best_cut_off, 1, 0)
      mean(pred.train == label.train$Label)

      #label.test contains UNK also. So can't change labels to 1 or 0      
      #           and therefore can't compare pred and Label
      pred_prob <- predict(model, newx = as.matrix(data.test), s = lambda_1se, type = 'response')
      pred <- ifelse(pred_prob > best_cut_off, 1, 0)
      # mean(pred == label.test$Label)
      
      
      result_df.train <- data.frame("TrueLabel" = ifelse(label.train$Label == 0, classes[1], classes[2]),
                                    "Pred_prob" = pred_prob.train[,1],
                                    "PredictedLabel" = pred.train[,1],
                                    "Type" = "train",
                                    "cutoff" = best_cut_off)
      
      result_df.test <- data.frame("TrueLabel" = label.test$Label,
                                   "Pred_prob" = pred_prob[,1],
                                   "PredictedLabel" = pred[,1],
                                   "Type" = "test",
                                   "cutoff" = best_cut_off)
    }
    else {
      model <- glm(label.train$Label ~., data = data.train, family = binomial)
      
      best_acc <- -1
      best_cut_off <- 0.5
      for(cut_off in c(0.5, seq(0.2, 0.8, 0.01))){
        print(cut_off)
        pred_prob.train <- predict(model, newdata = data.train, type = 'response')
        pred.train <- ifelse(pred_prob.train > cut_off, 1, 0)
        acc <- mean(pred.train == label.train$Label)
        if(acc > best_acc){
          best_acc <- acc
          best_cut_off <- cut_off
        }
      }

      pred_prob.train <- predict(model, newdata = data.train, type = 'response')
      pred.train <- ifelse(pred_prob.train > best_cut_off, 1, 0)
            
      pred_prob <- predict(model, newdata = data.test, type = 'response')
      pred <- ifelse(pred_prob > best_cut_off, 1, 0)
      
      result_df.train <- data.frame("TrueLabel" = ifelse(label.train$Label == 0, classes[1], classes[2]),
                                    "Pred_prob" = pred_prob.train,
                                    "PredictedLabel" = pred.train,
                                    "Type" = "train",
                                    "cutoff" = best_cut_off)
      
      result_df.test <- data.frame("TrueLabel" = label.test$Label,
                                   "Pred_prob" = pred_prob,
                                   "PredictedLabel" = pred,
                                   "Type" = "test",
                                   "cutoff" = best_cut_off)
    }
    

    
    result_df <- rbind(result_df.train %>%
                         rownames_to_column("sample"), 
                       result_df.test %>%
                         rownames_to_column("sample"))
    
    result_df <- result_df %>%
      mutate(Pred_prob = as.double(Pred_prob)) %>%
      mutate(PredictedLabel = ifelse(PredictedLabel == 0, classes[1], classes[2])) %>%
      mutate(model = model_name)
    
    return (result_df)
  })
  
  
}


