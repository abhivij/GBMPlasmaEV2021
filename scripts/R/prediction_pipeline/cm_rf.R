# library(randomForest)
# source("R/metrics/compute_metrics.R")

# data.train, output_labels, data.test, output_labels_test
# 
# x.train <- data.train
# y.train <- output_labels
# x.test <- data.test
# y.test <- output_labels_test
# random_seed <- 1000

rf_model <- function(x.train, y.train, x.test, y.test, classes, 
                     random_seed = 1000, ...){
  model_name <- "Random Forest"
  # #setting default value for metrics, to handle case where unable to train / execute classification model
  # metrics.test <- c(0, 0, 0, 0, 0, 0, 0, 0)
  # metrics.train <- c(0, 0, 0, 0, 0, 0, 0, 0)
  
  try({
    set.seed(random_seed)
    if(sum(x.train - colMeans(x.train)) != 0){
      #ensure that atleast one column is not a constant vector
      #all columns constant causes the below line to run forever
      model <- randomForest::randomForest(x = x.train, y = factor(y.train$Label, levels = classes))
      
      pred_prob.train <- predict(model, x.train, type="prob")
      pred_prob.train <- data.frame(pred_prob.train)[classes[2]]
      
      
      #identify best cutoff
      #cutoff such that max accuracy
      # if same max accuracy then min(cutoff - max(pred class 1))
      max_acc <- 0
      min_cutoff_maxpred_diff <- -1
      best_cutoff <- 0.1
      for(cutoff in seq(0.1, 0.8, by = 0.05)){
        print(cutoff)
        pred.train <- ifelse(pred_prob.train > cutoff, classes[2], classes[1])
        acc <- mean(pred.train[, 1] == y.train$Label)  
        if(max_acc < acc){
          best_cutoff <- cutoff
          max_acc <- acc
          cutoff_maxpred_diff <- cutoff - max(pred_prob.train[! pred_prob.train > cutoff])
          min_cutoff_maxpred_diff <- cutoff_maxpred_diff 
          print("*******")
          print(cutoff)
          print(acc)
        } else if(max_acc == acc){
          cutoff_maxpred_diff <- cutoff - max(pred_prob.train[! pred_prob.train > cutoff])
          if(min_cutoff_maxpred_diff == -1 || cutoff_maxpred_diff < min_cutoff_maxpred_diff){
            best_cutoff <- cutoff
            min_cutoff_maxpred_diff <- cutoff_maxpred_diff 
            print("-------")
            print(cutoff)
            print(cutoff_maxpred_diff)
          }
        }
      }
      print(best_cutoff)
      print(max_acc)
      print(min_cutoff_maxpred_diff)
      
      pred.train <- ifelse(pred_prob.train > best_cutoff, classes[2], classes[1])
      
      
      result_df <- data.frame("Sample" = y.train$Sample, "TrueLabel" = y.train$Label,
                              "PredictedLabel" = pred.train[, 1],
                              "Probability" = pred_prob.train[, 1],
                              "Type" = "train", 
                              "cutoff" = best_cutoff)

      pred_prob.test <- predict(model, x.test, type="prob")
      pred_prob.test <- data.frame(pred_prob.test)[classes[2]]
      pred.test <- ifelse(pred_prob.test > best_cutoff, classes[2], classes[1])
      
      result_df_test <- data.frame("Sample" = y.test$Sample, "TrueLabel" = y.test$Label,
                                   "PredictedLabel" = pred.test[, 1],
                                   "Probability" = pred_prob.test[, 1],
                                   "Type" = "test", "cutoff" = best_cutoff)
      
    } else{
      print("data to RF : all fields constant!")
      print(dim(x.train))
    }
    
  })
  
  return (rbind(result_df, result_df_test))
}
