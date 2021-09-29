library(FEMPipeline)
source("dataset_pipeline_arguments.R")
# base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
# setwd(base_dir)
setwd("../..")

args = commandArgs(trailingOnly = TRUE)
if (length(args) > 1) {
  print(paste('Executing pipeline on dataset', args[2]))
  dparg = dataset_pipeline_arguments[[strtoi(args[2])]]
  do.call(execute_pipeline, dparg)
} else {
  print('Executing pipeline on all datasets')
  for (dparg in dataset_pipeline_arguments) {
    do.call(execute_pipeline, dparg)
  }
}



