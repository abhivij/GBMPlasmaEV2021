library(FEMPipeline)
source("scripts/R/dataset_pipeline_arguments_transcriptomic.R")

args = commandArgs(trailingOnly = TRUE)
if (length(args) > 1) {
  print(paste('Executing pipeline on dataset', args[2]))
  dparg = dataset_pipeline_arguments_transcriptomic[[strtoi(args[2])]]
  do.call(execute_pipeline, dparg)
} else {
  print('Executing pipeline on all datasets')
  for (dparg in dataset_pipeline_arguments_transcriptomic) {
    do.call(execute_pipeline, dparg)
  }
}