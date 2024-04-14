source("scripts/R/RFE_execution_args.R")
source("scripts/R/de_ml_biomarker_discovery.R")

args = commandArgs(trailingOnly = TRUE)
if(length(args) > 1) {
  print(paste('RFE execution with arg', args[2]))
  arg_vals <- RFE_execution_args[[strtoi(args[2])]]
  do.call(RFE_from_ranked_list, arg_vals)
} else {
  print('arg index not specified !')
}
