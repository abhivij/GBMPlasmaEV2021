source("perform_de_comparisons_arglist.R")
source("perform_de_comparisons.R")

args = commandArgs(trailingOnly = TRUE)
if(length(args) > 1) {
  print(paste('executing perform de comparison with args', args[2]))
  arg_vals <- perform_de_comparisons_arglist[[strtoi(args[2])]]
  do.call(perform_de_comparisons, arg_vals)
} else {
  print('arg index not specified !')
}