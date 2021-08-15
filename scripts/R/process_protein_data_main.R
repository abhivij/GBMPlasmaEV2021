source("process_protein_data.R")
source("process_protein_arglist.R")

args = commandArgs(trailingOnly = TRUE)
if(length(args) > 1) {
  print(paste('executing process protein with args', args[2]))
  arg_vals <- process_protein_arglist[[strtoi(args[2])]]
  do.call(process_protein_data, arg_vals)
} else {
  print('arg index not specified !')
}


