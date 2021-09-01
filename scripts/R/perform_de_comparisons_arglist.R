perform_de_comparisons_arglist <- list(
  #1
  list(),
  
  #2
  list(file_name = "data_process_output_annotatedQ1-6_NA_equalizeMedians.rds", 
       comparison_list = list(c("PREOPE", "MET"),
                              c("PREOPE", "HC"),
                              c("MET", "HC")), 
       comparison_num = "2"),  

  #3
  list(file_name = "data_process_output_annotatedQ1-6_NA_equalizeMedians.rds", 
       comparison_list = list(c("PREOPE", "POSTOPE-T"),
                              c("PREOPE", "POSTOPE-P"),
                              c("POSTOPE-T", "POSTOPE-P")), 
       comparison_num = "3"), 
  
  #4
  list(file_name = "data_process_output_annotatedQ1-6_NA_equalizeMedians.rds", 
       comparison_list = list(c("POSTOPE-T", "REC-T")), 
       comparison_num = "4"),   
  
  #5
  list(file_name = "data_process_output_annotatedQ1-6_NA_equalizeMedians.rds", 
       comparison_list = list(c("POSTOPE-P", "REC-P")), 
       comparison_num = "5"), 

  #6
  list(file_name = "data_process_output_annotatedQ1-6_NA_equalizeMedians.rds", 
       comparison_list = list(c("POSTOPE-T", "PREREC")), 
       comparison_num = "6"),     
  
  #7
  list(file_name = "data_process_output_annotatedQ7_NA_equalizeMedians.rds", 
       comparison_list = list(c("PREOPE", "REC-TP")), 
       comparison_num = "7")  
)