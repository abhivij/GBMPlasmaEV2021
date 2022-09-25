process_protein_arglist <- list(
  #1
  list(data_dir = "annotatedQ1-6",
       condition_type = NA,
       norm = "equalizeMedians"),
  #2
  list(data_dir = "annotatedQ1-6",
       condition_type = NA,
       norm = FALSE),
  #3
  list(data_dir = "annotatedQ1-6",
       condition_type = "column",
       norm = "equalizeMedians"),
  #4
  list(data_dir = "annotatedQ1-6",
       condition_type = "column",
       norm = FALSE),
  #5
  list(data_dir = "annotatedQ1-6",
       condition_type = "disease",
       norm = "equalizeMedians"),
  #6
  list(data_dir = "annotatedQ1-6",
       condition_type = "disease",
       norm = FALSE),  
  
  #7
  list(data_dir = "annotatedQ7",
       condition_type = NA,
       norm = "equalizeMedians"),
  #8
  list(data_dir = "annotatedQ7",
       condition_type = NA,
       norm = FALSE),
  #9
  list(data_dir = "annotatedQ7",
       condition_type = "column",
       norm = "equalizeMedians"),
  #10
  list(data_dir = "annotatedQ7",
       condition_type = "column",
       norm = FALSE),
  #11
  list(data_dir = "annotatedQ7",
       condition_type = "disease",
       norm = "equalizeMedians"),
  #12
  list(data_dir = "annotatedQ7",
       condition_type = "disease",
       norm = FALSE),
  
  #13
  list(data_dir = "unannotated",
       condition_type = "column",
       norm = "equalizeMedians"),
  #14
  list(data_dir = "unannotated",
       condition_type = "column",
       norm = FALSE),
  #15
  list(data_dir = "unannotated",
       condition_type = "disease",
       norm = "equalizeMedians"),
  #16
  list(data_dir = "unannotated",
       condition_type = "disease",
       norm = FALSE),
  
  #17
  list(data_dir = "annotatedQ1-6",
       condition_type = NA,
       norm = "quantile"),  
  #18
  list(data_dir = "annotatedQ1-6",
       condition_type = "column",
       norm = "quantile"),
  #19
  list(data_dir = "annotatedQ1-6",
       condition_type = "disease",
       norm = "quantile"), 
  #20
  list(data_dir = "annotatedQ7",
       condition_type = NA,
       norm = "quantile"),  
  
  
  #21
  list(data_dir = "annotatedQ1-6",
       condition_type = NA,
       norm = "equalizeMedians",
       remove_50_missing = TRUE),   
  #22
  list(data_dir = "annotatedQ1-6",
       condition_type = NA,
       norm = FALSE,
       remove_50_missing = TRUE),
  #23
  list(data_dir = "annotatedQ1-6",
       condition_type = NA,
       norm = "quantile",
       remove_50_missing = TRUE),
  
  
  #24
  list(data_dir = "_newcohort_SBGNPlasma",
       condition_type = NA,
       norm = FALSE,
       remove_50_missing = FALSE)  
)
