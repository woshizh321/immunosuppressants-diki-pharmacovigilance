CFG <- list(
  case_features = "PATH/TO/FAERS_case_features_molecule_level.parquet",
  
  ## LASSO
  alpha = 1,
  nfolds = 10,
  lambda_rule = "lambda.1se",
  
  ## Feature space
  molecule_scope = c("CNI", "APA"),
  min_exposed_cases = 50,
  
  ## Covariates
  include_sex = TRUE,
  include_age = TRUE,
  include_indication = TRUE,
  
  ## Exclusions
  include_region = FALSE,
  include_tto = FALSE,
  
  ## Output
  out_dir = "05_molecule_level_lasso/outputs"
)
