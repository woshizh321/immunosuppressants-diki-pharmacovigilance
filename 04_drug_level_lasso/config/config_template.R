CFG <- list(
  faers_case_features = "PATH/TO/FAERS_case_features_drug_level.parquet",
  
  ## LASSO
  alpha = 1,
  nfolds = 10,
  lambda_rule = "lambda.1se",
  
  ## Model options
  include_region = FALSE,
  include_tto = FALSE,
  
  ## Output
  out_dir = "04_drug_level_lasso/outputs"
)
