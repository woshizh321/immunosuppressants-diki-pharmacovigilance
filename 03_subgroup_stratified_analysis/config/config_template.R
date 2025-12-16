CFG <- list(
  faers_master = "PATH/TO/FAERS_MASTER.parquet",
  out_dir = "03_subgroup_stratified_analysis/outputs",
  
  ## General filters
  roles_keep = c("PS","SS","C"),
  min_cell_count = 3,
  
  ## Subgroup definitions
  sex_levels = c("Male","Female"),
  age_levels = c("<18","18-39","40-64","65+"),
  seriousness_levels = c("Serious","Non-serious"),
  
  ## Region handling
  region_level = "COUNTRY",  # or REGION_GROUP
  
  ## TTO
  tto_bins = c("0-7","8-30","31-90","91-365",">365")
)
