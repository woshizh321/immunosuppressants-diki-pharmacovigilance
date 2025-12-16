CFG <- list(
  locked_model_json = "05_molecule_level_lasso/outputs/models/LockedModelSpec_molecule_level.json",
  
  ## Population grid
  sex_levels = c("Female", "Male"),
  age_groups = c("18-39", "40-64", "65+"),
  indications = c("Transplantation", "Autoimmune diseases"),
  
  ## Molecules
  molecules = c("Cyclosporine", "MPA", "MMF", "Azathioprine"),
  
  ## Baseline prevalence (used for absolute risk calibration)
  baseline_risk = 0.05,
  
  ## Output
  out_dir = "06_risk_portraits/outputs"
)
