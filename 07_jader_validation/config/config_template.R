CFG <- list(
  jader_master_parquet = "D:/JADER/MASTER/JADER_MASTER_PT_English_SUPERMASTER_v3_patched_clean.parquet",
  
  ## molecule scope
  molecules = c("Cyclosporine", "Tacrolimus", "MMF", "MPA", "Azathioprine"),
  
  ## outcome
  diki_definition = "16 MedDRA PT codes (same as FAERS)",
  
  ## role filter (if applicable)
  role_included = c("PS","SS"),  # adjust if JADER uses different coding
  
  ## external validation
  locked_model_json = "05_molecule_level_lasso/outputs/models/LockedModelSpec_molecule_level.json",
  
  ## exclusions
  include_tto = FALSE,
  include_region = FALSE,
  
  ## outputs
  out_dir = "07_jader_validation/outputs"
)
