CFG <- list(
  faers_master = "PATH/TO/FAERS_MASTER.parquet",
  out_dir = "02_hierarchical_signal_detection/outputs",
  ## optional: restrict roles if you do so
  roles_keep = c("PS","SS","C"),
  ## signal metrics
  min_cell_count = 3,
  ## hierarchy levels to compute
  levels = c("SOC","HLGT","HLT","PT")
)
