#!/usr/bin/env Rscript
## ============================================================================
## 07_fig7_part2_build_case_features_jader.R
##
## Module A: Generate LockedModelSpec JSON from molecule-level Refit OR results
## Module B: Generate case_features (molecule-level CNI+APA) from FAERS / JADER MASTER
## Module C: Perform external validation (Performance + Calibration) on FAERS / JADER 
##           using the LockedModel
##
## Notes:
##   1) This script uses the same DIKI definitions and indication keywords as 
##      the LASSO analysis.
##   2) Column names for JADER MASTER are assumed to be consistent with standardized 
##      FAERS MASTER (e.g., sex_std, age_years, indi_pt_norm). If inconsistent, 
##      manual adjustment may be required.
## ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(readr)
  library(arrow)
  library(duckdb)
  library(DBI)
  library(ggplot2)
  library(pROC)
  library(PRROC)
  library(scales)
  library(Hmisc)
  library(tibble)
})

options(arrow.use_threads = TRUE)
data.table::setDTthreads(6)

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

## ========================== Global Parameters ==========================

## --- A. Refit OR File & LockedModelSpec Output ---
REFIT_OR_FILE <- "D:/FAERS/Inhibitors/Derived_molecule_withIndi_noyear/lasso_refit_logistic_results_molecule_CNI_APA_withIndication_noRegion_noTTO_noSerious_noyear.csv"
LOCKED_JSON   <- "D:/FAERS/Inhibitors/Derived_molecule_withIndi_noyear/LockedModelSpec_molecule_CNI_APA_noyear.json"

## --- B. MASTER Files & case_features Output Paths ---
FAERS_MASTER_PARQUET <- "D:/FAERS/MASTER/FAERS_MASTER_FILE_2004-2024_with_serious.parquet"

## Updated to v4: Includes DRUG_INN_EN / DRUG_MOLECULE / DRUG_CLASS_STD / DRUG_CODE_STD
JADER_MASTER_PARQUET <- "D:/JADER/MASTER/JADER_MASTER_PT_English_SUPERMASTER_v4_withDrugCode.parquet"

OUT_CASE_FAERS <- "D:/FAERS/ML_outputs_duckdb/FAERS_case_features_molecule_CNI_APA_v1.parquet"
OUT_CASE_JADER <- "D:/FAERS/ML_outputs_duckdb/JADER_case_features_molecule_CNI_APA_v2_fromMASTERv4.parquet"


## --- C. External Validation Output ---
BASE_EXTVAL_DIR <- "D:/JADER/ExternalValidation_molecule_CNI_APA"

dir.create(dirname(LOCKED_JSON), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(OUT_CASE_FAERS), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(OUT_CASE_JADER), showWarnings = FALSE, recursive = TRUE)
dir.create(BASE_EXTVAL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(BASE_EXTVAL_DIR, "predictions"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(BASE_EXTVAL_DIR, "Evaluation"),  showWarnings = FALSE, recursive = TRUE)

## ===================== Utility Functions: General =====================

must_exist <- function(path, label = deparse(substitute(path))) {
  if (!file.exists(path)) stop(sprintf("[%s] File does not exist: %s", label, path))
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

preview_parquet <- function(path, n_head = 5L) {
  path <- must_exist(path)
  cat("\n==== Previewing File:", path, "====\n")
  DT <- as.data.table(arrow::read_parquet(path))
  cat("Total Rows:", nrow(DT), "  Total Columns:", ncol(DT), "\n")
  cat("First 30 Column Names:\n")
  print(head(names(DT), 30))
  
  if ("DIKI" %in% names(DT)) {
    cat("\nDIKI Distribution:\n")
    print(table(DT$DIKI, useNA = "ifany"))
  }
  
  cat("\nFirst", n_head, "rows preview:\n")
  print(head(DT, n_head))
  invisible(DT)
}

## ===================== Module A: Lock Model =====================

lock_model_from_refit <- function(refit_csv, locked_json) {
  refit_csv  <- must_exist(refit_csv, "REFIT_OR_FILE")
  dt_or <- fread(refit_csv)
  
  if (!all(c("term","OR") %in% names(dt_or))) {
    stop("Refit OR file must contain at least 'term' and 'OR' columns.")
  }
  
  ## Intercept
  OR_intercept <- dt_or[term == "(Intercept)", OR][1]
  if (is.na(OR_intercept)) stop("Could not find '(Intercept)' row in Refit OR file.")
  beta0 <- log(OR_intercept)
  
  ## Covariates: Remove intercept and id_case
  dt_beta <- dt_or[!(term %in% c("(Intercept)", "id_case"))]
  
  if (nrow(dt_beta) == 0L) {
    stop("No usable covariates found in Refit OR file (excluding intercept and id_case).")
  }
  
  dt_beta[, beta := log(OR)]
  coef_vec <- dt_beta$beta
  names(coef_vec) <- dt_beta$term
  
  cat("\n>>> Covariates to be locked (term, OR, beta):\n")
  print(dt_beta[, .(term, OR, beta)])
  
  locked_spec <- list(
    model_id   = "FAERS_Refit_v1_Molecule_CNI_APA_noYear_noRegion_noTTO_noSerious",
    family     = "binomial",
    link       = "logit",
    intercept  = unname(beta0),
    coefficients = coef_vec,
    training_meta = list(
      source_dataset = "FAERS",
      years          = "2004Q1–2024Q4",
      outcome        = "DIKI (16 PT codes)",
      exposure       = "CNI + APA (molecule-level, PS/SS only)",
      covariates     = names(coef_vec),
      refit_type     = "logistic_refit_at_lambda.1se_nonzero_features",
      note           = "LASSO used for variable selection; predictions use refit logistic coefficients."
    )
  )
  
  write_json(
    locked_spec,
    path       = locked_json,
    pretty     = TRUE,
    auto_unbox = TRUE,
    digits     = 10
  )
  
  cat("\n✅ LockedModelSpec generated:\n  ", locked_json, "\n")
  cat("Intercept (beta0) =", beta0, "\n")
  invisible(locked_spec)
}

## ===================== Module B: Generate case_features =====================

## 16 DIKI PT codes (consistent with LASSO script)
diki_pt_codes <- c(
  10001580, # Albuminuria
  10002847, # Anuria
  10005483, # Blood creatinine increased
  10018358, # Glomerular filtration rate decreased
  10018867, # Haematuria
  10027525, # Microalbuminuria
  10029155, # Nephropathy toxic
  10030302, # Oliguria
  10037032, # Proteinuria
  10038428, # Renal disorder
  10038435, # Renal failure
  10061480, # Renal function test abnormal
  10062237, # Renal impairment
  10062747, # Hypercreatininemia
  10064848, # Chronic kidney disease
  10069339  # Acute kidney injury
)

## Molecule dictionary: CNI + APA only
molecule_dict <- list(
  "Cyclosporine" = c("CYCLOSPORINE","NEORAL","SANDIMMUNE"),
  "MMF"          = c("MYCOPHENOLATE MOFETIL","CELLCEPT"),
  "MPA"          = c("MYCOPHENOLIC ACID","MYFORTIC"),
  "Azathioprine" = c("AZATHIOPRINE","IMURAN")
)

mol_flag_names <- c(
  "flag_mol_cyclosporine",
  "flag_mol_mmf",
  "flag_mol_mpa",
  "flag_mol_azathioprine"
)

## Indication keywords
kw_transplant <- c(
  "TRANSPLANT","TRANSPLANTATION","GRAFT","RENAL TRANSPLANT","KIDNEY TRANSPLANT",
  "HEART TRANSPLANT","LIVER TRANSPLANT","PANCREAS TRANSPLANT","BONE MARROW TRANSPLANT",
  "HSCT","SCT"
)
kw_autoimmune <- c(
  "SYSTEMIC LUPUS","LUPUS","SLE","RHEUMATOID ARTHRITIS","RA","PSORIASIS","PSORIATIC",
  "ANKYLOSING SPONDYLITIS","CROHN","ULCERATIVE COLITIS","IBD","INFLAMMATORY BOWEL",
  "VASCULITIS","ANCA","MULTIPLE SCLEROSIS","MYASTHENIA","SJOGREN","SCLERODERMA",
  "DERMATOMYOSITIS","POLYMYOSITIS","AUTOIMMUNE","IGA NEPHROPATHY","LUPUS NEPHRITIS"
)
pat_tx <- paste0("(", paste(kw_transplant, collapse = "|"), ")")
pat_ai <- paste0("(", paste(kw_autoimmune, collapse = "|"), ")")

detect_column <- function(cols, candidates, required = TRUE, label = "column") {
  hit <- intersect(candidates, cols)
  if (length(hit) >= 1L) return(hit[1L])
  if (required) stop(sprintf("Could not find any %s candidate columns in MASTER: %s", label, paste(candidates, collapse = ", ")))
  NA_character_
}

build_case_features_from_master <- function(
    master_parquet,
    out_parquet,
    source_tag = c("FAERS","JADER")
) {
  source_tag <- match.arg(source_tag)
  master_parquet <- must_exist(master_parquet, paste0(source_tag, "_MASTER_PARQUET"))
  
  cat("\n================ Building case_features (", source_tag, ") ================\n")
  con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:", read_only = TRUE)
  on.exit(try(dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)
  
  cols0 <- dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s') LIMIT 0", master_parquet))
  all_cols <- names(cols0)
  cat("MASTER Column Count:", length(all_cols), "\n")
  cat("First 40 Column Names:\n")
  print(head(all_cols, 40))
  
  ## Case ID
  CASE_ID_COL <- detect_column(
    all_cols,
    c("id_case","caseid","case_id","PRIMARYID","primaryid","ID"),
    required = TRUE,
    label = "Case ID"
  )
  cat("Case ID Column:", CASE_ID_COL, "\n")
  
  ## pt_code
  PT_COL <- detect_column(all_cols, c("pt_code","PT_CODE","ptcode"), required = TRUE, label = "pt_code")
  cat("PT Column:", PT_COL, "\n")
  
  ## Drug Name Selection
  ## For FAERS: Continue using DRUGNAME / DRUGNAME_EN
  ## For JADER MASTER v4: Prioritize DRUG_MOLECULE (standardized name)
  if (source_tag == "JADER" && "DRUG_MOLECULE" %in% all_cols) {
    DRUG_COL_RAW <- "DRUG_MOLECULE"
    cat("Drug Name Column: Using JADER standardized molecule column DRUG_MOLECULE\n")
  } else {
    DRUG_COL_RAW <- detect_column(
      all_cols,
      c("drugname","DRUGNAME","drugname_en","DRUGNAME_EN"),
      required = TRUE,
      label = "drugname"
    )
    cat("Drug Name Column:", DRUG_COL_RAW, "\n")
  }
  
  ## Sex
  SEX_COL_RAW <- detect_column(
    all_cols,
    c("sex_raw","sex_std","sex","SEX_STD","SEX"),
    required = TRUE,
    label = "sex"
  )
  cat("Sex Column:", SEX_COL_RAW, "\n")
  
  ## Continuous Age
  AGE_YEARS_COL <- detect_column(
    all_cols,
    c("age_years_raw","age_years","AGE_YEARS","age","AGE"),
    required = FALSE,
    label = "age_years"
  )
  if (!is.na(AGE_YEARS_COL)) cat("Continuous Age Column:", AGE_YEARS_COL, "\n") else cat("No continuous age column detected. Granular stratification unavailable.\n")
  
  ## Age Group (Backup)
  AGE_GROUP_COL <- detect_column(
    all_cols,
    c("age_group_raw","age_group","AGE_GROUP"),
    required = FALSE,
    label = "age_group"
  )
  if (!is.na(AGE_GROUP_COL)) cat("Age Group Column:", AGE_GROUP_COL, "\n")
  
  ## Year (For QC, not necessarily in model)
  YEAR_COL <- detect_column(
    all_cols,
    c("year_raw","year","YEAR"),
    required = FALSE,
    label = "year"
  )
  if (!is.na(YEAR_COL)) cat("Year Column:", YEAR_COL, "\n")
  
  ## Serious Column (Filter for non-serious if exists)
  SERIOUS_COL_CAND <- grep("serious", all_cols, ignore.case = TRUE, value = TRUE)
  SERIOUS_COL <- NA_character_
  if ("serious_flag" %in% all_cols) {
    SERIOUS_COL <- "serious_flag"
  } else if (length(SERIOUS_COL_CAND) == 1L) {
    SERIOUS_COL <- SERIOUS_COL_CAND[1L]
  }
  cat("Seriousness Column:", SERIOUS_COL, "\n")
  
  ## Indication Text
  INDI_COL <- detect_column(
    all_cols,
    c("indi_pt_norm_raw","indi_pt_norm","INDI_PT_NORM","INDI_PT_NORM_RAW",
      "INDI_PT_NAME_EN","indi_pt_name_en","INDICATION_GROUP"),
    required = TRUE,
    label = "indication"
  )
  cat("Indication Text Column:", INDI_COL, "\n")
  
  ## Role Column (FAERS uses PS/SS)
  ROLE_COL <- detect_column(
    all_cols,
    c("role_cod","ROLE_COD","ROLE_STD"),
    required = FALSE,
    label = "role"
  )
  cat("Role Column:", ROLE_COL, "\n")
  
  ## ========== SQL Extraction (Target DIKI PTs or Molecule Drugs Only) ==========
  ## Optimization: Fetch only relevant records to avoid memory overflow
  
  drug_all_upper <- unique(unlist(molecule_dict))
  drug_all_upper <- sprintf("'%s'", drug_all_upper)
  drug_list_sql  <- paste(drug_all_upper, collapse = ", ")
  
  diki_list_sql  <- paste(diki_pt_codes, collapse = ", ")
  
  select_fields <- c(
    sprintf("%s AS id_case", CASE_ID_COL),
    sprintf("%s AS pt_code", PT_COL),
    sprintf("%s AS sex_raw", SEX_COL_RAW),
    sprintf("%s AS indi_raw", INDI_COL),
    sprintf("UPPER(%s) AS drug_upper", DRUG_COL_RAW)
  )
  if (!is.na(AGE_YEARS_COL)) {
    select_fields <- c(select_fields, sprintf("%s AS age_years_raw", AGE_YEARS_COL))
  }
  if (!is.na(AGE_GROUP_COL)) {
    select_fields <- c(select_fields, sprintf("%s AS age_group_raw", AGE_GROUP_COL))
  }
  if (!is.na(YEAR_COL)) {
    select_fields <- c(select_fields, sprintf("%s AS year_raw", YEAR_COL))
  }
  if (!is.na(SERIOUS_COL)) {
    select_fields <- c(select_fields, sprintf("%s AS serious_raw", SERIOUS_COL))
  }
  if (!is.na(ROLE_COL)) {
    select_fields <- c(select_fields, sprintf("%s AS role_raw", ROLE_COL))
  }
  
  where_parts <- c(
    sprintf("(%s IN (%s) OR UPPER(%s) IN (%s))", PT_COL, diki_list_sql, DRUG_COL_RAW, drug_list_sql)
  )
  if (!is.na(ROLE_COL)) {
    if (ROLE_COL %in% c("role_cod","ROLE_COD")) {
      where_parts <- c(sprintf("%s IN ('PS','SS')", ROLE_COL), where_parts)
    } else if (ROLE_COL == "ROLE_STD") {
      where_parts <- c(sprintf("%s IN ('PS','C')", ROLE_COL), where_parts)
    }
  }
  
  sql_main <- sprintf("
    SELECT
      %s
    FROM read_parquet('%s')
    WHERE %s
  ",
                      paste(select_fields, collapse = ",\n      "),
                      master_parquet,
                      paste(where_parts, collapse = " AND ")
  )
  
  cat("\nExecuting SQL (Preview first 400 chars):\n")
  cat(substr(sql_main, 1, 400), "...\n")
  
  dt_raw <- as.data.table(dbGetQuery(con, sql_main))
  cat("\nExtraction complete. dt_raw dimensions:", nrow(dt_raw), "rows x", ncol(dt_raw), "cols\n")
  if (nrow(dt_raw) == 0L) stop("No relevant records found. Check DIKI PT or Molecule dictionary settings.")
  
  ## Row-level DIKI flag
  dt_raw[, flag_diki := as.integer(pt_code %in% diki_pt_codes)]
  cat("Row-level flag_diki distribution:\n")
  print(table(dt_raw$flag_diki, useNA = "ifany"))
  
  ## ========== Case-level Aggregation ==========
  cat("\n>>> Building case-level case_features...\n")
  
  dt_case <- dt_raw[, {
    out <- list()
    
    out[["caseid"]] <- as.character(id_case[1])
    
    out[["DIKI"]] <- as.integer(any(flag_diki == 1L, na.rm = TRUE))
    
    out[["sex_raw"]]  <- sex_raw[which.max(!is.na(sex_raw))]
    out[["indi_raw"]] <- indi_raw[which.max(!is.na(indi_raw))]
    
    if ("age_years_raw" %in% names(.SD)) {
      ay <- age_years_raw[which.max(!is.na(age_years_raw))]
      out[["age_years_raw"]] <- suppressWarnings(as.numeric(ay))
    }
    if ("age_group_raw" %in% names(.SD)) {
      out[["age_group_raw"]] <- age_group_raw[which.max(!is.na(age_group_raw))]
    }
    if ("year_raw" %in% names(.SD)) {
      out[["year_raw"]] <- year_raw[which.max(!is.na(year_raw))]
    }
    if ("serious_raw" %in% names(.SD)) {
      out[["serious_flag"]] <- max(
        as.integer(serious_raw %in% c(1, "1", "Y", "YES", TRUE)),
        na.rm = TRUE
      )
    }
    
    out[["flag_mol_cyclosporine"]] <- as.integer(any(drug_upper %in% molecule_dict$Cyclosporine, na.rm = TRUE))
    out[["flag_mol_mmf"]]          <- as.integer(any(drug_upper %in% molecule_dict$MMF,          na.rm = TRUE))
    out[["flag_mol_mpa"]]          <- as.integer(any(drug_upper %in% molecule_dict$MPA,          na.rm = TRUE))
    out[["flag_mol_azathioprine"]] <- as.integer(any(drug_upper %in% molecule_dict$Azathioprine, na.rm = TRUE))
    
    out
  }, by = id_case]
  
  cat("Case-level dt_case dimensions:", nrow(dt_case), "rows x", ncol(dt_case), "cols\n")
  cat("Case-level DIKI distribution:\n")
  print(table(dt_case$DIKI, useNA = "ifany"))
  
  ## Filter for exposure to at least 1 target molecule
  expo_any <- dt_case[, rowSums(.SD, na.rm = TRUE), .SDcols = mol_flag_names]
  dt_case <- dt_case[expo_any > 0]
  cat("Cases exposed to at least one CNI+APA molecule:", nrow(dt_case), "\n")
  
  ## Non-serious filtering (if serious_flag exists)
  if ("serious_flag" %in% names(dt_case)) {
    cat("serious_flag distribution:\n")
    print(table(dt_case$serious_flag, useNA = "ifany"))
    dt_case <- dt_case[is.na(serious_flag) | serious_flag == 0L]
    cat("Case count after non-serious filtering:", nrow(dt_case), "\n")
  } else {
    cat("⚠ serious_flag not detected. case_features includes all cases (no non-serious filter).\n")
  }
  
  ## Sex dummy variables: sexMale / sexUnknown
  dt_case[, sex_std := toupper(trimws(as.character(sex_raw)))]
  dt_case[sex_std %in% c("M","MALE"),   sexMale := 1L]
  dt_case[sex_std %in% c("F","FEMALE"), sexMale := 0L]
  dt_case[is.na(sexMale),               sexMale := 0L]  ## Baseline: Female/Other
  
  dt_case[, sexUnknown := as.integer(!sex_std %in% c("M","MALE","F","FEMALE") & !is.na(sex_std))]
  
  cat("\nsex_std distribution:\n")
  print(table(dt_case$sex_std, useNA = "ifany"))
  cat("sexMale distribution:\n");    print(table(dt_case$sexMale, useNA = "ifany"))
  cat("sexUnknown distribution:\n"); print(table(dt_case$sexUnknown, useNA = "ifany"))
  
  ## Age Categorization & age_cat40.64
  if ("age_years_raw" %in% names(dt_case)) {
    dt_case[, age_years_raw := suppressWarnings(as.numeric(age_years_raw))]
    dt_case[is.na(age_years_raw) | age_years_raw <= 0, age_cat := NA_character_]
    dt_case[!is.na(age_years_raw) & age_years_raw < 18,                   age_cat := "<18"]
    dt_case[!is.na(age_years_raw) & age_years_raw >= 18 & age_years_raw <= 39, age_cat := "18-39"]
    dt_case[!is.na(age_years_raw) & age_years_raw >= 40 & age_years_raw <= 64, age_cat := "40-64"]
    dt_case[!is.na(age_years_raw) & age_years_raw >= 65,                age_cat := "65+"]
    
    dt_case[, age_cat := factor(age_cat, levels = c("<18","18-39","40-64","65+"))]
  } else if ("age_group_raw" %in% names(dt_case)) {
    dt_case[, age_cat := factor(age_group_raw)]
  } else {
    dt_case[, age_cat := NA_character_]
  }
  
  cat("\nage_cat distribution:\n")
  print(table(dt_case$age_cat, useNA = "ifany"))
  
  dt_case[, age_cat40.64 := as.integer(age_cat == "40-64")]
  cat("age_cat40.64 distribution:\n")
  print(table(dt_case$age_cat40.64, useNA = "ifany"))
  
  ## indication_grp & dummy variables
  dt_case[, indi_up := toupper(trimws(as.character(indi_raw)))]
  dt_case[, indication_grp := "Other"]
  dt_case[grepl(pat_tx, indi_up), indication_grp := "Transplantation"]
  dt_case[indication_grp == "Other" & grepl(pat_ai, indi_up), indication_grp := "Autoimmune diseases"]
  dt_case[, indication_grp := factor(indication_grp,
                                     levels = c("Other","Transplantation","Autoimmune diseases"))]
  cat("\nindication_grp distribution:\n")
  print(table(dt_case$indication_grp, useNA = "ifany"))
  
  dt_case[, indication_grpTransplantation := as.integer(indication_grp == "Transplantation")]
  cat("indication_grpTransplantation distribution:\n")
  print(table(dt_case$indication_grpTransplantation, useNA = "ifany"))
  
  ## Final case_features: Retain columns required by Locked model + basic IDs
  dt_out <- dt_case[, .(
    caseid,
    DIKI = as.integer(DIKI),
    sexMale       = as.integer(sexMale),
    sexUnknown    = as.integer(sexUnknown),
    age_cat40.64  = as.integer(age_cat40.64),
    indication_grpTransplantation = as.integer(indication_grpTransplantation),
    flag_mol_cyclosporine  = as.integer(flag_mol_cyclosporine),
    flag_mol_mmf           = as.integer(flag_mol_mmf),
    flag_mol_mpa           = as.integer(flag_mol_mpa),
    flag_mol_azathioprine  = as.integer(flag_mol_azathioprine),
    source = source_tag
  )]
  
  cat("\n>>> Final case_features (", source_tag, ") dimensions:", nrow(dt_out), "rows x", ncol(dt_out), "cols\n")
  cat("DIKI distribution:\n")
  print(table(dt_out$DIKI, useNA = "ifany"))
  
  arrow::write_parquet(dt_out, out_parquet)
  cat("✅ case_features written to:\n  ", out_parquet, "\n")
  
  invisible(dt_out)
}

## ===================== Module C: External Validation =====================

read_locked_spec <- function(path_json){
  path_json <- must_exist(path_json, "locked_model_spec")
  js <- jsonlite::fromJSON(readr::read_file(path_json), simplifyVector = TRUE)
  if (is.null(js$coefficients)) stop("LockedModelSpec missing 'coefficients' field.")
  list(
    intercept = unname(js$intercept %||% 0),
    coefs     = unlist(js$coefficients)
  )
}

evaluate_and_plot <- function(idy, pred, out_png, tag = "") {
  y <- as.integer(idy$DIKI)
  p <- pmin(pmax(pred, 1e-6), 1 - 1e-6)
  
  ## Handle all NA predictions
  if (all(is.na(p))) {
    warning(sprintf("[%s] All predictions are NA. Evaluation impossible.", tag))
    calib_sum <- data.table(bin = factor(), pred = numeric(), obs = numeric(), n = integer())
    brier <- NA_real_; auroc <- NA_real_; auprc <- NA_real_
    
    g <- ggplot(calib_sum, aes(x = pred, y = obs)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      labs(
        title    = "Calibration (deciles)",
        subtitle = "All predictions are NA",
        x        = "Predicted risk",
        y        = "Observed risk"
      )
    ggplot2::ggsave(out_png, g, width = 5.5, height = 5.5, dpi = 300)
    
    return(list(AUROC = auroc, AUPRC = auprc, Brier = brier, CalibrationTable = calib_sum))
  }
  
  y_valid <- y[is.finite(y)]
  unique_y <- unique(y_valid)
  
  ## Brier Score
  brier <- mean((p - y)^2, na.rm = TRUE)
  
  ## Calibration Binning
  calib <- data.table(y = y, p = p)
  calib <- calib[is.finite(y) & is.finite(p)]
  
  calib[, bin := tryCatch(
    Hmisc::cut2(p, g = 10),
    error = function(e) {
      message(sprintf("[%s] cut2 failed: %s; falling back to quantile binning.", tag, e$message))
      brks <- quantile(p, probs = seq(0, 1, length.out = 11), na.rm = TRUE)
      cut(p, breaks = brks, include.lowest = TRUE, dig.lab = 6)
    }
  )]
  
  calib_sum <- calib[, .(
    pred = mean(p),
    obs  = mean(y),
    n    = .N
  ), by = bin]
  
  ## ROC/AUC + PR Curves
  auroc <- NA_real_; auprc <- NA_real_
  if (length(unique_y) >= 2L) {
    roc_obj <- pROC::roc(y, p, quiet = TRUE, direction = "<")
    auroc   <- as.numeric(pROC::auc(roc_obj))
    pr      <- PRROC::pr.curve(scores.class0 = p[y == 1], scores.class1 = p[y == 0], curve = TRUE)
    auprc   <- pr$auc.integral
    subtitle_txt <- sprintf("AUROC=%.3f  AUPRC=%.3f  Brier=%.3f", auroc, auprc, brier)
  } else {
    msg <- sprintf("[%s] DIKI has only one value (all 0 or 1). AUROC/AUPRC unavailable. Reporting Brier/Calibration only.", tag)
    message(msg)
    subtitle_txt <- sprintf("AUROC=NA  AUPRC=NA  Brier=%.3f", brier)
  }
  
  g <- ggplot(calib_sum, aes(x = pred, y = obs)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point(size = 2) +
    geom_line() +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    coord_fixed() +
    labs(
      title    = "Calibration (deciles)",
      subtitle = subtitle_txt,
      x        = "Predicted risk",
      y        = "Observed risk"
    )
  
  ggplot2::ggsave(out_png, g, width = 5.5, height = 5.5, dpi = 300)
  
  list(AUROC = auroc, AUPRC = auprc, Brier = brier, CalibrationTable = calib_sum)
}

score_one_dataset <- function(case_parquet, locked, tag, base_dir) {
  case_parquet <- must_exist(case_parquet, paste0(tag, "_case_features"))
  
  DT <- as.data.table(arrow::read_parquet(case_parquet))
  cat("\n================ External validation on", tag, "================\n")
  cat("Input case_features dimensions:", nrow(DT), "rows x", ncol(DT), "cols\n")
  cat("Column names:\n"); print(names(DT))
  
  ## Standardize ID column
  if ("ID" %in% names(DT)) {
    # ok
  } else if ("caseid" %in% names(DT)) {
    setnames(DT, "caseid", "ID")
  } else {
    stop(sprintf("[%s] Could not find ID/caseid column.", tag))
  }
  
  if (!"DIKI" %in% names(DT)) {
    stop(sprintf("[%s] DIKI column missing.", tag))
  }
  
  coef_names <- names(locked$coefs)
  
  ## Impute missing features with 0
  miss_cols <- setdiff(coef_names, names(DT))
  if (length(miss_cols)) {
    message(sprintf("[%s] The following Locked features are missing in data; imputing 0: %s",
                    tag, paste(miss_cols, collapse = ", ")))
    for (mc in miss_cols) DT[, (mc) := 0.0]
  }
  
  ## Clean feature columns (NA -> 0, convert to numeric)
  DT[, (coef_names) := lapply(.SD, function(x) {
    x[is.na(x)] <- 0
    as.numeric(x)
  }), .SDcols = coef_names]
  
  X   <- as.matrix(DT[, ..coef_names])
  eta <- as.numeric(X %*% as.numeric(locked$coefs)) + locked$intercept
  prob <- plogis(eta)
  
  pred_out <- DT[, .(ID, DIKI)]
  pred_out[, `:=`(prob = prob, eta = eta)]
  
  out_pred <- file.path(base_dir, "predictions", paste0(tag, "_Predictions.csv"))
  readr::write_csv(pred_out, out_pred)
  
  eval_res <- evaluate_and_plot(
    idy    = pred_out[, .(DIKI)],
    pred   = pred_out$prob,
    out_png = file.path(base_dir, "Evaluation", paste0(tag, "_Calibration.png")),
    tag    = tag
  )
  
  ## Metrics summary table
  readr::write_csv(
    tibble(
      Dataset = tag,
      Metric  = c("AUROC","AUPRC","Brier"),
      Value   = c(eval_res$AUROC, eval_res$AUPRC, eval_res$Brier)
    ),
    file.path(base_dir, "Evaluation", paste0(tag, "_metrics_summary.csv"))
  )
  
  data.table::fwrite(
    eval_res$CalibrationTable,
    file.path(base_dir, "Evaluation", paste0(tag, "_calibration_table.csv"))
  )
  
  invisible(list(pred = pred_out, metrics = eval_res))
}

## ===================== Main Execution Workflow =====================

cat("=========== Module A: Locking Molecule-level CNI+APA Model ===========\n")
locked_spec <- lock_model_from_refit(REFIT_OR_FILE, LOCKED_JSON)

cat("\n=========== Module B: Generating case_features (FAERS & JADER) ===========\n")
cf_faers <- build_case_features_from_master(
  master_parquet = FAERS_MASTER_PARQUET,
  out_parquet    = OUT_CASE_FAERS,
  source_tag     = "FAERS"
)

cf_jader <- build_case_features_from_master(
  master_parquet = JADER_MASTER_PARQUET,
  out_parquet    = OUT_CASE_JADER,
  source_tag     = "JADER"
)

cat("\n=========== Module C: External Validation using LockedModel ===========\n")

## Use locked_spec generated in Module A directly
locked <- list(
  intercept = locked_spec$intercept,
  coefs     = locked_spec$coefficients
)

cat("LockedModel loaded from memory object 'locked_spec'.\n")
cat("Intercept =", locked$intercept, "\n")
cat("Coefficient names:\n"); print(names(locked$coefs))

faers_eval <- score_one_dataset(
  case_parquet = OUT_CASE_FAERS,
  locked       = locked,
  tag          = "FAERS_molecule_CNI_APA",
  base_dir     = BASE_EXTVAL_DIR
)

jader_eval <- score_one_dataset(
  case_parquet = OUT_CASE_JADER,
  locked       = locked,
  tag          = "JADER_molecule_CNI_APA",
  base_dir     = BASE_EXTVAL_DIR
)

## Consolidated Metrics Summary
metrics_all <- rbind(
  readr::read_csv(file.path(BASE_EXTVAL_DIR, "Evaluation", "FAERS_molecule_CNI_APA_metrics_summary.csv"), show_col_types = FALSE),
  readr::read_csv(file.path(BASE_EXTVAL_DIR, "Evaluation", "JADER_molecule_CNI_APA_metrics_summary.csv"), show_col_types = FALSE)
)
readr::write_csv(metrics_all, file.path(BASE_EXTVAL_DIR, "Evaluation", "metrics_summary_ALL.csv"))

run_manifest <- list(
  timestamp           = as.character(Sys.time()),
  base_dir            = BASE_EXTVAL_DIR,
  locked_model_spec   = LOCKED_JSON,
  faers_case_parquet  = OUT_CASE_FAERS,
  jader_case_parquet  = OUT_CASE_JADER,
  n_rows_scored = list(
    FAERS = nrow(faers_eval$pred),
    JADER = nrow(jader_eval$pred)
  ),
  features_used        = length(locked$coefs),
  locked_feature_names = names(locked$coefs)
)
jsonlite::write_json(run_manifest,
                     file.path(BASE_EXTVAL_DIR, "run_manifest.json"),
                     pretty = TRUE, auto_unbox = TRUE)

cat("\n✅ Pipeline complete. Output directory:", BASE_EXTVAL_DIR, "\n")
cat("- predictions/*_molecule_CNI_APA_Predictions.csv\n")
cat("- Evaluation/*_molecule_CNI_APA_* and metrics_summary_ALL.csv\n")
cat("- run_manifest.json\n")