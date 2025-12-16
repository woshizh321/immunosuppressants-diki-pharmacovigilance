#!/usr/bin/env Rscript
## ============================================================================
## 05_LASSO_Molecule_CNI_APA_withIndication_noRegion_noTTO_noSerious_noyear.R
##
## Objectives:
##   - Based on FAERS MASTER (2004–2024)
##   - Construct a molecule-level DIKI LASSO model exclusively within the 
##     CNI + APA exposure population.
##   - Covariates are consistent with the final drug-level model:
##       * Gender (sex)
##       * Age stratification (age_cat: <18, 18–39, 40–64, 65+)
##       * Indication group (indication_grp: Other / Transplantation / Autoimmune)
##       * Molecular exposure: CNI + APA molecules only
##   - Exclude year / region / TTO / serious as covariates
##   - Outputs:
##       * CV curves / coefficient paths / OOF / Brier / HL / DCA / SHAP-like / Refit, etc.
##
## Results Directory:
##   - D:/FAERS/Inhibitors/Derived_molecule_withIndi_noyear
## Temp Directory:
##   - D:/FAERS/Inhibitors/Derived_molecule_withIndi_noyear/Temp
##
## ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(Matrix)
  library(pROC)
  library(PRROC)
  library(broom)
  library(sandwich)
  library(lmtest)
  library(ResourceSelection)
  library(duckdb)
  library(DBI)
})

## ---------------------------- CONFIG ---------------------------------------

SEED           <- 20251121

MASTER_PARQUET <- "D:/FAERS/MASTER/FAERS_MASTER_FILE_2004-2024_with_serious.parquet"

DATA_PATH      <- "D:/FAERS/Inhibitors"
OUT_MAIN       <- file.path(DATA_PATH, "Derived_molecule_withIndi_noyear")
OUT_PATH       <- file.path(OUT_MAIN, "LASSO_path")
TEMP_DIR       <- file.path(OUT_MAIN, "Temp")

dir.create(OUT_MAIN, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_PATH, showWarnings = FALSE, recursive = TRUE)
dir.create(TEMP_DIR, showWarnings = FALSE, recursive = TRUE)

## Outcome and ID
TARGET_COL  <- "y"
ID_COL      <- "id_case"   # Case-level ID

## LASSO Parameters
K_FOLDS     <- 10
FAMILY      <- "binomial"
ALPHA       <- 1
N_LAMBDA    <- 100
CAL_BINS    <- 10

set.seed(SEED)

## ------------------------- Utility Functions --------------------------------

to_sparse_model_matrix <- function(dt, y_col, drop_cols = NULL) {
  stopifnot(y_col %in% names(dt))
  dt <- as.data.table(copy(dt))
  
  ## 1. Remove columns not used in the model (e.g., ID)
  if (!is.null(drop_cols)) {
    drop_cols <- intersect(drop_cols, names(dt))
    if (length(drop_cols) > 0) {
      dt[, (drop_cols) := NULL]
    }
  }
  
  ## 2. Extract y
  y <- as.integer(dt[[y_col]])
  dt[[y_col]] <- NULL
  
  ## 3. Perform complete.cases on independent variables only
  cc <- complete.cases(dt)
  n_all  <- nrow(dt)
  n_keep <- sum(cc)
  n_drop <- n_all - n_keep
  
  cat("to_sparse_model_matrix(): Total cases =", n_all,
      "; Complete cases =", n_keep,
      "; Cases dropped due to NA =", n_drop, "\n")
  
  if (n_keep <= 0) {
    stop("All rows contain NA in predictors; cannot construct model matrix.")
  }
  
  dt <- dt[cc]
  y  <- y[cc]
  
  ## 4. Convert character/factor -> factor
  for (nm in names(dt)) {
    if (is.character(dt[[nm]]) || is.factor(dt[[nm]])) {
      dt[[nm]] <- as.factor(dt[[nm]])
    }
  }
  
  ## 5. Construct sparse design matrix
  mm <- sparse.model.matrix(~ . - 1, data = dt)
  
  list(
    X          = mm,
    y          = y,
    feat_names = colnames(mm)
  )
}

cv_lambda_summary <- function(cvfit) {
  data.table(lambda = cvfit$lambda, cvm = cvfit$cvm, cvsd = cvfit$cvsd)
}

coef_paths_export <- function(cvfit, feat_names, out_prefix) {
  betas   <- as.matrix(cvfit$glmnet.fit$beta)
  lambdas <- cvfit$glmnet.fit$lambda
  stopifnot(ncol(betas) == length(lambdas))
  
  n_feat <- length(feat_names)
  dt_path <- data.table(
    feature    = rep(feat_names, times = length(lambdas)),
    lambda     = rep(lambdas,   each  = n_feat),
    log_lambda = rep(log(lambdas), each = n_feat),
    coef       = as.vector(betas)
  )
  fwrite(
    dt_path,
    file.path(OUT_PATH, sprintf("coef_paths_%s.csv", out_prefix))
  )
  
  for (s in c("lambda.min", "lambda.1se")) {
    b <- as.matrix(coef(cvfit, s = s))
    terms <- rownames(b)
    dt <- data.table(term = terms, coef = as.numeric(b))
    dt <- dt[term != "(Intercept)"]
    fwrite(
      dt,
      file.path(OUT_PATH,
                sprintf("coef_at_%s_%s.csv",
                        sub("\\.", "_", s), out_prefix))
    )
  }
}

kfold_oof_predict <- function(X, y, k = 10, seed = 1,
                              alpha = 1, family = "binomial") {
  set.seed(seed)
  n     <- length(y)
  folds <- sample(rep(1:k, length.out = n))
  oof   <- data.table(id = seq_len(n), fold = folds, y = y, p_hat = NA_real_)
  aucs  <- numeric(k)
  for (i in seq_len(k)) {
    tr_idx <- which(folds != i)
    te_idx <- which(folds == i)
    cvfit_i <- cv.glmnet(
      x = X[tr_idx, , drop = FALSE],
      y = y[tr_idx],
      family        = family,
      alpha         = alpha,
      type.measure = "deviance",
      nlambda      = N_LAMBDA
    )
    ph <- as.numeric(predict(
      cvfit_i,
      newx = X[te_idx, , drop = FALSE],
      s    = "lambda.1se",
      type = "response"
    ))
    oof$p_hat[te_idx] <- ph
    aucs[i] <- as.numeric(
      auc(pROC::roc(y[te_idx], ph, quiet = TRUE, direction = "<"))
    )
  }
  list(
    oof      = oof,
    fold_auc = aucs,
    mean_auc = mean(aucs, na.rm = TRUE)
  )
}

brier_score <- function(y, p) {
  mean((p - y)^2, na.rm = TRUE)
}

hl_test_result <- function(y, p, g = 10) {
  suppressWarnings(capture.output({
    ht <- hoslem.test(y, p, g = g)
  }))
  data.table(
    statistic = unname(ht$statistic),
    df        = unname(ht$parameter),
    p_value   = ht$p.value
  )
}

calibration_table <- function(y, p, bins = 10) {
  dt <- data.table(y = y, p = p)
  dt <- dt[is.finite(y) & is.finite(p)]
  if (nrow(dt) == 0) {
    stop("calibration_table(): no valid rows after filtering.")
  }
  
  if (sd(dt$p) == 0) {
    return(data.table(
      bin       = "[unique]",
      n         = nrow(dt),
      mean_pred = mean(dt$p),
      obs_rate  = mean(dt$y),
      sd_pred   = 0
    ))
  }
  
  q <- quantile(dt$p, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE)
  q <- unique(q)
  if (length(q) < 2) {
    return(data.table(
      bin       = "[all]",
      n         = nrow(dt),
      mean_pred = mean(dt$p),
      obs_rate  = mean(dt$y),
      sd_pred   = sd(dt$p)
    ))
  }
  
  dt[, bin := cut(p, breaks = q, include.lowest = TRUE, dig.lab = 6)]
  
  dt[, .(
    n         = .N,
    mean_pred = mean(p),
    obs_rate  = mean(y),
    sd_pred   = sd(p)
  ), by = bin]
}

decision_curve_table <- function(y, p,
                                 thresholds = seq(0.01, 0.99, by = 0.01)) {
  dt <- data.table(y = y, p = p)
  N  <- nrow(dt)
  rbindlist(lapply(thresholds, function(pt) {
    pred_pos <- as.integer(dt$p >= pt)
    TP <- sum(pred_pos == 1 & dt$y == 1)
    FP <- sum(pred_pos == 1 & dt$y == 0)
    NB <- TP / N - FP / N * (pt / (1 - pt))
    data.table(threshold = pt, net_benefit = NB)
  }))
}

shap_like_global_importance <- function(X, feat_names, coef_vec) {
  stopifnot(length(coef_vec) == length(feat_names))
  contrib_mean <- numeric(length(feat_names))
  for (j in seq_along(feat_names)) {
    xj <- X[, j]
    bj <- coef_vec[j]
    contrib_mean[j] <- mean(abs(as.numeric(xj) * bj))
  }
  data.table(feature = feat_names, shap_value = contrib_mean)
}

## ----- 16 DIKI PT codes (consistent with previous) -------------------------

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

## ---------------- CNI + APA Molecule Dictionary & Indication Definitions ---

## Retain only molecules involved in CNI + APA
## Each key is a "molecule-level" concept; value is a set of synonyms in UPPER(drugname)
molecule_dict <- list(
  ## CNI
  "Tacrolimus"   = c("TACROLIMUS", "FK506", "PROGRAF", "ADVAGRAF"),
  "Cyclosporine" = c("CYCLOSPORINE", "NEORAL", "SANDIMMUNE"),
  
  ## APA
  "Mycophenolate_mofetil" = c("MYCOPHENOLATE MOFETIL", "CELLCEPT"),
  "Mycophenolic_acid"      = c("MYCOPHENOLIC ACID", "MYFORTIC"),
  "Azathioprine"           = c("AZATHIOPRINE", "IMURAN")
)

mol_flags <- c(
  "flag_mol_tacrolimus",
  "flag_mol_cyclosporine",
  "flag_mol_mmf",
  "flag_mol_mpa",
  "flag_mol_azathioprine"
)

## 🔍 ADDED: Combine all CNI+APA synonyms into a SQL 'IN' clause to avoid pulling full database
mol_all_keywords <- unique(unlist(molecule_dict))
drug_filter_sql <- sprintf(
  "AND UPPER(drugname) IN (%s)",
  paste(sprintf("'%s'", mol_all_keywords), collapse = ",")
)

## Indication Keywords (consistent with previous)
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

## ============================================================================
## Step 1. Connect to DuckDB and inspect MASTER structure
## ============================================================================

stopifnot(file.exists(MASTER_PARQUET))

cat(">>> Step 1: Connecting to DuckDB and reading MASTER column names (LIMIT 0)...\n")

con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:", read_only = TRUE)
on.exit({
  try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
}, add = TRUE)

cols0 <- dbGetQuery(con,
                    sprintf("SELECT * FROM read_parquet('%s') LIMIT 0",
                            MASTER_PARQUET))
all_cols <- names(cols0)
cat("MASTER Column Count:", length(all_cols), "\n")

## Case ID Column
if ("id_case" %in% all_cols) {
  CASE_ID_COL <- "id_case"
} else if ("primaryid" %in% all_cols) {
  CASE_ID_COL <- "primaryid"
} else if ("caseid" %in% all_cols) {
  CASE_ID_COL <- "caseid"
} else {
  stop("id_case / primaryid / caseid not found in MASTER.")
}
cat("Case ID Column:", CASE_ID_COL, "\n")

## Drug Name Column
if (!"drugname" %in% all_cols) {
  stop("drugname column not found in MASTER.")
}

## Gender Column
if ("sex_raw" %in% all_cols) {
  SEX_COL_RAW <- "sex_raw"
} else if ("sex_std" %in% all_cols) {
  SEX_COL_RAW <- "sex_std"
} else if ("sex" %in% all_cols) {
  SEX_COL_RAW <- "sex"
} else {
  stop("sex_raw/sex_std/sex not found in MASTER.")
}
cat("Gender Column:", SEX_COL_RAW, "\n")

## Age Column (Continuous or Grouped)
AGE_YEARS_COL <- NA_character_
AGE_GROUP_COL <- NA_character_
if ("age_years_raw" %in% all_cols) {
  AGE_YEARS_COL <- "age_years_raw"
} else if ("age_years" %in% all_cols) {
  AGE_YEARS_COL <- "age_years"
} else if ("age" %in% all_cols) {
  AGE_YEARS_COL <- "age"
}
if ("age_group_raw" %in% all_cols) {
  AGE_GROUP_COL <- "age_group_raw"
} else if ("age_group" %in% all_cols) {
  AGE_GROUP_COL <- "age_group"
}

if (is.na(AGE_YEARS_COL) && is.na(AGE_GROUP_COL)) {
  stop("Age related columns not found in MASTER; cannot construct age categories.")
}
cat("Continuous Age Column:", AGE_YEARS_COL, "; Grouped Age Column:", AGE_GROUP_COL, "\n")

## Year Column (Description only, not included in model)
YEAR_COL <- NA_character_
if ("year_raw" %in% all_cols) {
  YEAR_COL <- "year_raw"
} else if ("year" %in% all_cols) {
  YEAR_COL <- "year"
}
cat("Year Column:", YEAR_COL, "\n")

## Severity Column
ser_cand <- grep("serious", all_cols, ignore.case = TRUE, value = TRUE)
if ("serious_flag" %in% all_cols) {
  SERIOUS_COL <- "serious_flag"
} else if (length(ser_cand) == 1) {
  SERIOUS_COL <- ser_cand[1]
} else if (length(ser_cand) == 0) {
  SERIOUS_COL <- NA_character_
  cat("serious* column not found in MASTER; unable to distinguish serious/non-serious (overall only).\n")
} else {
  stop("Multiple serious* columns detected; unable to uniquely identify severity variable.")
}
cat("Severity Column:", SERIOUS_COL, "\n")

## Indication Text
INDI_COL <- NA_character_
if ("indi_pt_norm_raw" %in% all_cols) {
  INDI_COL <- "indi_pt_norm_raw"
} else if ("indi_pt_norm" %in% all_cols) {
  INDI_COL <- "indi_pt_norm"
}
if (is.na(INDI_COL)) {
  stop("indi_pt_norm_raw / indi_pt_norm not found; cannot construct indication_grp.")
}
cat("Indication Text Column:", INDI_COL, "\n")

## pt_code Column (for constructing flag_diki)
if (!"pt_code" %in% all_cols) {
  stop("pt_code not found; cannot define flag_diki based on 16 DIKI PTs.")
}
cat("Constructing flag_diki dynamically using pt_code based on 16 DIKI PTs.\n")

## Role Code Column
if (!"role_cod" %in% all_cols) {
  stop("role_cod column not found in MASTER.")
}

## ============================================================================
## Step 2. Extract Row-level Records (PS/SS) using DuckDB, Build Exposure Info
## ============================================================================

cat(">>> Step 2: Extracting PS/SS row-level records from MASTER...\n")

## Required raw columns
needed <- c(CASE_ID_COL, "pt_code", "drugname", "role_cod",
            INDI_COL, SEX_COL_RAW)
if (!is.na(AGE_YEARS_COL)) needed <- c(needed, AGE_YEARS_COL)
if (!is.na(AGE_GROUP_COL)) needed <- c(needed, AGE_GROUP_COL)
if (!is.na(YEAR_COL))      needed <- c(needed, YEAR_COL)
if (!is.na(SERIOUS_COL))   needed <- c(needed, SERIOUS_COL)
needed <- unique(needed)

select_fields <- c(
  sprintf("%s AS id_case", CASE_ID_COL),
  "pt_code",
  sprintf("%s AS indi_pt_norm_raw", INDI_COL),
  sprintf("%s AS sex_raw", SEX_COL_RAW),
  "UPPER(drugname) AS drug_upper",
  "role_cod"
)

if (!is.na(AGE_YEARS_COL)) {
  select_fields <- c(select_fields,
                     sprintf("%s AS age_years_raw", AGE_YEARS_COL))
}
if (!is.na(AGE_GROUP_COL)) {
  select_fields <- c(select_fields,
                     sprintf("%s AS age_group_raw", AGE_GROUP_COL))
}
if (!is.na(YEAR_COL)) {
  select_fields <- c(select_fields,
                     sprintf("%s AS year_raw", YEAR_COL))
}
if (!is.na(SERIOUS_COL)) {
  select_fields <- c(select_fields,
                     sprintf("%s AS serious_raw", SERIOUS_COL))
}

sql_main <- sprintf("
  SELECT
    %s
  FROM read_parquet('%s')
  WHERE role_cod IN ('PS','SS')
  %s
",
                    paste(select_fields, collapse = ",\n    "),
                    MASTER_PARQUET,
                    drug_filter_sql
)

cat("Executing SQL (truncated display):\n")
cat(substr(sql_main, 1, 400), "...\n")

dt_raw <- as.data.table(dbGetQuery(con, sql_main))
cat("dt_raw Dimensions:", nrow(dt_raw), "rows x", ncol(dt_raw), "cols\n")

if (nrow(dt_raw) == 0L) {
  stop("No PS/SS records filtered from MASTER.")
}

## flag_diki
dt_raw[, flag_diki := as.integer(pt_code %in% diki_pt_codes)]
cat("Row-level flag_diki distribution:\n")
print(table(dt_raw$flag_diki, useNA = "ifany"))

## ============================================================================
## Step 3. Build Case-level Data: CNI+APA flags + flag_diki + covariates
## ============================================================================

cat(">>> Step 3.1: Building Case-level CNI+APA molecule exposure flags...\n")

dt_case <- dt_raw[, {
  out <- list()
  
  ## Case ID
  out[[ID_COL]] <- id_case[1]
  
  ## flag_diki: Case-level, 1 if any record is 1
  out[["flag_diki"]] <- as.integer(any(flag_diki == 1, na.rm = TRUE))
  
  ## Sex, Indication, Year, Age, etc.: Take the first non-NA
  out[["sex_raw"]]          <- sex_raw[which.max(!is.na(sex_raw))]
  out[["indi_pt_norm_raw"]] <- indi_pt_norm_raw[which.max(!is.na(indi_pt_norm_raw))]
  
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
  
  ## Molecular Exposure: 1 if any drug_upper record is in the dictionary
  du <- drug_upper
  
  out[["flag_mol_tacrolimus"]]   <- as.integer(any(du %in% molecule_dict[["Tacrolimus"]],   na.rm = TRUE))
  out[["flag_mol_cyclosporine"]] <- as.integer(any(du %in% molecule_dict[["Cyclosporine"]], na.rm = TRUE))
  out[["flag_mol_mmf"]]          <- as.integer(any(du %in% molecule_dict[["Mycophenolate_mofetil"]], na.rm = TRUE))
  out[["flag_mol_mpa"]]          <- as.integer(any(du %in% molecule_dict[["Mycophenolic_acid"]],     na.rm = TRUE))
  out[["flag_mol_azathioprine"]] <- as.integer(any(du %in% molecule_dict[["Azathioprine"]],          na.rm = TRUE))
  
  out
}, by = id_case]

cat("Case-level dt_case Dimensions:", nrow(dt_case), "rows x", ncol(dt_case), "\n")
cat("Case-level flag_diki distribution:\n")
print(table(dt_case$flag_diki, useNA = "ifany"))

## Exposure to at least one CNI/APA molecule
expo_any <- dt_case[, rowSums(.SD, na.rm = TRUE),
                    .SDcols = mol_flags]
dt_case <- dt_case[expo_any > 0]
cat("Cases exposed to at least one CNI/APA molecule:", nrow(dt_case), "\n")

## Severity: Description only, no filtering
if ("serious_flag" %in% names(dt_case)) {
  cat("serious_flag distribution (0=non-serious, 1=serious, NA=missing):\n")
  print(table(dt_case$serious_flag, useNA = "ifany"))
} else {
  cat("serious_flag not detected; model is 'overall' (severity not distinguished).\n")
}

## ============================================================================
## Step 4. Construct LASSO input data raw_lasso_dt (including indication_grp)
## ============================================================================

cat(">>> Step 4: Constructing LASSO input data (including indication_grp)...\n")

raw_lasso_dt <- copy(dt_case)

## Outcome Column
raw_lasso_dt[, (TARGET_COL) := as.integer(flag_diki)]

## Gender Factor
raw_lasso_dt[, sex := factor(sex_raw)]

## Age Stratification
if ("age_years_raw" %in% names(raw_lasso_dt)) {
  raw_lasso_dt[, age_years_raw := suppressWarnings(as.numeric(age_years_raw))]
  raw_lasso_dt[is.na(age_years_raw) | age_years_raw <= 0,
               age_cat := NA_character_]
  raw_lasso_dt[!is.na(age_years_raw) & age_years_raw < 18,
               age_cat := "<18"]
  raw_lasso_dt[!is.na(age_years_raw) & age_years_raw >= 18 & age_years_raw <= 39,
               age_cat := "18-39"]
  raw_lasso_dt[!is.na(age_years_raw) & age_years_raw >= 40 & age_years_raw <= 64,
               age_cat := "40-64"]
  raw_lasso_dt[!is.na(age_years_raw) & age_years_raw >= 65,
               age_cat := "65+"]
  
  raw_lasso_dt[, age_cat := factor(
    age_cat,
    levels = c("<18", "18-39", "40-64", "65+")
  )]
} else if ("age_group_raw" %in% names(raw_lasso_dt)) {
  raw_lasso_dt[, age_cat := factor(age_group_raw)]
} else {
  stop("Missing age_years_raw or age_group_raw; cannot construct age stratification.")
}

cat("age_cat distribution:\n")
print(table(raw_lasso_dt$age_cat, useNA = "ifany"))

## indication_grp Three-category
if (!"indi_pt_norm_raw" %in% names(raw_lasso_dt)) {
  stop("indi_pt_norm_raw not found in raw_lasso_dt; cannot construct indication_grp.")
}

raw_lasso_dt[, indi_up := toupper(trimws(indi_pt_norm_raw))]
raw_lasso_dt[, indication_grp := "Other"]

raw_lasso_dt[grepl(pat_tx, indi_up),
             indication_grp := "Transplantation"]

raw_lasso_dt[indication_grp == "Other" &
               grepl(pat_ai, indi_up),
             indication_grp := "Autoimmune diseases"]

raw_lasso_dt[, indication_grp := factor(
  indication_grp,
  levels = c("Other", "Transplantation", "Autoimmune diseases")
)]

cat("indication_grp distribution:\n")
print(table(raw_lasso_dt$indication_grp, useNA = "ifany"))

## Force molecule flags to 0/1
for (cl in mol_flags) {
  raw_lasso_dt[, (cl) :=
                 as.integer(get(cl) %in% c(1, "1", "Y", "YES", TRUE))]
}
cat("molecule-level flag variables:\n")
print(mol_flags)

## Remove auxiliary columns not for modeling (year_raw is excluded)
drop_aux <- c(
  "flag_diki",
  "sex_raw",
  "age_years_raw",
  "age_group_raw",
  "year_raw",
  "indi_pt_norm_raw",
  "indi_up",
  "serious_raw",
  "serious_flag",
  "drug_upper"   # Already gone at case-level, but just in case
)
drop_aux <- intersect(drop_aux, names(raw_lasso_dt))
if (length(drop_aux) > 0) {
  cat("Removing auxiliary columns from predictors:\n")
  print(drop_aux)
  raw_lasso_dt[, (drop_aux) := NULL]
}

## Outcome Distribution QC
cat(">>> LASSO input data outcome distribution (y):\n")
print(table(raw_lasso_dt[[TARGET_COL]], useNA = "ifany"))
if (length(unique(raw_lasso_dt[[TARGET_COL]])) < 2) {
  stop("Outcome y has only one unique value; cannot perform modeling.")
}

## ============================================================================
## Step 5. Build Sparse Model Matrix, Run LASSO + OOF + Calibration + Refit
## ============================================================================

cat(">>> Step 5: Building sparse model matrix and running LASSO...\n")

prep <- to_sparse_model_matrix(
  dt        = raw_lasso_dt,
  y_col     = TARGET_COL,
  drop_cols = ID_COL   # id_case does not enter model
)
X     <- prep$X
y     <- prep$y
FEATS <- prep$feat_names

cat("Sparse matrix Dimensions:", dim(X)[1], "x", dim(X)[2], "\n")

set.seed(SEED)
cvfit <- cv.glmnet(
  x = X, y = y,
  family        = FAMILY,
  alpha         = ALPHA,
  type.measure = "deviance",
  nlambda      = N_LAMBDA
)

## CV Curve
cv_curve <- cv_lambda_summary(cvfit)
fwrite(
  cv_curve,
  file.path(OUT_MAIN,
            "lasso_molecule_CNI_APA_cv_curve_withIndication_noRegion_noTTO_noSerious_noyear.csv")
)

## Coefficient Path
coef_paths_export(
  cvfit,
  FEATS,
  out_prefix = "lasso_molecule_CNI_APA_cvfit_withIndication_noRegion_noTTO_noSerious_noyear"
)

## Model Lock Object
model_lock <- list(
  cvfit             = cvfit,
  seed              = SEED,
  family            = FAMILY,
  alpha             = ALPHA,
  nlambda           = N_LAMBDA,
  features          = FEATS,
  target            = TARGET_COL,
  id_col            = ID_COL,
  molecule_dict     = molecule_dict,
  mol_flags         = mol_flags,
  indication_levels = levels(raw_lasso_dt$indication_grp)
)
saveRDS(
  model_lock,
  file = file.path(
    OUT_MAIN,
    "LASSO_Molecule_CNI_APA_Model_withIndication_noRegion_noTTO_noSerious_noyear.rds"
  )
)
cat("Model lock file saved to:\n  ",
    file.path(OUT_MAIN,
              "LASSO_Molecule_CNI_APA_Model_withIndication_noRegion_noTTO_noSerious_noyear.rds"),
    "\n")

## --------------------------- OOF Predictions and Metrics -------------------

cat(">>> Step 6: Calculating OOF predictions and performance metrics...\n")

oof <- kfold_oof_predict(
  X, y,
  k      = K_FOLDS,
  seed   = SEED,
  alpha  = ALPHA,
  family = FAMILY
)

fwrite(
  oof$oof,
  file.path(OUT_MAIN,
            "cv_oof_predictions_molecule_CNI_APA_withIndication_noRegion_noTTO_noSerious_noyear.csv")
)
fwrite(
  data.table(
    fold     = 1:K_FOLDS,
    AUC      = oof$fold_auc,
    mean_AUC = oof$mean_auc
  ),
  file.path(OUT_MAIN,
            "OOF_auc_summary_molecule_CNI_APA_withIndication_noRegion_noTTO_noSerious_noyear.csv")
)

fwrite(
  data.table(
    Brier = brier_score(oof$oof$y, oof$oof$p_hat)
  ),
  file.path(OUT_MAIN,
            "brier_score_molecule_CNI_APA_withIndication_noRegion_noTTO_noSerious_noyear.csv")
)

cal_tbl <- calibration_table(oof$oof$y, oof$oof$p_hat, bins = CAL_BINS)
fwrite(
  cal_tbl,
  file.path(OUT_MAIN,
            "calibration_table_oof_molecule_CNI_APA_withIndication_noRegion_noTTO_noSerious_noyear.csv")
)

fwrite(
  hl_test_result(oof$oof$y, oof$oof$p_hat, g = CAL_BINS),
  file.path(OUT_MAIN,
            "hl_test_molecule_CNI_APA_withIndication_noRegion_noTTO_noSerious_noyear.csv")
)

pr <- PRROC::pr.curve(
  scores.class0 = oof$oof$p_hat[oof$oof$y == 1],
  scores.class1 = oof$oof$p_hat[oof$oof$y == 0]
)
fwrite(
  data.table(AUPRC = pr$auc.integral),
  file.path(OUT_MAIN,
            "auprc_molecule_CNI_APA_withIndication_noRegion_noTTO_noSerious_noyear.csv")
)

fwrite(
  decision_curve_table(oof$oof$y, oof$oof$p_hat),
  file.path(OUT_MAIN,
            "dca_table_molecule_CNI_APA_withIndication_noRegion_noTTO_noSerious_noyear.csv")
)

## --------------------------- SHAP-like Global Importance -------------------

cat(">>> Step 7: Calculating SHAP-like feature importance...\n")

b_1se <- as.matrix(coef(cvfit, s = "lambda.1se"))
b_1se <- b_1se[rownames(b_1se) != "(Intercept)", , drop = FALSE]
coef_vec <- as.numeric(b_1se)
names(coef_vec) <- rownames(b_1se)

idx  <- match(FEATS, names(coef_vec))
keep <- which(!is.na(idx))

imp  <- shap_like_global_importance(
  X[, keep, drop = FALSE],
  FEATS[keep],
  coef_vec[idx[keep]]
)
imp[, beta := coef_vec[idx[keep]]]
imp[, OR   := exp(beta)]

fwrite(
  imp,
  file.path(OUT_MAIN,
            "lasso_shap_like_importance_molecule_CNI_APA_withIndication_noRegion_noTTO_noSerious_noyear.csv")
)

## --------------------------- Refit (Standard + Robust) ---------------------

cat(">>> Step 8: Unpenalized logit Refit using lambda.1se non-zero variables...\n")

sel <- names(coef_vec)[abs(coef_vec) > 0]
sel <- sel[!is.na(sel) & nzchar(sel)]
cat("Count of non-zero features under lambda.1se:", length(sel), "\n")
if (length(sel) == 0) stop("No non-zero features found under lambda.1se.")

df_refit <- as.data.frame(as.matrix(X[, sel, drop = FALSE]))
colnames(df_refit) <- make.names(colnames(df_refit), unique = TRUE)
df_refit[[TARGET_COL]] <- y

form <- as.formula(paste(TARGET_COL, "~ ."))
cat("Refit formula used:\n"); print(form)

## Standard Logistic Refit
fit_std  <- glm(form, data = df_refit, family = binomial(link = "logit"))
tidy_std <- broom::tidy(fit_std, conf.int = TRUE, exponentiate = TRUE)

std_out  <- tidy_std[, c("term","estimate","conf.low","conf.high","p.value")]
setDT(std_out)
setnames(
  std_out,
  c("estimate","conf.low","conf.high","p.value"),
  c("OR","CI_low","CI_high","p_value")
)

fwrite(
  std_out,
  file.path(OUT_MAIN,
            "lasso_refit_logistic_results_molecule_CNI_APA_withIndication_noRegion_noTTO_noSerious_noyear.csv")
)

tab3 <- std_out[, .(
  Variable               = term,
  Group                  = NA_character_,
  `Adjusted OR (95% CI)` = sprintf("%.2f [%.2f–%.2f]", OR, CI_low, CI_high),
  `p value`              = sprintf("%.3g", p_value)
)]
fwrite(
  tab3,
  file.path(OUT_MAIN,
            "Table3_Refit_Logistic_OR_95CI_molecule_CNI_APA_withIndication_noRegion_noTTO_noSerious_noyear.csv")
)

## Sandwich Robust Refit
vc <- sandwich::vcovHC(fit_std, type = "HC0")
ct <- lmtest::coeftest(fit_std, vcov = vc)
rb <- data.table(
  term    = rownames(ct),
  beta    = ct[, "Estimate"],
  se      = ct[, "Std. Error"],
  z       = ct[, "z value"],
  p_value = ct[, "Pr(>|z|)"]
)
rb[, `:=`(
  OR      = exp(beta),
  CI_low  = exp(beta - 1.96 * se),
  CI_high = exp(beta + 1.96 * se)
)]
rb <- rb[term != "(Intercept)"]

fwrite(
  rb,
  file.path(OUT_MAIN,
            "lasso_refit_logistic_robust_molecule_CNI_APA_withIndication_noRegion_noTTO_noSerious_noyear.csv")
)

cat(">>> CNI + APA molecule-level + indication + noRegion + noTTO + noSerious-as-covariate LASSO pipeline complete.\n",
    "All outputs located in:", OUT_MAIN, " and ", OUT_PATH, "\n")

## ============================================================================
## Step 9. Save key objects to RDS (using Temp directory)
## ============================================================================

## 1. LASSO Input Data
saveRDS(
  raw_lasso_dt,
  file = file.path(TEMP_DIR,
                   "raw_lasso_dt_molecule_CNI_APA_withIndication_noRegion_noTTO_noSerious_noyear.rds")
)

## 2. Design Matrix Object
design_obj <- list(
  X     = X,
  y     = y,
  feats = FEATS
)
saveRDS(
  design_obj,
  file = file.path(TEMP_DIR,
                   "design_matrix_molecule_CNI_APA_withIndication_noRegion_noTTO_noSerious_noyear.rds")
)

## 3. cv.glmnet fit object
saveRDS(
  cvfit,
  file = file.path(TEMP_DIR,
                   "cvfit_molecule_CNI_APA_withIndication_noRegion_noTTO_noSerious_noyear.rds")
)

## 4. Model Lock File
saveRDS(
  model_lock,
  file = file.path(TEMP_DIR,
                   "model_lock_molecule_CNI_APA_withIndication_noRegion_noTTO_noSerious_noyear.rds")
)

## 5. OOF Results
saveRDS(
  oof,
  file = file.path(TEMP_DIR,
                   "oof_results_molecule_CNI_APA_withIndication_noRegion_noTTO_noSerious_noyear.rds")
)

## 6. Refit Logistic Model
refit_obj <- list(
  fit_std   = fit_std,
  tidy_std  = std_out,
  robust_rb = rb
)
saveRDS(
  refit_obj,
  file = file.path(TEMP_DIR,
                   "refit_logit_molecule_CNI_APA_withIndication_noRegion_noTTO_noSerious_noyear.rds")
)

## 7. SHAP-like Results
shap_obj <- list(
  shap_orig_lambda1se = imp
)
saveRDS(
  shap_obj,
  file = file.path(TEMP_DIR,
                   "shap_like_molecule_CNI_APA_withIndication_noRegion_noTTO_noSerious_noyear.rds")
)

cat("Key LASSO objects saved to Temp directory:\n  ", TEMP_DIR, "\n")