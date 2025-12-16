#!/usr/bin/env Rscript
## ============================================================================
## 04_LASSO_Drug_withIndication_noRegion_noTTO_noSerious_from_FAERS.R
## Objectives:
##   - Build drug-class level case data from FAERS MASTER (parquet).
##   - Drug classes based on 'drug_class' column in drug_class_map.csv:
##       CNI, mTOR Inhibitors, Antiproliferative Agents, Corticosteroids, Biologics
##   - Include indications (3 categories: Transplantation / Autoimmune / Other).
##   - Exclude region / TTO / serious as covariates (seriousness used for description only).
##   - Perform LASSO + OOF (Out-of-Fold) + Calibration + Refit + SHAP-like importance.
## Output Directory:
##   D:/FAERS/Inhibitors/Derived_drug_withIndi
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

## ---------------------------- CONFIG ----------------------------

SEED        <- 20251121

MASTER_PARQUET  <- "D:/FAERS/MASTER/FAERS_MASTER_FILE_2004-2024_with_serious.parquet"
DRUG_CLASS_FILE <- "D:/FAERS/MASTER/drug_class_map.csv"

DATA_PATH   <- "D:/FAERS/Inhibitors"
OUT_MAIN    <- file.path(DATA_PATH, "Derived_drug_withIndi")
OUT_PATH    <- file.path(OUT_MAIN, "LASSO_path")
dir.create(OUT_MAIN, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_PATH, showWarnings = FALSE, recursive = TRUE)

TARGET_COL  <- "y"
ID_COL      <- "id_case"       # Case-level ID
K_FOLDS     <- 10
FAMILY      <- "binomial"
ALPHA       <- 1
N_LAMBDA    <- 100
CAL_BINS    <- 10

set.seed(SEED)

## 16 DIKI PT codes (Drug-Induced Kidney Injury)
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

## ------------------------- UTILITIES ----------------------------

to_sparse_model_matrix <- function(dt, y_col, drop_cols = NULL) {
  stopifnot(y_col %in% names(dt))
  dt <- as.data.table(copy(dt))
  
  ## 1. Remove columns not used in modeling (e.g., id_case)
  if (!is.null(drop_cols)) {
    drop_cols <- intersect(drop_cols, names(dt))
    if (length(drop_cols) > 0) {
      dt[, (drop_cols) := NULL]
    }
  }
  
  ## 2. Extract y
  y <- as.integer(dt[[y_col]])
  dt[[y_col]] <- NULL
  
  ## 3. Independent variables complete.cases
  cc <- complete.cases(dt)
  n_all  <- nrow(dt)
  n_keep <- sum(cc)
  n_drop <- n_all - n_keep
  
  cat("to_sparse_model_matrix(): Total cases =", n_all,
      "; Complete cases =", n_keep,
      "; Cases excluded due to NA =", n_drop, "\n")
  
  if (n_keep <= 0) {
    stop("All rows contain NA in independent variables. Cannot construct model matrix.")
  }
  
  dt <- dt[cc]
  y  <- y[cc]
  
  ## 4. Convert Character/Factor to Factor
  for (nm in names(dt)) {
    if (is.character(dt[[nm]]) || is.factor(dt[[nm]])) {
      dt[[nm]] <- as.factor(dt[[nm]])
    }
  }
  
  ## 5. Sparse design matrix
  mm <- sparse.model.matrix(~ . - 1, data = dt)
  
  list(
    X          = mm,
    y          = y,
    feat_names = colnames(mm)
  )
}

cv_lambda_summary <- function(cvfit) {
  data.table(lambda = cvfit$lambda,
             cvm    = cvfit$cvm,
             cvsd   = cvfit$cvsd)
}

coef_paths_export <- function(cvfit, feat_names, out_prefix) {
  betas   <- as.matrix(cvfit$glmnet.fit$beta)
  lambdas <- cvfit$glmnet.fit$lambda
  stopifnot(ncol(betas) == length(lambdas))
  dt_path <- data.table(
    feature    = rep(feat_names, times = length(lambdas)),
    lambda     = rep(lambdas, each = length(feat_names)),
    log_lambda = rep(log(lambdas), each = length(feat_names)),
    coef       = as.vector(betas)
  )
  fwrite(dt_path,
         file.path(OUT_PATH,
                   sprintf("coef_paths_%s.csv", out_prefix)))
  
  for (s in c("lambda.min", "lambda.1se")) {
    b <- as.matrix(coef(cvfit, s = s))
    terms <- rownames(b)
    dt <- data.table(term = terms, coef = as.numeric(b))
    dt <- dt[term != "(Intercept)"]
    fwrite(dt,
           file.path(OUT_PATH,
                     sprintf("coef_at_%s_%s.csv",
                             sub("\\.", "_", s), out_prefix)))
  }
}

kfold_oof_predict <- function(X, y, k = 10, seed = 1,
                              alpha = 1, family = "binomial") {
  set.seed(seed)
  n     <- length(y)
  folds <- sample(rep(1:k, length.out = n))
  oof   <- data.table(id = seq_len(n), fold = folds, y = y, p_hat = NA_real_)
  aucs  <- numeric(k)
  for (i in 1:k) {
    cat("  [OOF] fold", i, "of", k, "...\n")
    tr_idx <- which(folds != i)
    te_idx <- which(folds == i)
    cvfit_i <- cv.glmnet(
      x = X[tr_idx, , drop = FALSE],
      y = y[tr_idx],
      family = family,
      alpha  = alpha,
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

brier_score <- function(y, p) mean((p - y)^2, na.rm = TRUE)

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
  dt[, bin := cut(
    p,
    breaks = quantile(p,
                      probs = seq(0, 1, length.out = bins + 1),
                      na.rm = TRUE),
    include.lowest = TRUE, dig.lab = 6
  )]
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

## ------------------------ Step 1. Define Drug Classes ----------------------------

stopifnot(file.exists(MASTER_PARQUET))
stopifnot(file.exists(DRUG_CLASS_FILE))

cat(">>> Step 1: Reading drug_class_map (Drug Class Definitions)...\n")
drug_map <- fread(DRUG_CLASS_FILE)
cat("drug_class_map column names:\n"); print(names(drug_map))

if (!all(c("drug_u", "drug_class") %in% names(drug_map))) {
  stop("drug_class_map.csv must contain 'drug_u' and 'drug_class' columns.")
}

## Retain only the 5 target classes
target_classes <- c("CNI",
                    "mTOR Inhibitors",
                    "Antiproliferative Agents",
                    "Corticosteroids",
                    "Biologics")

drug_map <- drug_map[drug_class %in% target_classes]
if (nrow(drug_map) == 0L) {
  stop("No matching target drug_classes found in drug_class_map.csv.")
}

drug_map[, drug_upper := toupper(drug_u)]
drug_map[, class      := as.character(drug_class)]

cat("Target drug classes and examples:\n")
print(unique(drug_map[, .(class)])[order(class)])

## Build case-level flag names for the 5 classes
class_levels <- target_classes
flag_names   <- paste0(
  "flag_",
  gsub("[^A-Za-z0-9]+", "_",
       tolower(class_levels))
)
drug_class_dict <- data.table(
  class = class_levels,
  flag  = flag_names
)
cat("Case-level drug class flags:\n")
print(drug_class_dict)

## ------------------------ Step 2. Extract from MASTER via DuckDB -----------------

cat(">>> Step 2: Connecting to DuckDB and inspecting MASTER columns...\n")
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:", read_only = TRUE)

cols0 <- dbGetQuery(con, sprintf(
  "SELECT * FROM read_parquet('%s') LIMIT 0", MASTER_PARQUET
))
all_cols <- names(cols0)
cat("Preview of MASTER column names (first 40):\n")
print(head(all_cols, 40))

## Case ID column
if ("primaryid" %in% all_cols) {
  CASE_ID_COL <- "primaryid"
} else if ("caseid" %in% all_cols) {
  CASE_ID_COL <- "caseid"
} else {
  stop("primaryid or caseid not found in MASTER.")
}
cat("MASTER Case ID column:", CASE_ID_COL, "\n")

## pt_code
if (!"pt_code" %in% all_cols) {
  stop("pt_code column not found in MASTER. Cannot define flag_diki.")
}

## Indication raw column: assumes indi_pt_norm_raw or equivalent exists
if ("indi_pt_norm" %in% all_cols) {
  INDI_COL <- "indi_pt_norm"
} else if ("indi_pt" %in% all_cols) {
  INDI_COL <- "indi_pt"
} else {
  stop("indi_pt_norm or indi_pt not found in MASTER. Cannot construct indication.")
}
cat("MASTER Indication column:", INDI_COL, "\n")

## Gender column
if ("sex_std" %in% all_cols) {
  SEX_COL_RAW <- "sex_std"
} else if ("sex" %in% all_cols) {
  SEX_COL_RAW <- "sex"
} else {
  stop("sex_std or sex not found in MASTER.")
}
cat("MASTER Sex column:", SEX_COL_RAW, "\n")

## Age
if ("age_years" %in% all_cols) {
  AGE_YEARS_COL_RAW <- "age_years"
  AGE_GROUP_COL_RAW <- NA_character_
  cat("Using continuous age column: age_years.\n")
} else if ("age_group" %in% all_cols) {
  AGE_YEARS_COL_RAW <- NA_character_
  AGE_GROUP_COL_RAW <- "age_group"
  cat("Using pre-defined age group column: age_group.\n")
} else {
  stop("age_years or age_group not found in MASTER.")
}

## Year
if ("year" %in% all_cols) {
  YEAR_COL_RAW <- "year"
} else {
  YEAR_COL_RAW <- NA_character_
  cat("year column not found in MASTER. Year variable will be excluded.\n")
}

## Seriousness (For description only, not for model input)
ser_cand <- grep("serious", all_cols, ignore.case = TRUE, value = TRUE)
if ("serious_flag" %in% all_cols) {
  SERIOUS_COL_RAW <- "serious_flag"
} else if (length(ser_cand) == 1) {
  SERIOUS_COL_RAW <- ser_cand[1]
} else if (length(ser_cand) == 0) {
  SERIOUS_COL_RAW <- NA_character_
  cat("Seriousness column not found in MASTER. serious_flag will not be constructed.\n")
} else {
  stop("Multiple serious* columns detected in MASTER. Cannot determine unique variable.")
}
cat("MASTER Seriousness raw column:", SERIOUS_COL_RAW, "\n")

## drugname & role_cod
if (!"drugname" %in% all_cols) {
  stop("drugname column not found in MASTER. Cannot match with drug_class_map.")
}
if (!"role_cod" %in% all_cols) {
  stop("role_cod column not found in MASTER. Cannot filter PS/SS records.")
}

## Write drug_map to DuckDB
dbWriteTable(
  con, "drug_class_map_tmp",
  drug_map[, .(drug_upper, class)],
  temporary = TRUE, overwrite = TRUE
)

## Build SELECT fields
select_fields <- c(
  sprintf("%s AS id_case", CASE_ID_COL),
  "pt_code",
  sprintf("%s AS indi_pt_norm_raw", INDI_COL),
  sprintf("%s AS sex_raw", SEX_COL_RAW),
  "UPPER(drugname) AS drug_upper",
  "role_cod"
)

if (!is.na(AGE_YEARS_COL_RAW)) {
  select_fields <- c(select_fields,
                     sprintf("%s AS age_years_raw", AGE_YEARS_COL_RAW))
}
if (!is.na(AGE_GROUP_COL_RAW)) {
  select_fields <- c(select_fields,
                     sprintf("%s AS age_group_raw", AGE_GROUP_COL_RAW))
}
if (!is.na(YEAR_COL_RAW)) {
  select_fields <- c(select_fields,
                     sprintf("%s AS year_raw", YEAR_COL_RAW))
}
if (!is.na(SERIOUS_COL_RAW)) {
  select_fields <- c(select_fields,
                     sprintf("%s AS serious_raw", SERIOUS_COL_RAW))
}

sql_main <- sprintf("
  WITH base AS (
    SELECT %s
    FROM read_parquet('%s')
    WHERE role_cod IN ('PS','SS')
  )
  SELECT
    b.*,
    dcm.class AS drug_class
  FROM base b
  JOIN drug_class_map_tmp dcm
    ON b.drug_upper = dcm.drug_upper
", paste(select_fields, collapse = ",\n        "),
                    MASTER_PARQUET)

cat("Executing SQL (Truncated view):\n")
cat(substr(sql_main, 1, 400), "...\n")

dt_raw <- as.data.table(dbGetQuery(con, sql_main))
cat("dt_raw rows:", nrow(dt_raw), " columns:", ncol(dt_raw), "\n")

if (nrow(dt_raw) == 0L) {
  stop("No PS/SS records matching drug class definitions found in MASTER.")
}

## Construct row-level flag_diki based on pt_code
dt_raw[, flag_diki := as.integer(pt_code %in% diki_pt_codes)]
cat("Row-level flag_diki distribution:\n")
print(table(dt_raw$flag_diki, useNA = "ifany"))

## ------------------------ Step 3. Aggregate to Case-level (dt_case) --------------------

cat(">>> Step 3: Aggregating to case-level (one row per id_case)...\n")

dt_case <- dt_raw[, {
  out <- list()
  
  out[[ID_COL]] <- id_case[1]
  
  ## sex_raw: first non-NA
  if (any(!is.na(sex_raw))) {
    out[["sex_raw"]] <- sex_raw[which.max(!is.na(sex_raw))]
  } else {
    out[["sex_raw"]] <- NA_character_
  }
  
  ## Age
  if ("age_years_raw" %in% names(.SD)) {
    age_tmp <- age_years_raw[which.max(!is.na(age_years_raw))]
    out[["age_years_raw"]] <- age_tmp
  }
  if ("age_group_raw" %in% names(.SD)) {
    out[["age_group_raw"]] <- age_group_raw[which.max(!is.na(age_group_raw))]
  }
  
  ## Year
  if ("year_raw" %in% names(.SD)) {
    out[["year_raw"]] <- year_raw[which.max(!is.na(year_raw))]
  }
  
  ## Indication text (single representative per case)
  if ("indi_pt_norm_raw" %in% names(.SD)) {
    out[["indi_pt_norm_raw"]] <- indi_pt_norm_raw[which.max(!is.na(indi_pt_norm_raw))]
  }
  
  ## Seriousness: Flag 1 if any record is serious
  if ("serious_raw" %in% names(.SD)) {
    out[["serious_raw"]] <- max(
      as.integer(serious_raw %in% c(1, "1", "Y", "YES", TRUE)),
      na.rm = TRUE
    )
  }
  
  ## Drug class flags (5 classes)
  for (i in seq_along(class_levels)) {
    cls  <- class_levels[i]
    vnm  <- flag_names[i]
    out[[vnm]] <- as.integer(any(drug_class == cls, na.rm = TRUE))
  }
  
  ## Case-level DIKI outcome
  out[["flag_diki"]] <- as.integer(any(flag_diki == 1, na.rm = TRUE))
  
  out
}, by = id_case]

cat("dt_case rows:", nrow(dt_case), " columns:", ncol(dt_case), "\n")
cat("Case-level DIKI distribution:\n")
print(table(dt_case$flag_diki, useNA = "ifany"))

if (length(unique(dt_case$flag_diki)) < 2) {
  stop("Case-level flag_diki has only one unique value. Modeling impossible.")
}

## ------------------------ Step 4. Build LASSO input (raw_lasso_dt) -----------

cat(">>> Step 4: Constructing LASSO input data (with indication, no region/TTO/serious)...\n")

raw_lasso_dt <- copy(dt_case)

## Outcome y
raw_lasso_dt[, (TARGET_COL) := flag_diki]

## Gender factor
raw_lasso_dt[, sex := factor(sex_raw)]

## Age Categorization
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
  stop("Missing age_years_raw or age_group_raw. Cannot build age categories.")
}

cat("age_cat distribution:\n")
print(table(raw_lasso_dt$age_cat, useNA = "ifany"))

## Year
if ("year_raw" %in% names(raw_lasso_dt)) {
  raw_lasso_dt[, year := factor(year_raw)]
}

## Seriousness: Kept for description only
if ("serious_raw" %in% names(raw_lasso_dt)) {
  raw_lasso_dt[, serious_flag :=
                 as.integer(serious_raw %in% c(1, "1", "Y", "YES", TRUE))]
  cat("Seriousness (serious_flag) distribution (Descriptive only):\n")
  print(table(raw_lasso_dt$serious_flag, useNA = "ifany"))
}

## Build indication_grp (3 Categories)
cat(">>> Building indication_grp (Other / Transplantation / Autoimmune diseases)...\n")
if (!"indi_pt_norm_raw" %in% names(raw_lasso_dt)) {
  stop("indi_pt_norm_raw not found in raw_lasso_dt. Cannot construct indication_grp.")
}
raw_lasso_dt[, indi_up := toupper(trimws(indi_pt_norm_raw))]

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
pat_tx <- paste0("(", paste(kw_transplant,  collapse = "|"), ")")
pat_ai <- paste0("(", paste(kw_autoimmune, collapse = "|"), ")")

raw_lasso_dt[, indication_grp := "Other"]
raw_lasso_dt[grepl(pat_tx, indi_up), indication_grp := "Transplantation"]
raw_lasso_dt[indication_grp == "Other" &
               grepl(pat_ai, indi_up),
             indication_grp := "Autoimmune diseases"]

raw_lasso_dt[, indication_grp := factor(
  indication_grp,
  levels = c("Other", "Transplantation", "Autoimmune diseases")
)]

cat("indication_grp distribution:\n")
print(table(raw_lasso_dt$indication_grp, useNA = "ifany"))

## Drug class flags (0/1)
drug_flags <- drug_class_dict$flag
cat("Drug class flag variables:\n"); print(drug_flags)
for (cl in drug_flags) {
  raw_lasso_dt[, (cl) :=
                 as.integer(get(cl) %in% c(1, "1", "Y", "YES", TRUE))]
}

## Remove auxiliary columns not included in modeling
drop_aux <- c(
  "flag_diki",
  "sex_raw",
  "age_years_raw",
  "age_group_raw",
  "year_raw",
  "indi_pt_norm_raw",
  "indi_up",
  "serious_raw",
  "serious_flag"
)
drop_aux <- intersect(drop_aux, names(raw_lasso_dt))
if (length(drop_aux) > 0) {
  cat("Removing auxiliary columns from independent variables:\n")
  print(drop_aux)
  raw_lasso_dt[, (drop_aux) := NULL]
}

## Outcome Distribution QC
cat(">>> LASSO input data outcome distribution (y):\n")
print(table(raw_lasso_dt[[TARGET_COL]], useNA = "ifany"))
if (length(unique(raw_lasso_dt[[TARGET_COL]])) < 2) {
  stop("LASSO input y has only one unique value. Modeling impossible.")
}

## ------------------------ Step 5. LASSO + OOF + Refit + SHAP ---------------

cat(">>> Step 5: Building sparse model matrix and running LASSO...\n")

prep <- to_sparse_model_matrix(
  dt        = raw_lasso_dt,
  y_col     = TARGET_COL,
  drop_cols = ID_COL
)
X     <- prep$X
y     <- prep$y
FEATS <- prep$feat_names

cat("Sparse Matrix Dimensions:", dim(X)[1], "x", dim(X)[2], "\n")

set.seed(SEED)
cvfit <- cv.glmnet(
  x = X, y = y,
  family       = FAMILY,
  alpha        = ALPHA,
  type.measure = "deviance",
  nlambda      = N_LAMBDA
)

## CV Curve & Coef Paths
cv_curve <- cv_lambda_summary(cvfit)
fwrite(
  cv_curve,
  file.path(OUT_MAIN,
            "lasso_drug_cv_curve_withIndication_noRegion_noTTO_noSerious.csv")
)
coef_paths_export(
  cvfit,
  FEATS,
  out_prefix = "lasso_drug_cvfit_withIndication_noRegion_noTTO_noSerious"
)

## Save model lock file
model_lock <- list(
  cvfit             = cvfit,
  seed              = SEED,
  family            = FAMILY,
  alpha             = ALPHA,
  nlambda           = N_LAMBDA,
  features          = FEATS,
  target            = TARGET_COL,
  id_col            = ID_COL,
  drug_class_dict   = drug_class_dict,
  indication_levels = levels(raw_lasso_dt$indication_grp)
)
saveRDS(
  model_lock,
  file = file.path(
    OUT_MAIN,
    "LASSO_Drug_Model_withIndication_noRegion_noTTO_noSerious.rds"
  )
)
cat("Model lock file saved to:\n  ",
    file.path(OUT_MAIN,
              "LASSO_Drug_Model_withIndication_noRegion_noTTO_noSerious.rds"),
    "\n")

## OOF Predictions and Performance
cat(">>> Calculating OOF predictions and performance metrics...\n")
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
            "cv_oof_predictions_drug_withIndication_noRegion_noTTO_noSerious.csv")
)
fwrite(
  data.table(
    fold     = 1:K_FOLDS,
    AUC      = oof$fold_auc,
    mean_AUC = oof$mean_auc
  ),
  file.path(OUT_MAIN,
            "OOF_auc_summary_drug_withIndication_noRegion_noTTO_noSerious.csv")
)

fwrite(
  data.table(
    Brier = brier_score(oof$oof$y, oof$oof$p_hat)
  ),
  file.path(OUT_MAIN,
            "brier_score_drug_withIndication_noRegion_noTTO_noSerious.csv")
)
fwrite(
  hl_test_result(oof$oof$y, oof$oof$p_hat, g = CAL_BINS),
  file.path(OUT_MAIN,
            "hl_test_drug_withIndication_noRegion_noTTO_noSerious.csv")
)
fwrite(
  calibration_table(oof$oof$y, oof$oof$p_hat, bins = CAL_BINS),
  file.path(OUT_MAIN,
            "calibration_table_oof_drug_withIndication_noRegion_noTTO_noSerious.csv")
)

pr <- PRROC::pr.curve(
  scores.class0 = oof$oof$p_hat[oof$oof$y == 1],
  scores.class1 = oof$oof$p_hat[oof$oof$y == 0]
)
fwrite(
  data.table(AUPRC = pr$auc.integral),
  file.path(OUT_MAIN,
            "auprc_drug_withIndication_noRegion_noTTO_noSerious.csv")
)

fwrite(
  decision_curve_table(oof$oof$y, oof$oof$p_hat),
  file.path(OUT_MAIN,
            "dca_table_drug_withIndication_noRegion_noTTO_noSerious.csv")
)

## SHAP-like importance (based on lambda.1se)
cat(">>> Calculating SHAP-like feature importance (lambda.1se)...\n")
b_1se <- as.matrix(coef(cvfit, s = "lambda.1se"))
b_1se <- b_1se[rownames(b_1se) != "(Intercept)", , drop = FALSE]
coef_vec <- as.numeric(b_1se)
names(coef_vec) <- rownames(b_1se)

idx   <- match(FEATS, names(coef_vec))
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
            "lasso_shap_like_importance_drug_withIndication_noRegion_noTTO_noSerious.csv")
)

## Refit logistic (based on lambda.1se non-zero features)
cat(">>> Refitting unpenalized logistic regression with lambda.1se non-zero features...\n")
sel <- names(coef_vec)[abs(coef_vec) > 0]
sel <- sel[!is.na(sel) & nzchar(sel)]
cat("Number of non-zero features at lambda.1se:", length(sel), "\n")
if (length(sel) == 0) {
  stop("No non-zero features found at lambda.1se.")
}

df_refit <- as.data.frame(as.matrix(X[, sel, drop = FALSE]))
colnames(df_refit) <- make.names(colnames(df_refit), unique = TRUE)
df_refit[[TARGET_COL]] <- y

form <- as.formula(paste(TARGET_COL, "~ ."))
cat("Refit formula:\n"); print(form)

fit_std  <- glm(form, data = df_refit, family = binomial(link = "logit"))
tidy_std <- broom::tidy(fit_std, conf.int = TRUE, exponentiate = TRUE)

std_out <- as.data.table(tidy_std)[term != "(Intercept)"]
setnames(
  std_out,
  c("estimate","conf.low","conf.high","p.value"),
  c("OR","CI_low","CI_high","p_value")
)

fwrite(
  std_out,
  file.path(OUT_MAIN,
            "lasso_refit_logistic_results_drug_withIndication_noRegion_noTTO_noSerious.csv")
)

## Sandwich Robust Estimation
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
            "lasso_refit_logistic_robust_drug_withIndication_noRegion_noTTO_noSerious.csv")
)

cat(">>> Drug-level + indication LASSO workflow complete.\n",
    "Outputs located in: ", OUT_MAIN, "\n")

## ----------------- Step 5d. SHAP-like: lambda.1se & lambda.min -----------------

cat(">>> Calculating SHAP-like feature importance (lambda.1se)...\n")

# ---- 1) lambda.1se version ----
b_1se <- as.matrix(coef(cvfit, s = "lambda.1se"))
b_1se <- b_1se[rownames(b_1se) != "(Intercept)", , drop = FALSE]

coef_vec_1se <- as.numeric(b_1se)
names(coef_vec_1se) <- rownames(b_1se)

idx_1se  <- match(FEATS, names(coef_vec_1se))
keep_1se <- which(!is.na(idx_1se))

imp_1se <- shap_like_global_importance(
  X[, keep_1se, drop = FALSE],
  FEATS[keep_1se],
  coef_vec_1se[idx_1se[keep_1se]]
)
imp_1se[, beta := coef_vec_1se[idx_1se[keep_1se]]]
imp_1se[, OR   := exp(beta)]

fwrite(
  imp_1se,
  file.path(
    OUT_MAIN,
    "lasso_shap_like_importance_drug_withIndication_noRegion_noTTO_noSerious_lambda1se.csv"
  )
)

# ---- 2) lambda.min version (including per-SD SHAP) ----
cat(">>> Calculating SHAP-like feature importance (lambda.min + per-SD)...\n")

b_min <- as.matrix(coef(cvfit, s = "lambda.min"))
b_min <- b_min[rownames(b_min) != "(Intercept)", , drop = FALSE]

coef_vec_min <- as.numeric(b_min)
names(coef_vec_min) <- rownames(b_min)

idx_min  <- match(FEATS, names(coef_vec_min))
keep_min <- which(!is.na(idx_min))

X_min    <- X[, keep_min, drop = FALSE]
feat_min <- FEATS[keep_min]
beta_min <- coef_vec_min[idx_min[keep_min]]

## Original scale SHAP-like (mean |x*beta|)
imp_min <- shap_like_global_importance(
  X_min,
  feat_min,
  beta_min
)
imp_min[, beta := beta_min]
imp_min[, OR   := exp(beta)]

fwrite(
  imp_min,
  file.path(
    OUT_MAIN,
    "lasso_shap_like_importance_drug_withIndication_noRegion_noTTO_noSerious_lambdaMIN.csv"
  )
)

## per-SD SHAP-like: |beta| * SD(x)
sd_x_min <- apply(X_min, 2, function(z) sd(as.numeric(z)))
sd_x_min[!is.finite(sd_x_min) | sd_x_min == 0] <- NA_real_

imp_min_sd <- data.table(
  feature     = feat_min,
  beta        = beta_min,
  OR          = exp(beta_min),
  sd_x        = sd_x_min,
  shap_per_sd = abs(beta_min) * sd_x_min
)

fwrite(
  imp_min_sd,
  file.path(
    OUT_MAIN,
    "lasso_shap_like_importance_perSD_drug_withIndication_noRegion_noTTO_noSerious_lambdaMIN.csv"
  )
)

cat("SHAP-like output files:\n")
cat("  lambda.1se original: ", 
    file.path(OUT_MAIN, "lasso_shap_like_importance_drug_withIndication_noRegion_noTTO_noSerious_lambda1se.csv"), "\n")
cat("  lambda.min original: ", 
    file.path(OUT_MAIN, "lasso_shap_like_importance_drug_withIndication_noRegion_noTTO_noSerious_lambdaMIN.csv"), "\n")
cat("  lambda.min per-SD: ", 
    file.path(OUT_MAIN, "lasso_shap_like_importance_perSD_drug_withIndication_noRegion_noTTO_noSerious_lambdaMIN.csv"), "\n")


## ----------------- Step 5e. Refit logistic (OR Table) ---------------------------

cat(">>> Refitting unpenalized logistic regression with lambda.1se features to export OR/CI...\n")

sel_1se <- names(coef_vec_1se)[abs(coef_vec_1se) > 0]
sel_1se <- sel_1se[!is.na(sel_1se) & nzchar(sel_1se)]

cat("Number of non-zero features at lambda.1se:", length(sel_1se), "\n")
if (length(sel_1se) == 0) {
  stop("No non-zero features at lambda.1se. Cannot perform Refit.")
}

## Build refit dataframe
df_refit <- as.data.frame(as.matrix(X[, sel_1se, drop = FALSE]))
colnames(df_refit) <- make.names(colnames(df_refit), unique = TRUE)
df_refit[[TARGET_COL]] <- y

form_refit <- as.formula(paste(TARGET_COL, "~ ."))
cat("Refit formula:\n")
print(form_refit)

## Standard logistic fit
fit_std  <- glm(form_refit, data = df_refit, family = binomial(link = "logit"))
tidy_std <- broom::tidy(fit_std, conf.int = TRUE, exponentiate = TRUE)

std_out <- as.data.table(tidy_std)[term != "(Intercept)"]
setnames(
  std_out,
  c("estimate", "conf.low", "conf.high", "p.value"),
  c("OR",       "CI_low",   "CI_high",   "p_value")
)

fwrite(
  std_out,
  file.path(OUT_MAIN, "lasso_refit_logistic_results_drug_withIndication_noRegion_noTTO_noSerious.csv")
)

## Sandwich Robust OR
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
  file.path(OUT_MAIN, "lasso_refit_logistic_robust_drug_withIndication_noRegion_noTTO_noSerious.csv")
)

cat("Refit ORs exported to:\n")
cat("  Standard logit: ", 
    file.path(OUT_MAIN, "lasso_refit_logistic_results_drug_withIndication_noRegion_noTTO_noSerious.csv"), "\n")
cat("  Robust logit: ", 
    file.path(OUT_MAIN, "lasso_refit_logistic_robust_drug_withIndication_noRegion_noTTO_noSerious.csv"), "\n")

cat(">>> Drug-level LASSO + SHAP + Refit process complete.\n")