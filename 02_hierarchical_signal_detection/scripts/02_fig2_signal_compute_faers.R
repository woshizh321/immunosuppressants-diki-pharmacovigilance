#!/usr/bin/env Rscript
# =============================================================================
# 02_hierarchical_signals_FAERS.R
# Hierarchical signal detection (ROR, PRR, chi-square) at SOC / HLGT / HLT / PT
# - Exposure: immunosuppressants (PS or SS) vs all other drugs
# - Outcome: presence of MedDRA term(s) in a case at each level
# - Levels: SOC, HLGT, HLT, PT (overall)
# - Stats: ROR (95% CI), PRR, chi-square (Haldane–Anscombe 0.5 correction)
# - Engine: DuckDB for memory-efficient case-term aggregation
# =============================================================================

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(data.table)
  library(fs)
  library(stringr)
})

# ------------------------------- Configuration -------------------------------

FAERS_MASTER <- "./data/FAERS_MASTER_FILE_2004-2024.parquet"  # <- change to your path
OUT_DIR       <- "./outputs/signals"
TMP_DUCK_DIR  <- "./.duck_tmp"
if (!dir_exists(OUT_DIR)) dir_create(OUT_DIR, recurse = TRUE)
if (!dir_exists(TMP_DUCK_DIR)) dir_create(TMP_DUCK_DIR, recurse = TRUE)

# Consistent ordering for all downstream reporting/plots
drug_order <- c("Calcineurin inhibitors",
                "mTOR inhibitors",
                "Antiproliferative agents",
                "Corticosteroids",
                "Biologics")

# Immunosuppressant mapping (uppercase, generic-wide). Extend as needed.
immu_dict <- list(
  "Calcineurin inhibitors" = c("TACROLIMUS","FK506","PROGRAF","ADVAGRAF",
                               "CYCLOSPORINE","NEORAL","SANDIMMUNE"),
  "mTOR inhibitors"        = c("EVEROLIMUS","AFINITOR","CERTICAN",
                               "SIROLIMUS","RAPAMUNE",
                               "TEMSIROLIMUS","TORISEL"),
  "Antiproliferative agents" = c("MYCOPHENOLATE MOFETIL","CELLCEPT",
                                 "MYCOPHENOLIC ACID","MYFORTIC",
                                 "AZATHIOPRINE","IMURAN"),
  "Corticosteroids"        = c("PREDNISONE","DELTASONE",
                               "METHYLPREDNISOLONE","MEDROL","SOLU-MEDROL",
                               "DEXAMETHASONE","DECADRON","PREDNISOLONE"),
  "Biologics"              = c("RITUXIMAB","MABTHERA","RITUXAN",
                               "BASILIXIMAB","SIMULECT","BELATACEPT","NULOJIX")
)

# ----------------------------- DuckDB connection -----------------------------

con <- dbConnect(duckdb::duckdb(), dbdir=":memory:", read_only=TRUE)
DBI::dbExecute(con, "PRAGMA threads=4;")
DBI::dbExecute(con, "PRAGMA memory_limit='8GB';")
DBI::dbExecute(con, sprintf("PRAGMA temp_directory='%s';", gsub("\\\\","/", TMP_DUCK_DIR)))
DBI::dbExecute(con, "PRAGMA preserve_insertion_order=false;")

# ------------------------------ Helper macros --------------------------------

DBI::dbExecute(con, "
  CREATE OR REPLACE MACRO NORM(s) AS
    TRIM(REGEXP_REPLACE(
      REGEXP_REPLACE(UPPER(COALESCE(s,'')), '[^A-Z0-9 \\-\\+]', ' ', 'g'
    ), '\\s+', ' ', 'g'));
")

# MASTER minimal projection (adjust if your schema differs)
DBI::dbExecute(con, sprintf("
  CREATE OR REPLACE VIEW v_master AS
  SELECT
    CAST(primaryid AS VARCHAR) AS primaryid,
    CAST(caseid    AS VARCHAR) AS caseid,
    NORM(role_cod)            AS role_cod,
    NORM(drugname)            AS drugname,
    NORM(pt)                  AS pt,
    NORM(hlt_name)            AS hlt_name,
    NORM(hlgt_name)           AS hlgt_name,
    NORM(soc_name)            AS soc_name
  FROM parquet_scan('%s');
", gsub("\\\\","/", FAERS_MASTER)))

# ----------------------- Build immunosuppressant exposure --------------------

# PS/SS drug rows only
DBI::dbExecute(con, "
  CREATE OR REPLACE VIEW v_drug_psss AS
  SELECT DISTINCT caseid, drugname
  FROM v_master
  WHERE role_cod IN ('PS','SS') AND LENGTH(drugname) > 0;
")

# Map drugname → drug_class / drug_generic (exact or contains).
# We use CONTAINS (LIKE) to be robust to brand→generic mixtures.
DBI::dbExecute(con, "DROP TABLE IF EXISTS immu_map; CREATE TEMP TABLE immu_map(drug TEXT, drug_class TEXT);")
for (cls in names(immu_dict)) {
  vals <- sprintf("('%s','%s')", gsub("'", "''", toupper(immu_dict[[cls]])), cls)
  DBI::dbExecute(con, paste0("INSERT INTO immu_map VALUES ", paste(vals, collapse=","), ";"))
}

# Exposure label at case-level: per drug_class and per drug (generic string)
DBI::dbExecute(con, "
  CREATE OR REPLACE VIEW v_exposure_class AS
  SELECT DISTINCT d.caseid, m.drug_class
  FROM v_drug_psss d
  JOIN immu_map m ON d.drugname LIKE ('%' || m.drug || '%');
")
DBI::dbExecute(con, "
  CREATE OR REPLACE VIEW v_exposure_drug AS
  SELECT DISTINCT d.caseid, d.drugname AS drug_generic
  FROM v_drug_psss d
  JOIN immu_map m ON d.drugname LIKE ('%' || m.drug || '%');
")

# Universe of cases for denominator (all cases that appear in MASTER)
DBI::dbExecute(con, "CREATE OR REPLACE VIEW v_all_cases AS SELECT DISTINCT caseid FROM v_master;")

# ---------------------- Case-term presence at each level ---------------------

DBI::dbExecute(con, "
  CREATE OR REPLACE VIEW v_case_pt AS
  SELECT DISTINCT caseid, pt
  FROM v_master WHERE LENGTH(pt) > 0;
")
DBI::dbExecute(con, "
  CREATE OR REPLACE VIEW v_case_hlt AS
  SELECT DISTINCT caseid, hlt_name AS hlt
  FROM v_master WHERE LENGTH(hlt_name) > 0;
")
DBI::dbExecute(con, "
  CREATE OR REPLACE VIEW v_case_hlgt AS
  SELECT DISTINCT caseid, hlgt_name AS hlgt
  FROM v_master WHERE LENGTH(hlgt_name) > 0;
")
DBI::dbExecute(con, "
  CREATE OR REPLACE VIEW v_case_soc AS
  SELECT DISTINCT caseid, soc_name AS soc
  FROM v_master WHERE LENGTH(soc_name) > 0;
")

# --------------------------- 2×2 counts generator ----------------------------
# For each exposure (class or drug) and each term, compute:
#   a = exposed & event
#   b = non-exposed & event
#   c = exposed & no-event
#   d = non-exposed & no-event

mk_counts_sql <- function(level = c("SOC","HLGT","HLT","PT"),
                          mode  = c("CLASS","DRUG")) {
  level <- match.arg(level)
  mode  <- match.arg(mode)
  
  if (level == "SOC")  term_tbl <- "v_case_soc";  term_col <- "soc"
  if (level == "HLGT") term_tbl <- "v_case_hlgt"; term_col <- "hlgt"
  if (level == "HLT")  term_tbl <- "v_case_hlt";  term_col <- "hlt"
  if (level == "PT")   term_tbl <- "v_case_pt";   term_col <- "pt"
  
  if (mode == "CLASS") {
    expo_tbl <- "v_exposure_class"; expo_col <- "drug_class"; level_tag <- "Class"
  } else {
    expo_tbl <- "v_exposure_drug";  expo_col <- "drug_generic"; level_tag <- "Drug"
  }
  
  sprintf("
    WITH
    universe AS (SELECT c.caseid FROM v_all_cases c),
    expo AS (SELECT DISTINCT caseid, %s AS exposure FROM %s),
    term AS (SELECT DISTINCT caseid, %s AS term FROM %s),

    join1 AS (
      SELECT u.caseid,
             e.exposure,
             CASE WHEN t.caseid IS NOT NULL THEN 1 ELSE 0 END AS event
      FROM universe u
      LEFT JOIN expo e ON u.caseid = e.caseid
      LEFT JOIN term t ON u.caseid = t.caseid
    ),

    agg AS (
      SELECT exposure,
             SUM(CASE WHEN exposure IS NOT NULL AND event=1 THEN 1 ELSE 0 END) AS a,
             SUM(CASE WHEN exposure IS NULL  AND event=1 THEN 1 ELSE 0 END) AS b,
             SUM(CASE WHEN exposure IS NOT NULL AND event=0 THEN 1 ELSE 0 END) AS c,
             SUM(CASE WHEN exposure IS NULL  AND event=0 THEN 1 ELSE 0 END) AS d
      FROM join1
      GROUP BY exposure
    ),

    # split by term: we need per term-exposure; do it explicitly
    term_expo AS (
      SELECT e.exposure, t.term, j.event
      FROM universe u
      LEFT JOIN (SELECT DISTINCT caseid, %s AS term FROM %s) t ON u.caseid = t.caseid
      LEFT JOIN (SELECT DISTINCT caseid, %s AS exposure FROM %s) e ON u.caseid = e.caseid
    ),

    tab AS (
      SELECT
        '%s' AS exposure_level,
        exposure,
        term AS signal_term,
        SUM(CASE WHEN exposure IS NOT NULL AND term IS NOT NULL THEN 1 ELSE 0 END) AS a,
        SUM(CASE WHEN exposure IS NULL  AND term IS NOT NULL THEN 1 ELSE 0 END) AS b,
        SUM(CASE WHEN exposure IS NOT NULL AND term IS NULL  THEN 1 ELSE 0 END) AS c,
        SUM(CASE WHEN exposure IS NULL  AND term IS NULL  THEN 1 ELSE 0 END) AS d
      FROM term_expo
      GROUP BY exposure, term
    )

    SELECT * FROM tab;
  ",
          expo_col, expo_tbl, term_col, term_tbl,
          term_col, term_tbl, expo_col, expo_tbl, level_tag
  )
}

# -------------------------- Pull counts for all levels ------------------------

get_counts_dt <- function(level, mode) {
  sql <- mk_counts_sql(level = level, mode = mode)
  as.data.table(DBI::dbGetQuery(con, sql))
}

counts_class_soc  <- get_counts_dt("SOC",  "CLASS")
counts_class_hlgt <- get_counts_dt("HLGT", "CLASS")
counts_class_hlt  <- get_counts_dt("HLT",  "CLASS")
counts_class_pt   <- get_counts_dt("PT",   "CLASS")

counts_drug_soc   <- get_counts_dt("SOC",  "DRUG")
counts_drug_hlgt  <- get_counts_dt("HLGT", "DRUG")
counts_drug_hlt   <- get_counts_dt("HLT",  "DRUG")
counts_drug_pt    <- get_counts_dt("PT",   "DRUG")

# ----------------------------- Metrics (R side) -------------------------------

calc_metrics <- function(dt) {
  # Haldane–Anscombe correction
  safe <- function(x) ifelse(is.na(x), 0, x)
  a <- safe(dt$a); b <- safe(dt$b); c <- safe(dt$c); d <- safe(dt$d)
  a <- ifelse(a==0, a+0.5, a); b <- ifelse(b==0, b+0.5, b)
  c <- ifelse(c==0, c+0.5, c); d <- ifelse(d==0, d+0.5, d)
  N <- a+b+c+d
  
  ror <- (a/c)/(b/d)
  se  <- sqrt(1/a + 1/b + 1/c + 1/d)
  lcl <- exp(log(ror) - 1.96*se)
  ucl <- exp(log(ror) + 1.96*se)
  
  prr <- (a/(a+c)) / (b/(b+d))
  chi2 <- (N * (a*d - b*c)^2) / ((a+b)*(c+d)*(a+c)*(b+d))
  
  out <- data.table(N=N, ROR=ror, LCL_ROR=lcl, UCL_ROR=ucl, PRR=prr, chi2=chi2)
  out[]
}

summarise_counts <- function(dt) {
  if (!nrow(dt)) return(dt)
  res <- cbind(
    dt[, .(exposure_level, exposure, signal_term, a,b,c,d)],
    calc_metrics(dt)
  )
  setcolorder(res, c("exposure_level","exposure","signal_term","N","a","b","c","d",
                     "ROR","LCL_ROR","UCL_ROR","PRR","chi2"))
  res[]
}

# Class-level results (respect fixed exposure order)
apply_order <- function(dt) {
  if (!nrow(dt)) return(dt)
  dt[, exposure := as.character(exposure)]
  dt[exposure %in% drug_order, exposure := factor(exposure, levels = drug_order)]
  setorder(dt, exposure, -ROR, -chi2, signal_term)
  dt[]
}

res_class_soc  <- apply_order(summarise_counts(counts_class_soc))
res_class_hlgt <- apply_order(summarise_counts(counts_class_hlgt))
res_class_hlt  <- apply_order(summarise_counts(counts_class_hlt))
res_class_pt   <- apply_order(summarise_counts(counts_class_pt))

# Drug-level results
res_drug_soc   <- summarise_counts(counts_drug_soc)
res_drug_hlgt  <- summarise_counts(counts_drug_hlgt)
res_drug_hlt   <- summarise_counts(counts_drug_hlt)
res_drug_pt    <- summarise_counts(counts_drug_pt)

# ------------------------------- Write outputs --------------------------------

fwrite(res_class_soc,  file.path(OUT_DIR, "SOC_signals_class_ROR_PRR_chi2.csv"))
fwrite(res_class_hlgt, file.path(OUT_DIR, "HLGT_signals_class_ROR_PRR_chi2.csv"))
fwrite(res_class_hlt,  file.path(OUT_DIR, "HLT_signals_class_ROR_PRR_chi2.csv"))
fwrite(res_class_pt,   file.path(OUT_DIR, "PT_signals_class_ROR_PRR_chi2.csv"))

fwrite(res_drug_soc,   file.path(OUT_DIR, "SOC_signals_drug_ROR_PRR_chi2.csv"))
fwrite(res_drug_hlgt,  file.path(OUT_DIR, "HLGT_signals_drug_ROR_PRR_chi2.csv"))
fwrite(res_drug_hlt,   file.path(OUT_DIR, "HLT_signals_drug_ROR_PRR_chi2.csv"))
fwrite(res_drug_pt,    file.path(OUT_DIR, "PT_signals_drug_ROR_PRR_chi2.csv"))

# Combined index (useful for quick browsing)
add_level <- function(dt, L) { if (!nrow(dt)) return(dt); dt[, level := L][] }
index_all <- rbindlist(list(
  add_level(res_class_soc,  "SOC-Class"),
  add_level(res_class_hlgt, "HLGT-Class"),
  add_level(res_class_hlt,  "HLT-Class"),
  add_level(res_class_pt,   "PT-Class"),
  add_level(res_drug_soc,   "SOC-Drug"),
  add_level(res_drug_hlgt,  "HLGT-Drug"),
  add_level(res_drug_hlt,   "HLT-Drug"),
  add_level(res_drug_pt,    "PT-Drug")
), use.names = TRUE, fill = TRUE)

setcolorder(index_all, c("level","exposure_level","exposure","signal_term","N","a","b","c","d",
                         "ROR","LCL_ROR","UCL_ROR","PRR","chi2"))
fwrite(index_all, file.path(OUT_DIR, "ALL_levels_signals_index.csv"))

message("✅ Hierarchical signal detection completed. Outputs => ", normalizePath(OUT_DIR, winslash = "/"))
try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)