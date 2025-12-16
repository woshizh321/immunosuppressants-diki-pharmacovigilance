# ============================================================
# Figure 1. Overall characteristics of immunosuppressant-
# associated reports in FAERS
#
# Purpose:
#   To summarize baseline reporting characteristics including
#   reporting volume, temporal trends, and demographic
#   distributions across major immunosuppressant classes.
#
# Data source:
#   FAERS MASTER FILE (local parquet, not redistributed)
#
# Key outputs:
#   - Annual report counts by drug category
#   - Sex distribution table
#   - Age median and IQR table
#   - Top-reporting countries
#
# Notes:
#   - QC-first workflow
#   - No region or TTO filtering applied
# ============================================================
# ===== A0. environment =====
suppressPackageStartupMessages({
  library(DBI); library(duckdb); library(data.table); library(stringr)
})

# ===== A1. connect DuckDB =====
con <- dbConnect(duckdb::duckdb(), dbdir=":memory:", read_only=TRUE)
# Do not use threads=auto (DuckDB requires INT); use a conservative integer
DBI::dbExecute(con, "PRAGMA threads=4;")
# Moderate memory limit (adjust based on machine availability)
DBI::dbExecute(con, "PRAGMA memory_limit='8GB';")

# ===== A2. Helper: Normalize Windows paths to forward slashes =====
norm_path <- function(p) gsub("\\\\", "/", p)

# ===== A3. Extract column names (get cols) =====
get_cols <- function(con, parquet_path){
  sql <- sprintf("SELECT * FROM read_parquet('%s') LIMIT 0", norm_path(parquet_path))
  df0 <- dbGetQuery(con, sql)   # Returns 0 rows but with column headers
  names(df0)
}

# ===== A4. Select column (case-insensitive, returns first match) =====
pick_col <- function(cols, candidates){
  # cols: character vector (actual names); candidates: potential matches
  ci <- tolower(cols)
  for (cand in candidates){
    hit <- which(ci == tolower(cand))
    if (length(hit)) return(cols[hit[1]])
  }
  return(NA_character_)
}

# ===== A5. Uniformly convert AGE to years =====
age_to_years_sql <- function(age_col, age_cod_col){
  # Generates DuckDB SQL snippet (CASE WHEN) for robust unit conversion
  sprintf("
    CASE
      WHEN TRY_CAST(NULLIF(CAST(%1$s AS VARCHAR), '') AS DOUBLE) IS NULL THEN NULL
      WHEN UPPER(COALESCE(%2$s,'')) IN ('','YR','YRS','YEAR','YEARS') THEN TRY_CAST(NULLIF(CAST(%1$s AS VARCHAR), '') AS DOUBLE)
      WHEN UPPER(%2$s) IN ('MON','MONTH','MONTHS') THEN TRY_CAST(NULLIF(CAST(%1$s AS VARCHAR), '') AS DOUBLE)/12.0
      WHEN UPPER(%2$s) IN ('WK','WEEK','WEEKS')   THEN TRY_CAST(NULLIF(CAST(%1$s AS VARCHAR), '') AS DOUBLE)/52.0
      WHEN UPPER(%2$s) IN ('DY','DAY','DAYS')     THEN TRY_CAST(NULLIF(CAST(%1$s AS VARCHAR), '') AS DOUBLE)/365.25
      WHEN UPPER(%2$s) IN ('DEC','DECADE')        THEN TRY_CAST(NULLIF(CAST(%1$s AS VARCHAR), '') AS DOUBLE)*10.0
      ELSE NULL
    END
  ", age_col, age_cod_col)
}

# ===== A6. Convert text to uppercase (for LIKE matching) =====
upper_coalesce <- function(col) sprintf("UPPER(COALESCE(%s,''))", col)

# ===== A7. Critical: PS+SS Exposure source priority (Strict Criteria) =====
# Priority 1: Use pre-aggregated long-text fields containing only PS+SS (e.g., ALL_DRUGS_PS_SS).
# Priority 2: If unavailable, filter ROLE %in% ('PS','SS') from drug-level table and aggregate to case-level.
# If both missing: Raise QC error (Avoid using ALL_DRUGS as it includes C/I/OT which contaminates the PS+SS criteria).
build_exposure_sql <- function(cols, dict_vec, 
                               case_key, # 'primaryid' or 'ID'
                               all_drugs_psss_col = NA, # e.g., 'ALL_DRUGS_PS_SS'
                               drugname_col = NA, role_col = NA){
  # dict_vec: category keywords (already in uppercase)
  like_or <- paste(sprintf("%s LIKE '%%%s%%'", "drug_blob", gsub("'", "''", toupper(dict_vec))), collapse=" OR ")
  
  if (!is.na(all_drugs_psss_col)) {
    # Perform LIKE matching directly on the concatenated field (report/case level)
    sql <- sprintf("
      SELECT %1$s AS case_key,
             CASE WHEN (%2$s) THEN 1 ELSE 0 END AS flag
      FROM tmp_base
    ",
                   case_key,
                   paste(sprintf("%s LIKE '%%%s%%'", upper_coalesce(all_drugs_psss_col), gsub("'", "''", toupper(dict_vec))), collapse=" OR ")
    )
    return(sql)
  }
  
  # Fallback to drug-level table: requires drugname_col + role_col
  if (is.na(drugname_col) || is.na(role_col)) {
    stop("PS+SS criteria requires: ALL_DRUGS_PS_SS (or equivalent) OR drug-level drugname + role columns. Neither available.")
  }
  
  # Construct subquery for drug text where role %in% (PS, SS), then aggregate by case_key
  # To avoid multiple scans, tmp_drugs_psss is built beforehand (see Section B)
  sql <- sprintf("
    SELECT %1$s AS case_key,
           CASE WHEN (%2$s) THEN 1 ELSE 0 END AS flag
    FROM (
      SELECT %1$s, STRING_AGG(%3$s, '|') AS drug_blob
      FROM tmp_drugs_psss
      GROUP BY %1$s
    )
  ",
                 case_key,
                 like_or,
                 upper_coalesce(drugname_col)
  )
  sql
}

FAERS_FILE <- "D:/FAERS/MASTER files/FAERS_MASTER_FILE_2004-2024_with_serious.parquet"
OUT_DIR    <- "D:/FAERS/Inhibitors"; if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, TRUE)

# Thread/Memory/Temp settings
DBI::dbExecute(con, "PRAGMA threads=2;")
DBI::dbExecute(con, "PRAGMA memory_limit='4GB';")
DBI::dbExecute(con, "PRAGMA temp_directory='D:/duck_tmp';")  # Ensure directory exists
DBI::dbExecute(con, "PRAGMA preserve_insertion_order=false;")


# ===== B1. Base Table (FAERS) =====

# Tool: Replace missing columns with NULL
sql_or_null_as <- function(colname, alias) {
  if (is.na(colname)) sprintf("NULL AS %s", alias) else colname
}

sql_base <- sprintf("
  CREATE OR REPLACE VIEW tmp_base AS
  SELECT
    %s AS case_key,                       -- Primary Key
    %s AS event_dt,                       -- Event Date
    %s,                                   -- OCCR_COUNTRY or NULL
    %s,                                   -- REPORTER_COUNTRY or NULL
    %s AS sex_raw,                        -- Gender
    %s AS age_raw,                        -- Age
    %s                                    -- Age Unit/Code
  FROM read_parquet('%s')
",
                    col_id,
                    col_eventdt,
                    sql_or_null_as(col_occr,   "occr_country"),
                    sql_or_null_as(col_report, "reporter_country"),
                    col_sex,
                    col_age,
                    sql_or_null_as(col_agecod, "age_cod"),
                    norm_path(FAERS_FILE)
)

DBI::dbExecute(con, sql_base)

# QC: Verify tmp_base creation
DBI::dbGetQuery(con, "PRAGMA table_info('tmp_base')")
DBI::dbGetQuery(con, "SELECT COUNT(*) FROM tmp_base")

# ===== B2. Build PS+SS Drug View & Year/Country Representative View =====

# 1) Filter PS/SS from drug-level records (Select essential columns only; normalize case)
stopifnot(!is.na(col_drugnm), !is.na(col_role))
DBI::dbExecute(con, sprintf("
  CREATE OR REPLACE VIEW tmp_drugs_psss AS
  SELECT
    %1$s AS case_key,
    %2$s AS drugname,
    %3$s AS role_cod
  FROM read_parquet('%4$s')
  WHERE UPPER(COALESCE(%3$s,'')) IN ('PS','SS')
", col_id, col_drugnm, col_role, norm_path(FAERS_FILE)))

# 2) Derive year and country_use from tmp_base, retaining original sex/age columns
DBI::dbExecute(con, "
  CREATE OR REPLACE VIEW tmp_with_year_country AS
  SELECT
    case_key,
    CASE
      WHEN LENGTH(CAST(event_dt AS VARCHAR))=8
           AND regexp_matches(CAST(event_dt AS VARCHAR),'^[0-9]{8}$')
      THEN TRY_CAST(SUBSTR(CAST(event_dt AS VARCHAR),1,4) AS INTEGER)
      ELSE NULL
    END AS year,
    CASE
      WHEN occr_country IS NOT NULL AND occr_country <> '' THEN occr_country
      ELSE reporter_country
    END AS country_use,
    sex_raw,
    age_raw,
    age_cod
  FROM tmp_base
")

# ---- QC (Short) ----
DBI::dbGetQuery(con, "SELECT COUNT(*) AS n_psss FROM tmp_drugs_psss")
DBI::dbGetQuery(con, "SELECT COUNT(*) AS n_cases FROM tmp_with_year_country")
DBI::dbGetQuery(con, "SELECT * FROM tmp_with_year_country WHERE year IS NOT NULL LIMIT 3")

# Define Dictionary Vectors
CNI  <- c("TACROLIMUS","FK506","PROGRAF","ADVAGRAF","CYCLOSPORINE","NEORAL","SANDIMMUNE")
MTOR <- c("EVEROLIMUS","AFINITOR","CERTICAN","SIROLIMUS","RAPAMUNE","TEMSIROLIMUS","TORISEL")
ANTI <- c("MYCOPHENOLATE MOFETIL","CELLCEPT","MYCOPHENOLIC ACID","MYFORTIC","AZATHIOPRINE","IMURAN")
CS   <- c("PREDNISONE","DELTASONE","METHYLPREDNISOLONE","MEDROL","SOLU-MEDROL","DEXAMETHASONE","DECADRON","PREDNISOLONE")
BIO  <- c("RITUXIMAB","MABTHERA","RITUXAN","BASILIXIMAB","SIMULECT","BELATACEPT","NULOJIX")

like_or <- function(vec) paste(sprintf("UPPER(COALESCE(d.drugname,'')) LIKE '%%%s%%'", gsub("'", "''", vec)), collapse=" OR ")

# Create Flag Tables for each drug category
DBI::dbExecute(con, sprintf("
  CREATE OR REPLACE TABLE flag_cni AS
  SELECT DISTINCT c.case_key, 1 AS flag
  FROM tmp_with_year_country c
  WHERE EXISTS (SELECT 1 FROM tmp_drugs_psss d WHERE d.case_key=c.case_key AND (%s))
", like_or(CNI)))

DBI::dbExecute(con, sprintf("
  CREATE OR REPLACE TABLE flag_mtor AS
  SELECT DISTINCT c.case_key, 1 AS flag
  FROM tmp_with_year_country c
  WHERE EXISTS (SELECT 1 FROM tmp_drugs_psss d WHERE d.case_key=c.case_key AND (%s))
", like_or(MTOR)))

DBI::dbExecute(con, sprintf("
  CREATE OR REPLACE TABLE flag_anti AS
  SELECT DISTINCT c.case_key, 1 AS flag
  FROM tmp_with_year_country c
  WHERE EXISTS (SELECT 1 FROM tmp_drugs_psss d WHERE d.case_key=c.case_key AND (%s))
", like_or(ANTI)))

DBI::dbExecute(con, sprintf("
  CREATE OR REPLACE TABLE flag_cs AS
  SELECT DISTINCT c.case_key, 1 AS flag
  FROM tmp_with_year_country c
  WHERE EXISTS (SELECT 1 FROM tmp_drugs_psss d WHERE d.case_key=c.case_key AND (%s))
", like_or(CS)))

DBI::dbExecute(con, sprintf("
  CREATE OR REPLACE TABLE flag_bio AS
  SELECT DISTINCT c.case_key, 1 AS flag
  FROM tmp_with_year_country c
  WHERE EXISTS (SELECT 1 FROM tmp_drugs_psss d WHERE d.case_key=c.case_key AND (%s))
", like_or(BIO)))

# QC (Count should be > 0 for each category)
DBI::dbGetQuery(con, "SELECT (SELECT COUNT(*) FROM flag_cni) cni,
                             (SELECT COUNT(*) FROM flag_mtor) mtor,
                             (SELECT COUNT(*) FROM flag_anti) anti,
                             (SELECT COUNT(*) FROM flag_cs)   cs,
                             (SELECT COUNT(*) FROM flag_bio)  bio")

# Summary QC Table
DBI::dbGetQuery(con, "
  SELECT 'cni'  AS cat, COUNT(*) AS n FROM flag_cni
  UNION ALL
  SELECT 'mtor', COUNT(*) FROM flag_mtor
  UNION ALL
  SELECT 'anti', COUNT(*) FROM flag_anti
  UNION ALL
  SELECT 'cs',   COUNT(*) FROM flag_cs
  UNION ALL
  SELECT 'bio',  COUNT(*) FROM flag_bio
  ORDER BY cat
")

# Merging flags back into the case table
DBI::dbExecute(con, "
  CREATE OR REPLACE TABLE faers_case_psss AS
  SELECT c.*, COALESCE(f.flag,0) AS flag_cni
  FROM tmp_with_year_country c LEFT JOIN flag_cni f USING(case_key);
")
DBI::dbExecute(con, "
  CREATE OR REPLACE TABLE faers_case_psss AS
  SELECT c.*, COALESCE(m.flag,0) AS flag_mtor
  FROM faers_case_psss c LEFT JOIN flag_mtor m USING(case_key);
")
DBI::dbExecute(con, "
  CREATE OR REPLACE TABLE faers_case_psss AS
  SELECT c.*, COALESCE(a.flag,0) AS flag_anti
  FROM faers_case_psss c LEFT JOIN flag_anti a USING(case_key);
")
DBI::dbExecute(con, "
  CREATE OR REPLACE TABLE faers_case_psss AS
  SELECT c.*, COALESCE(s.flag,0) AS flag_cs
  FROM faers_case_psss c LEFT JOIN flag_cs s USING(case_key);
")
DBI::dbExecute(con, "
  CREATE OR REPLACE TABLE faers_case_psss AS
  SELECT c.*, COALESCE(b.flag,0) AS flag_bio
  FROM faers_case_psss c LEFT JOIN flag_bio b USING(case_key);
")

# QC Final Counts
DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM faers_case_psss")
DBI::dbGetQuery(con, "SELECT SUM(flag_cni) cni, SUM(flag_mtor) mtor, SUM(flag_anti) anti, SUM(flag_cs) cs, SUM(flag_bio) bio FROM faers_case_psss")

# ===== B6. Age Distribution: Dual-Track Output =====

# SQL snippet: Normalize age_raw + age_cod to years
age_sql <- "
  CASE
    WHEN TRY_CAST(NULLIF(CAST(age_raw AS VARCHAR), '') AS DOUBLE) IS NULL THEN NULL
    WHEN UPPER(COALESCE(age_cod,'')) IN ('','YR','YRS','YEAR','YEARS') 
         THEN TRY_CAST(NULLIF(CAST(age_raw AS VARCHAR), '') AS DOUBLE)
    WHEN UPPER(age_cod) IN ('MON','MONTH','MONTHS') 
         THEN TRY_CAST(NULLIF(CAST(age_raw AS VARCHAR), '') AS DOUBLE)/12.0
    WHEN UPPER(age_cod) IN ('WK','WEEK','WEEKS')   
         THEN TRY_CAST(NULLIF(CAST(age_raw AS VARCHAR), '') AS DOUBLE)/52.0
    WHEN UPPER(age_cod) IN ('DY','DAY','DAYS')     
         THEN TRY_CAST(NULLIF(CAST(age_raw AS VARCHAR), '') AS DOUBLE)/365.25
    WHEN UPPER(age_cod) IN ('DEC','DECADE')        
         THEN TRY_CAST(NULLIF(CAST(age_raw AS VARCHAR), '') AS DOUBLE)*10.0
    ELSE NULL
  END
"

# --- Table 1: Median + IQR ---
age_median <- DBI::dbGetQuery(con, sprintf("
  WITH base AS (
    SELECT %s AS age_years, flag_cni, flag_mtor, flag_anti, flag_cs, flag_bio
    FROM faers_case_psss
  )
  SELECT 'Calcineurin inhibitors' AS drug_category,
         COUNT(age_years) AS N_nonmissing,
         quantile_cont(age_years,0.5)  AS median_age,
         quantile_cont(age_years,0.25) AS q1,
         quantile_cont(age_years,0.75) AS q3
  FROM base WHERE flag_cni=1
  UNION ALL
  SELECT 'mTOR inhibitors', COUNT(age_years), quantile_cont(age_years,0.5), quantile_cont(age_years,0.25), quantile_cont(age_years,0.75)
  FROM base WHERE flag_mtor=1
  UNION ALL
  SELECT 'Antiproliferative agents', COUNT(age_years), quantile_cont(age_years,0.5), quantile_cont(age_years,0.25), quantile_cont(age_years,0.75)
  FROM base WHERE flag_anti=1
  UNION ALL
  SELECT 'Corticosteroids', COUNT(age_years), quantile_cont(age_years,0.5), quantile_cont(age_years,0.25), quantile_cont(age_years,0.75)
  FROM base WHERE flag_cs=1
  UNION ALL
  SELECT 'Biologics', COUNT(age_years), quantile_cont(age_years,0.5), quantile_cont(age_years,0.25), quantile_cont(age_years,0.75)
  FROM base WHERE flag_bio=1
", age_sql))
data.table::fwrite(age_median, "D:/FAERS/Inhibitors/age_median_IQR_by_category_PS+SS.csv")

# --- Table 2: Mean ± SD ---
age_mean <- DBI::dbGetQuery(con, sprintf("
  WITH base AS (
    SELECT %s AS age_years, flag_cni, flag_mtor, flag_anti, flag_cs, flag_bio
    FROM faers_case_psss
  )
  SELECT 'Calcineurin inhibitors' AS drug_category,
         COUNT(age_years) AS N_nonmissing,
         avg(age_years)   AS mean_age,
         stddev_samp(age_years) AS sd_age
  FROM base WHERE flag_cni=1
  UNION ALL
  SELECT 'mTOR inhibitors', COUNT(age_years), avg(age_years), stddev_samp(age_years)
  FROM base WHERE flag_mtor=1
  UNION ALL
  SELECT 'Antiproliferative agents', COUNT(age_years), avg(age_years), stddev_samp(age_years)
  FROM base WHERE flag_anti=1
  UNION ALL
  SELECT 'Corticosteroids', COUNT(age_years), avg(age_years), stddev_samp(age_years)
  FROM base WHERE flag_cs=1
  UNION ALL
  SELECT 'Biologics', COUNT(age_years), avg(age_years), stddev_samp(age_years)
  FROM base WHERE flag_bio=1
", age_sql))
data.table::fwrite(age_mean, "D:/FAERS/Inhibitors/age_mean_SD_by_category_PS+SS.csv")

# ===== B7. Age Stats Export (Combined) =====
age_expr <- age_to_years_sql("age_raw","age_cod")

age_stats_faers <- DBI::dbGetQuery(con, sprintf("
  WITH base AS (
    SELECT %s AS age_years, flag_cni, flag_mtor, flag_anti, flag_cs, flag_bio
    FROM faers_case_psss
  )
  SELECT 'Calcineurin inhibitors' AS drug_category, COUNT(age_years) AS N_nonmissing,
         quantile_cont(age_years,0.5)  AS median_age,
         quantile_cont(age_years,0.25) AS q1,
         quantile_cont(age_years,0.75) AS q3
  FROM base WHERE flag_cni=1
  UNION ALL
  SELECT 'mTOR inhibitors', COUNT(age_years),
         quantile_cont(age_years,0.5), quantile_cont(age_years,0.25), quantile_cont(age_years,0.75)
  FROM base WHERE flag_mtor=1
  UNION ALL
  SELECT 'Antiproliferative agents', COUNT(age_years),
         quantile_cont(age_years,0.5), quantile_cont(age_years,0.25), quantile_cont(age_years,0.75)
  FROM base WHERE flag_anti=1
  UNION ALL
  SELECT 'Corticosteroids', COUNT(age_years),
         quantile_cont(age_years,0.5), quantile_cont(age_years,0.25), quantile_cont(age_years,0.75)
  FROM base WHERE flag_cs=1
  UNION ALL
  SELECT 'Biologics', COUNT(age_years),
         quantile_cont(age_years,0.5), quantile_cont(age_years,0.25), quantile_cont(age_years,0.75)
  FROM base WHERE flag_bio=1
", age_expr))

data.table::fwrite(age_stats_faers, "D:/FAERS/Inhibitors/age_stats_by_category_PS+SS.csv")

# ===== B8. Sex Distribution =====
sex_dist_faers <- DBI::dbGetQuery(con, "
  WITH sex_std AS (
    SELECT CASE
             WHEN UPPER(COALESCE(sex_raw,'')) IN ('F','FEMALE','2') THEN 'F'
             WHEN UPPER(COALESCE(sex_raw,'')) IN ('M','MALE','1')   THEN 'M'
             ELSE 'Other' END AS sex3,
           flag_cni, flag_mtor, flag_anti, flag_cs, flag_bio
    FROM faers_case_psss
  )
  SELECT 'Calcineurin inhibitors' AS drug_category, sex3, COUNT(*) n FROM sex_std WHERE flag_cni=1 GROUP BY sex3
  UNION ALL SELECT 'mTOR inhibitors', sex3, COUNT(*) FROM sex_std WHERE flag_mtor=1 GROUP BY sex3
  UNION ALL SELECT 'Antiproliferative agents', sex3, COUNT(*) FROM sex_std WHERE flag_anti=1 GROUP BY sex3
  UNION ALL SELECT 'Corticosteroids', sex3, COUNT(*) FROM sex_std WHERE flag_cs=1 GROUP BY sex3
  UNION ALL SELECT 'Biologics', sex3, COUNT(*) FROM sex_std WHERE flag_bio=1 GROUP BY sex3
  ORDER BY drug_category, sex3
")
data.table::fwrite(sex_dist_faers, "D:/FAERS/Inhibitors/sex_distribution_by_category_PS+SS.csv")

# ===== B9. Country Distribution =====
country_dist_faers <- DBI::dbGetQuery(con, "
  SELECT 'Calcineurin inhibitors' AS drug_category, country_use, COUNT(*) case_count
    FROM faers_case_psss WHERE flag_cni=1 AND country_use IS NOT NULL AND country_use<>'' GROUP BY country_use
  UNION ALL SELECT 'mTOR inhibitors', country_use, COUNT(*) FROM faers_case_psss WHERE flag_mtor=1 AND country_use IS NOT NULL AND country_use<>'' GROUP BY country_use
  UNION ALL SELECT 'Antiproliferative agents', country_use, COUNT(*) FROM faers_case_psss WHERE flag_anti=1 AND country_use IS NOT NULL AND country_use<>'' GROUP BY country_use
  UNION ALL SELECT 'Corticosteroids', country_use, COUNT(*) FROM faers_case_psss WHERE flag_cs=1 AND country_use IS NOT NULL AND country_use<>'' GROUP BY country_use
  UNION ALL SELECT 'Biologics', country_use, COUNT(*) FROM faers_case_psss WHERE flag_bio=1 AND country_use IS NOT NULL AND country_use<>'' GROUP BY country_use
  ORDER BY drug_category, case_count DESC
")
data.table::fwrite(country_dist_faers, "D:/FAERS/Inhibitors/country_distribution_by_category_PS+SS.csv")