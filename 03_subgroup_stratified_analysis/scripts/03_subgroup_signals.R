#!/usr/bin/env Rscript
# =============================================================================
# 03_subgroup_signals_DIKI.R
# Subgroup disproportionality for DIKI (ROR, PRR, Chi2, IC)
# Subgroups: Sex, Age, Region(5), Indication(3), Seriousness(2),
#            and Indication × Seriousness.
# =============================================================================

suppressPackageStartupMessages({
  library(DBI); library(duckdb); library(data.table); library(glue); library(stringr); library(fs)
})

# ------------------------------ Configuration --------------------------------
# Master parquet (same as in prior steps; contains primaryid/drugname/pt_code/age/region/indi/serious*)
MASTER_PARQUET <- "D:/FAERS/MASTER files/FAERS_MASTER_FILE_2004-2024_with_serious.parquet"

OUT_DIR <- "D:/FAERS/Inhibitors/results"
TMP_DUCK <- "D:/duck_tmp"
dir_create(OUT_DIR, recurse = TRUE); dir_create(TMP_DUCK, recurse = TRUE)

# Fixed display order for drug classes (used across all outputs)
drug_order <- c("Calcineurin inhibitors",
                "mTOR inhibitors",
                "Antiproliferative agents",
                "Corticosteroids",
                "Biologics")

# Immunosuppressant dictionary (exact uppercase match against DRUGNAME)
drug_classes <- list(
  "Calcineurin inhibitors"   = c("TACROLIMUS","FK506","PROGRAF","ADVAGRAF","CYCLOSPORINE","NEORAL","SANDIMMUNE"),
  "mTOR inhibitors"          = c("EVEROLIMUS","AFINITOR","CERTICAN","SIROLIMUS","RAPAMUNE","TEMSIROLIMUS","TORISEL"),
  "Antiproliferative agents" = c("MYCOPHENOLATE MOFETIL","CELLCEPT","MYCOPHENOLIC ACID","MYFORTIC","AZATHIOPRINE","IMURAN"),
  "Corticosteroids"          = c("PREDNISONE","DELTASONE","METHYLPREDNISOLONE","MEDROL","SOLU-MEDROL","DEXAMETHASONE","DECADRON","PREDNISOLONE"),
  "Biologics"                = c("RITUXIMAB","MABTHERA","RITUXAN","BASILIXIMAB","SIMULECT","BELATACEPT","NULOJIX")
)

# 16 DIKI PTs (authoritative list you provided)
diki_pts_csv <- "
pt_code,pt_name
10001580,Albuminuria
10002847,Anuria
10005483,Blood creatinine increased
10018358,Glomerular filtration rate decreased
10018867,Haematuria
10027525,Microalbuminuria
10029155,Nephropathy toxic
10030302,Oliguria
10037032,Proteinuria
10038428,Renal disorder
10038435,Renal failure
10061480,Renal function test abnormal
10062237,Renal impairment
10062747,Hypercreatininaemia
10064848,Chronic kidney disease
10069339,Acute kidney injury
"

# ------------------------------- DuckDB --------------------------------------
con <- dbConnect(duckdb::duckdb(), dbdir=":memory:", read_only = FALSE)
DBI::dbExecute(con, "PRAGMA threads=4;")
DBI::dbExecute(con, "PRAGMA memory_limit='12GB';")
DBI::dbExecute(con, glue("PRAGMA temp_directory='{gsub('\\\\','/', TMP_DUCK)}';"))

# Minimal projection (uppercasing text for stable matching)
DBI::dbExecute(con, glue("
  CREATE OR REPLACE TEMP VIEW faers_min AS
  SELECT
    CAST(primaryid AS BIGINT)                                      AS primaryid,
    UPPER(CAST(drugname AS VARCHAR))                               AS drugname_up,
    CAST(pt_code AS BIGINT)                                        AS pt_code,
    -- Age (numeric if available)
    TRY_CAST(age AS DOUBLE)                                        AS age_raw,
    -- Region (you already pre-computed region; fallback to upper text)
    UPPER(CAST(region AS VARCHAR))                                 AS region_up,
    -- Indication (prefer normalized; fallback to original)
    UPPER(COALESCE(CAST(indi_pt_norm AS VARCHAR), CAST(indi_pt AS VARCHAR))) AS indi_up,
    -- Sex (prefer single-letter if present; fallback to upper)
    UPPER(CAST(sex AS VARCHAR))                                    AS sex_up,

    -- Seriousness (row-level flags; default 0 if null)
    COALESCE(CAST(serious_flag               AS INTEGER), 0) AS serious_flag_row,
    COALESCE(CAST(serious_death              AS INTEGER), 0) AS serious_death_row,
    COALESCE(CAST(serious_life_threat        AS INTEGER), 0) AS serious_life_threat_row,
    COALESCE(CAST(serious_hospitalization    AS INTEGER), 0) AS serious_hospitalization_row,
    COALESCE(CAST(serious_disability         AS INTEGER), 0) AS serious_disability_row,
    COALESCE(CAST(serious_congenital_anomaly AS INTEGER), 0) AS serious_congenital_anomaly_row,
    COALESCE(CAST(serious_other              AS INTEGER), 0) AS serious_other_row
  FROM read_parquet('{gsub(\"'\",\"''\", MASTER_PARQUET)}');
"))

# -------------------------- Case-level scaffolding ---------------------------
# Sex (F/M/Other)
DBI::dbExecute(con, "
  CREATE OR REPLACE TEMP VIEW case_sex AS
  SELECT primaryid,
         CASE
           WHEN sex_up IN ('F','FEMALE','2') THEN 'F'
           WHEN sex_up IN ('M','MALE','1')   THEN 'M'
           ELSE 'Other'
         END AS sex_grp
  FROM (SELECT DISTINCT primaryid, sex_up FROM faers_min);
")

# Age → representative per case, then group
DBI::dbExecute(con, "
  CREATE OR REPLACE TEMP VIEW case_age AS
  SELECT primaryid, ANY_VALUE(age_raw) AS age
  FROM faers_min GROUP BY primaryid;
")
DBI::dbExecute(con, "
  CREATE OR REPLACE TEMP VIEW case_age_grp AS
  SELECT primaryid,
         CASE
           WHEN age IS NULL OR age <= 0 THEN NULL
           WHEN age < 18 THEN '<18'
           WHEN age BETWEEN 18 AND 39 THEN '18-39'
           WHEN age BETWEEN 40 AND 64 THEN '40-64'
           WHEN age >= 65 THEN '65+'
         END AS age_group
  FROM case_age;
")

# Region (US/Europe/Japan/China/Other); assumes region_up already harmonized
DBI::dbExecute(con, "
  CREATE OR REPLACE TEMP VIEW case_region AS
  SELECT primaryid,
         CASE
           WHEN region_up IS NULL OR region_up='' THEN 'OTHER'
           WHEN region_up IN ('US','EUROPE','JAPAN','CHINA') THEN region_up
           ELSE 'OTHER'
         END AS region_grp
  FROM (SELECT DISTINCT primaryid, region_up FROM faers_min);
")

# Seriousness (case-level): “any sub-flag = 1” defines serious
DBI::dbExecute(con, "
  CREATE OR REPLACE TEMP VIEW case_serious AS
  SELECT
    primaryid,
    MAX(
      (serious_death_row
      + serious_life_threat_row
      + serious_hospitalization_row
      + serious_disability_row
      + serious_congenital_anomaly_row
      + serious_other_row) > 0
    )::INT AS serious_calc_case
  FROM faers_min
  GROUP BY primaryid;
")
DBI::dbExecute(con, "
  CREATE OR REPLACE TEMP VIEW case_serious_grp AS
  SELECT primaryid,
         CASE WHEN serious_calc_case=1 THEN 'Serious' ELSE 'Non-serious' END AS serious_grp
  FROM case_serious;
")

# Indication groups (3-class; keywords reviewed by you)
kw_transplant <- c("TRANSPLANT","TRANSPLANTATION","GRAFT","RENAL TRANSPLANT","KIDNEY TRANSPLANT",
                   "HEART TRANSPLANT","LIVER TRANSPLANT","PANCREAS TRANSPLANT","BONE MARROW TRANSPLANT",
                   "HSCT","SCT")
kw_autoimmune <- c("SYSTEMIC LUPUS","LUPUS","SLE","RHEUMATOID ARTHRITIS","RA","PSORIASIS","PSORIATIC",
                   "ANKYLOSING SPONDYLITIS","CROHN","ULCERATIVE COLITIS","IBD","INFLAMMATORY BOWEL",
                   "VASCULITIS","ANCA","MULTIPLE SCLEROSIS","MYASTHENIA","SJOGREN","SCLERODERMA",
                   "DERMATOMYOSITIS","POLYMYOSITIS","AUTOIMMUNE","IGA NEPHROPATHY","LUPUS NEPHRITIS")
pat_tx <- paste0("(", paste(kw_transplant, collapse="|"), ")")
pat_ai <- paste0("(", paste(kw_autoimmune, collapse="|"), ")")

DBI::dbExecute(con, "
  CREATE OR REPLACE TEMP VIEW case_indi AS
  SELECT primaryid, ANY_VALUE(indi_up) AS indi_up
  FROM faers_min GROUP BY primaryid;
")
DBI::dbExecute(con, glue("
  CREATE OR REPLACE TEMP VIEW case_indi_grp AS
  SELECT primaryid,
         CASE
           WHEN indi_up IS NULL OR indi_up='' THEN 'Other'
           WHEN REGEXP_MATCHES(indi_up, '{gsub(\"'\",\"''\", pat_tx)}') THEN 'Transplantation'
           WHEN REGEXP_MATCHES(indi_up, '{gsub(\"'\",\"''\", pat_ai)}') THEN 'Autoimmune diseases'
           ELSE 'Other'
         END AS indication_grp
  FROM case_indi;
"))

# Exposure (any record hits dictionary → exposed=1). Exact match on DRUGNAME_UP.
map_df <- rbindlist(lapply(names(drug_classes), function(cat)
  data.table(drugname_up = drug_classes[[cat]], category = cat)
))
DBI::dbExecute(con, "CREATE OR REPLACE TEMP TABLE drug_map(drugname_up VARCHAR, category VARCHAR);")
DBI::dbAppendTable(con, "drug_map", as.data.frame(map_df))

DBI::dbExecute(con, "
  CREATE OR REPLACE TEMP VIEW case_expo AS
  SELECT
    f.primaryid,
    MAX(CASE WHEN m.category='Calcineurin inhibitors'   THEN 1 ELSE 0 END) AS exposed_CNI,
    MAX(CASE WHEN m.category='mTOR inhibitors'          THEN 1 ELSE 0 END) AS exposed_mTOR,
    MAX(CASE WHEN m.category='Antiproliferative agents' THEN 1 ELSE 0 END) AS exposed_ANTI,
    MAX(CASE WHEN m.category='Corticosteroids'          THEN 1 ELSE 0 END) AS exposed_CS,
    MAX(CASE WHEN m.category='Biologics'                THEN 1 ELSE 0 END) AS exposed_BIO
  FROM faers_min f
  JOIN drug_map m ON f.drugname_up = m.drugname_up
  GROUP BY f.primaryid;
")

# DIKI endpoint: any of PT16
diki_pts <- fread(text = diki_pts_csv)
DBI::dbExecute(con, "CREATE OR REPLACE TEMP TABLE diki_pt_list(pt_code BIGINT, pt_name VARCHAR);")
DBI::dbAppendTable(con, "diki_pt_list", as.data.frame(diki_pts))

DBI::dbExecute(con, "
  CREATE OR REPLACE TEMP VIEW case_diki AS
  SELECT primaryid, 1 AS diki
  FROM faers_min
  WHERE pt_code IN (SELECT pt_code FROM diki_pt_list)
  GROUP BY primaryid;
")

# Universe join (include unexposed and non-DIKI)
DBI::dbExecute(con, "
  CREATE OR REPLACE TEMP VIEW case_universe AS
  SELECT DISTINCT
    i.primaryid,
    COALESCE(s.sex_grp,'Other')            AS sex_grp,
    ageg.age_group                          AS age_group,
    COALESCE(r.region_grp,'OTHER')         AS region_grp,
    COALESCE(ig.indication_grp,'Other')    AS indication_grp,
    COALESCE(sg.serious_grp,'Non-serious') AS serious_grp,

    COALESCE(e.exposed_CNI,0)  AS exposed_CNI,
    COALESCE(e.exposed_mTOR,0) AS exposed_mTOR,
    COALESCE(e.exposed_ANTI,0) AS exposed_ANTI,
    COALESCE(e.exposed_CS,0)   AS exposed_CS,
    COALESCE(e.exposed_BIO,0)  AS exposed_BIO
  FROM (SELECT DISTINCT primaryid FROM faers_min) i
  LEFT JOIN case_sex        s   USING(primaryid)
  LEFT JOIN case_age_grp    ageg USING(primaryid)
  LEFT JOIN case_region     r   USING(primaryid)
  LEFT JOIN case_indi_grp   ig  USING(primaryid)
  LEFT JOIN case_serious_grp sg USING(primaryid)
  LEFT JOIN case_expo       e   USING(primaryid);
")

DBI::dbExecute(con, "
  CREATE OR REPLACE TEMP VIEW case_diki_flag AS
  SELECT u.*, CASE WHEN d.primaryid IS NULL THEN 0 ELSE 1 END AS diki
  FROM case_universe u
  LEFT JOIN case_diki d USING(primaryid);
")

# ---------------------------- R helpers (metrics) -----------------------------
calc_metrics <- function(a,b,c,d){
  aa <- as.numeric(a); bb <- as.numeric(b); cc <- as.numeric(c); dd <- as.numeric(d)
  # Haldane–Anscombe continuity correction
  if (min(aa,bb,cc,dd) == 0) { aa<-aa+0.5; bb<-bb+0.5; cc<-cc+0.5; dd<-dd+0.5 }
  N <- aa+bb+cc+dd
  
  ROR <- (aa*dd)/(bb*cc)
  se  <- sqrt(1/aa + 1/bb + 1/cc + 1/dd)
  ROR_LCL <- exp(log(ROR) - 1.96*se)
  ROR_UCL <- exp(log(ROR) + 1.96*se)
  
  PRR <- (aa/(aa+bb)) / (cc/(cc+dd))
  Chi2 <- as.numeric(((aa*dd - bb*cc)^2 * N) / ((aa+bb)*(cc+dd)*(aa+cc)*(bb+dd)))
  
  # BCPNN Information Component (normal approximation for CI)
  IC <- log2((aa*N)/((aa+bb)*(aa+cc)))
  var_ic <- (1/log(2)^2) * (1/aa - 1/(aa+bb) - 1/(aa+cc) + 1/N)
  se_ic <- sqrt(pmax(var_ic, 0))
  IC025 <- IC - 1.96*se_ic
  IC975 <- IC + 1.96*se_ic
  
  data.table(ROR, ROR_LCL, ROR_UCL, PRR, Chi2, IC, IC025, IC975)
}

# Utility to enforce drug order and write CSV
finalize_and_write <- function(dt, outfile, subgroup_col, subgroup_levels = NULL) {
  if (!nrow(dt)) { fwrite(dt, outfile); message("Saved (empty): ", outfile); return(invisible(dt)) }
  # Exposure order
  dt[, exposure := factor(exposure, levels = drug_order)]
  # Optional subgroup ordering
  if (!is.null(subgroup_levels) && subgroup_col %in% names(dt)) {
    dt[[subgroup_col]] <- factor(dt[[subgroup_col]], levels = subgroup_levels)
  }
  setorder(dt, exposure, !!as.name(subgroup_col))
  fwrite(dt, outfile, bom = TRUE)
  message("Saved -> ", outfile)
  invisible(dt)
}

# ----------------------------- 1) Sex (F/M/Other) ----------------------------
sql_sex <- "
WITH base AS (SELECT * FROM case_diki_flag)
SELECT 'Calcineurin inhibitors' AS exposure, sex_grp,
       SUM(CASE WHEN exposed_CNI  = 1 AND diki = 1 THEN 1 ELSE 0 END) AS a,
       SUM(CASE WHEN exposed_CNI  = 1 AND diki = 0 THEN 1 ELSE 0 END) AS b,
       SUM(CASE WHEN exposed_CNI  = 0 AND diki = 1 THEN 1 ELSE 0 END) AS c,
       SUM(CASE WHEN exposed_CNI  = 0 AND diki = 0 THEN 1 ELSE 0 END) AS d
FROM base GROUP BY sex_grp
UNION ALL
SELECT 'mTOR inhibitors', sex_grp,
       SUM(CASE WHEN exposed_mTOR = 1 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_mTOR = 1 AND diki = 0 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_mTOR = 0 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_mTOR = 0 AND diki = 0 THEN 1 ELSE 0 END)
FROM base GROUP BY sex_grp
UNION ALL
SELECT 'Antiproliferative agents', sex_grp,
       SUM(CASE WHEN exposed_ANTI = 1 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_ANTI = 1 AND diki = 0 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_ANTI = 0 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_ANTI = 0 AND diki = 0 THEN 1 ELSE 0 END)
FROM base GROUP BY sex_grp
UNION ALL
SELECT 'Corticosteroids', sex_grp,
       SUM(CASE WHEN exposed_CS   = 1 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_CS   = 1 AND diki = 0 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_CS   = 0 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_CS   = 0 AND diki = 0 THEN 1 ELSE 0 END)
FROM base GROUP BY sex_grp
UNION ALL
SELECT 'Biologics', sex_grp,
       SUM(CASE WHEN exposed_BIO  = 1 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_BIO  = 1 AND diki = 0 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_BIO  = 0 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_BIO  = 0 AND diki = 0 THEN 1 ELSE 0 END)
FROM base GROUP BY sex_grp;
"
dt_sex <- as.data.table(DBI::dbGetQuery(con, sql_sex))
res_sex <- dt_sex[, cbind(.SD, calc_metrics(a,b,c,d)), by = .(exposure, sex_grp, a,b,c,d)]
res_sex <- finalize_and_write(res_sex,
                              file.path(OUT_DIR, "Signals_DIKI_bySex_byDrug.csv"),
                              subgroup_col = "sex_grp",
                              subgroup_levels = c("F","M","Other"))

# ------------------------------- 2) Age groups -------------------------------
sql_age <- "
WITH base AS (SELECT * FROM case_diki_flag WHERE age_group IS NOT NULL)
SELECT 'Calcineurin inhibitors' AS exposure, age_group,
       SUM(CASE WHEN exposed_CNI  = 1 AND diki = 1 THEN 1 ELSE 0 END) AS a,
       SUM(CASE WHEN exposed_CNI  = 1 AND diki = 0 THEN 1 ELSE 0 END) AS b,
       SUM(CASE WHEN exposed_CNI  = 0 AND diki = 1 THEN 1 ELSE 0 END) AS c,
       SUM(CASE WHEN exposed_CNI  = 0 AND diki = 0 THEN 1 ELSE 0 END) AS d
FROM base GROUP BY age_group
UNION ALL
SELECT 'mTOR inhibitors', age_group,
       SUM(CASE WHEN exposed_mTOR = 1 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_mTOR = 1 AND diki = 0 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_mTOR = 0 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_mTOR = 0 AND diki = 0 THEN 1 ELSE 0 END)
FROM base GROUP BY age_group
UNION ALL
SELECT 'Antiproliferative agents', age_group,
       SUM(CASE WHEN exposed_ANTI = 1 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_ANTI = 1 AND diki = 0 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_ANTI = 0 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_ANTI = 0 AND diki = 0 THEN 1 ELSE 0 END)
FROM base GROUP BY age_group
UNION ALL
SELECT 'Corticosteroids', age_group,
       SUM(CASE WHEN exposed_CS   = 1 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_CS   = 1 AND diki = 0 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_CS   = 0 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_CS   = 0 AND diki = 0 THEN 1 ELSE 0 END)
FROM base GROUP BY age_group
UNION ALL
SELECT 'Biologics', age_group,
       SUM(CASE WHEN exposed_BIO  = 1 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_BIO  = 1 AND diki = 0 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_BIO  = 0 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_BIO  = 0 AND diki = 0 THEN 1 ELSE 0 END)
FROM base GROUP BY age_group;
"
dt_age <- as.data.table(DBI::dbGetQuery(con, sql_age))
res_age <- dt_age[, cbind(.SD, calc_metrics(a,b,c,d)), by = .(exposure, age_group, a,b,c,d)]
res_age <- finalize_and_write(res_age,
                              file.path(OUT_DIR, "Signals_DIKI_byAgeGroup_byDrug.csv"),
                              subgroup_col = "age_group",
                              subgroup_levels = c("<18","18-39","40-64","65+"))

# ------------------------------- 3) Region (5) -------------------------------
sql_region <- "
WITH base AS (SELECT * FROM case_diki_flag)
SELECT 'Calcineurin inhibitors' AS exposure, region_grp,
       SUM(CASE WHEN exposed_CNI  = 1 AND diki = 1 THEN 1 ELSE 0 END) AS a,
       SUM(CASE WHEN exposed_CNI  = 1 AND diki = 0 THEN 1 ELSE 0 END) AS b,
       SUM(CASE WHEN exposed_CNI  = 0 AND diki = 1 THEN 1 ELSE 0 END) AS c,
       SUM(CASE WHEN exposed_CNI  = 0 AND diki = 0 THEN 1 ELSE 0 END) AS d
FROM base GROUP BY region_grp
UNION ALL
SELECT 'mTOR inhibitors', region_grp,
       SUM(CASE WHEN exposed_mTOR = 1 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_mTOR = 1 AND diki = 0 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_mTOR = 0 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_mTOR = 0 AND diki = 0 THEN 1 ELSE 0 END)
FROM base GROUP BY region_grp
UNION ALL
SELECT 'Antiproliferative agents', region_grp,
       SUM(CASE WHEN exposed_ANTI = 1 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_ANTI = 1 AND diki = 0 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_ANTI = 0 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_ANTI = 0 AND diki = 0 THEN 1 ELSE 0 END)
FROM base GROUP BY region_grp
UNION ALL
SELECT 'Corticosteroids', region_grp,
       SUM(CASE WHEN exposed_CS   = 1 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_CS   = 1 AND diki = 0 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_CS   = 0 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_CS   = 0 AND diki = 0 THEN 1 ELSE 0 END)
FROM base GROUP BY region_grp
UNION ALL
SELECT 'Biologics', region_grp,
       SUM(CASE WHEN exposed_BIO  = 1 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_BIO  = 1 AND diki = 0 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_BIO  = 0 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_BIO  = 0 AND diki = 0 THEN 1 ELSE 0 END)
FROM base GROUP BY region_grp;
"
dt_region <- as.data.table(DBI::dbGetQuery(con, sql_region))
res_region <- dt_region[, cbind(.SD, calc_metrics(a,b,c,d)), by = .(exposure, region_grp, a,b,c,d)]
res_region <- finalize_and_write(res_region,
                                 file.path(OUT_DIR, "Signals_DIKI_byRegion5_byDrug.csv"),
                                 subgroup_col = "region_grp",
                                 subgroup_levels = c("US","EUROPE","JAPAN","CHINA","OTHER"))

# ------------------------------- 4) Indication -------------------------------
sql_indi <- "
WITH base AS (SELECT * FROM case_diki_flag)
SELECT 'Calcineurin inhibitors' AS exposure, indication_grp,
       SUM(CASE WHEN exposed_CNI  = 1 AND diki = 1 THEN 1 ELSE 0 END) AS a,
       SUM(CASE WHEN exposed_CNI  = 1 AND diki = 0 THEN 1 ELSE 0 END) AS b,
       SUM(CASE WHEN exposed_CNI  = 0 AND diki = 1 THEN 1 ELSE 0 END) AS c,
       SUM(CASE WHEN exposed_CNI  = 0 AND diki = 0 THEN 1 ELSE 0 END) AS d
FROM base GROUP BY indication_grp
UNION ALL
SELECT 'mTOR inhibitors', indication_grp,
       SUM(CASE WHEN exposed_mTOR = 1 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_mTOR = 1 AND diki = 0 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_mTOR = 0 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_mTOR = 0 AND diki = 0 THEN 1 ELSE 0 END)
FROM base GROUP BY indication_grp
UNION ALL
SELECT 'Antiproliferative agents', indication_grp,
       SUM(CASE WHEN exposed_ANTI = 1 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_ANTI = 1 AND diki = 0 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_ANTI = 0 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_ANTI = 0 AND diki = 0 THEN 1 ELSE 0 END)
FROM base GROUP BY indication_grp
UNION ALL
SELECT 'Corticosteroids', indication_grp,
       SUM(CASE WHEN exposed_CS   = 1 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_CS   = 1 AND diki = 0 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_CS   = 0 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_CS   = 0 AND diki = 0 THEN 1 ELSE 0 END)
FROM base GROUP BY indication_grp
UNION ALL
SELECT 'Biologics', indication_grp,
       SUM(CASE WHEN exposed_BIO  = 1 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_BIO  = 1 AND diki = 0 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_BIO  = 0 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_BIO  = 0 AND diki = 0 THEN 1 ELSE 0 END)
FROM base GROUP BY indication_grp;
"
dt_indi <- as.data.table(DBI::dbGetQuery(con, sql_indi))
res_indi <- dt_indi[, cbind(.SD, calc_metrics(a,b,c,d)), by = .(exposure, indication_grp, a,b,c,d)]
res_indi <- finalize_and_write(res_indi,
                               file.path(OUT_DIR, "Signals_DIKI_byIndication_byDrug.csv"),
                               subgroup_col = "indication_grp",
                               subgroup_levels = c("Transplantation","Autoimmune diseases","Other"))

# --------------------------- 5) Seriousness (2) ------------------------------
# (a) Indication × Seriousness grid (as you coded) → more granular table
sql_indi_ser <- "
WITH base AS (SELECT * FROM case_diki_flag)
SELECT 'Calcineurin inhibitors' AS exposure, indication_grp, serious_grp,
       SUM(CASE WHEN exposed_CNI  = 1 AND diki = 1 THEN 1 ELSE 0 END) AS a,
       SUM(CASE WHEN exposed_CNI  = 1 AND diki = 0 THEN 1 ELSE 0 END) AS b,
       SUM(CASE WHEN exposed_CNI  = 0 AND diki = 1 THEN 1 ELSE 0 END) AS c,
       SUM(CASE WHEN exposed_CNI  = 0 AND diki = 0 THEN 1 ELSE 0 END) AS d
FROM base GROUP BY indication_grp, serious_grp
UNION ALL
SELECT 'mTOR inhibitors', indication_grp, serious_grp,
       SUM(CASE WHEN exposed_mTOR = 1 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_mTOR = 1 AND diki = 0 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_mTOR = 0 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_mTOR = 0 AND diki = 0 THEN 1 ELSE 0 END)
FROM base GROUP BY indication_grp, serious_grp
UNION ALL
SELECT 'Antiproliferative agents', indication_grp, serious_grp,
       SUM(CASE WHEN exposed_ANTI = 1 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_ANTI = 1 AND diki = 0 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_ANTI = 0 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_ANTI = 0 AND diki = 0 THEN 1 ELSE 0 END)
FROM base GROUP BY indication_grp, serious_grp
UNION ALL
SELECT 'Corticosteroids', indication_grp, serious_grp,
       SUM(CASE WHEN exposed_CS   = 1 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_CS   = 1 AND diki = 0 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_CS   = 0 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_CS   = 0 AND diki = 0 THEN 1 ELSE 0 END)
FROM base GROUP BY indication_grp, serious_grp
UNION ALL
SELECT 'Biologics', indication_grp, serious_grp,
       SUM(CASE WHEN exposed_BIO  = 1 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_BIO  = 1 AND diki = 0 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_BIO  = 0 AND diki = 1 THEN 1 ELSE 0 END),
       SUM(CASE WHEN exposed_BIO  = 0 AND diki = 0 THEN 1 ELSE 0 END)
FROM base GROUP BY indication_grp, serious_grp;
"
dt_indi_ser <- as.data.table(DBI::dbGetQuery(con, sql_indi_ser))
res_indi_ser <- dt_indi_ser[, cbind(.SD, calc_metrics(a,b,c,d)), by=.(exposure, indication_grp, serious_grp, a,b,c,d)]
res_indi_ser <- finalize_and_write(res_indi_ser,
                                   file.path(OUT_DIR, "Signals_DIKI_byIndication_bySerious_byDrug.csv"),
                                   subgroup_col = "serious_grp",
                                   subgroup_levels = c("Serious","Non-serious"))

# (b) Collapse over indication → Serious vs Non-serious by class
res_ser_sum <- res_indi_ser[, .(a = sum(a), b = sum(b), c = sum(c), d = sum(d)), by = .(exposure, serious_grp)]
res_ser_fin <- res_ser_sum[, cbind(.SD, calc_metrics(a,b,c,d)), by=.(exposure, serious_grp, a,b,c,d)]
res_ser_fin <- finalize_and_write(res_ser_fin,
                                  file.path(OUT_DIR, "Signals_DIKI_bySerious_byDrug.csv"),
                                  subgroup_col = "serious_grp",
                                  subgroup_levels = c("Serious","Non-serious"))
# --------------------------- 6) TTO ------------------------------
suppressPackageStartupMessages({
  library(DBI); library(duckdb); library(data.table)
})

FAERS_FILE <- "D:/FAERS/MASTER/FAERS_MASTER_FILE_2004-2024_with_TTO_ALL.parquet"
DIKI_PT_CSV <- "D:/FAERS/MASTER/diki_pt_v28.0_final.csv"
OUT_DIR <- "D:/FAERS/Inhibitors"; if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, TRUE)

# ---- DuckDB connection (read-only parquet; temp DB is fine) ----
con <- dbConnect(duckdb::duckdb(), dbdir=":memory:", read_only=TRUE)
DBI::dbExecute(con, "PRAGMA threads=4;")
DBI::dbExecute(con, "PRAGMA memory_limit='8GB';")
DBI::dbExecute(con, "PRAGMA temp_directory='D:/duck_tmp';")

norm_path <- function(p) gsub("\\\\", "/", p)

# ==== 1) Base view: map your actual columns ====
# We prefer TTO computed on PS+SS pairs; fallback to *_ps then *_any if needed.
DBI::dbExecute(con, sprintf("
  CREATE OR REPLACE VIEW base_raw AS
  SELECT
    primaryid,
    UPPER(COALESCE(drugname_u,''))   AS drugname_std,
    UPPER(COALESCE(role_cod_u,''))   AS role_cod_std,
    UPPER(COALESCE(pt_u,''))         AS pt_std,
    meddra_pt_code                   AS pt_code,
    /* choose TTO in priority: psss -> ps -> any */
    COALESCE(tto_onset_days_psss,
             tto_onset_days_ps,
             tto_onset_days_any)     AS tto_days
  FROM read_parquet('%s')
", norm_path(FAERS_FILE)))
DBI::dbGetQuery(con, "
  SELECT 
    SUM(flag_cni)  AS cni,
    SUM(flag_mtor) AS mtor,
    SUM(flag_anti) AS anti_n,
    SUM(flag_cs)   AS cs,
    SUM(flag_bio)  AS bio,
    COUNT(*)       AS n
  FROM exposure_by_case;
")
# Quick check
print(DBI::dbGetQuery(con, "SELECT * FROM base_raw LIMIT 5"))

# ==== 2) DIKI PT list (join by code or PT name) ====
# Expecting CSV to have columns like: pt_code / pt / etc. We'll normalize.
diki <- fread(DIKI_PT_CSV)
# Standardize helper columns
if (!'pt_code' %in% names(diki)) {
  # try to find a code-like column
  cand <- grep("code|pt_code|meddra", names(diki), ignore.case=TRUE, value=TRUE)
  if (length(cand)) setnames(diki, cand[1], "pt_code")
}
if (!'pt' %in% names(diki)) {
  cand <- grep("pt$|pt_name|preferred", names(diki), ignore.case=TRUE, value=TRUE)
  if (length(cand)) setnames(diki, cand[1], "pt")
}
if ("pt" %in% names(diki)) diki[, pt_upper := toupper(trimws(pt))]
if (!"pt_code" %in% names(diki)) diki[, pt_code := NA_integer_]

# Register into DuckDB as temp table
DBI::dbWriteTable(con, "diki_ptlist_tmp", diki, temporary = TRUE, overwrite = TRUE)

# Build unified DIKI lookup with both code and name (where available)
DBI::dbExecute(con, "
  CREATE OR REPLACE VIEW diki_lookup AS
  SELECT DISTINCT
    CAST(pt_code AS BIGINT) AS pt_code_norm,
    UPPER(COALESCE(pt_upper,'')) AS pt_upper_norm
  FROM diki_ptlist_tmp
")

# ==== 3) Attach DIKI flag to rows (prefer code match; fallback name) ====
DBI::dbExecute(con, "
  CREATE OR REPLACE VIEW base_with_diki AS
  WITH by_code AS (
    SELECT b.*, CASE WHEN dl.pt_code_norm IS NOT NULL THEN 1 ELSE 0 END AS is_diki_by_code
    FROM base_raw b
    LEFT JOIN diki_lookup dl
      ON b.pt_code = dl.pt_code_norm
  ),
  final AS (
    SELECT bc.*,
           CASE
             WHEN is_diki_by_code = 1 THEN 1
             WHEN dl2.pt_upper_norm IS NOT NULL THEN 1
             ELSE 0
           END AS flag_diki
    FROM by_code bc
    LEFT JOIN diki_lookup dl2
      ON bc.pt_std = dl2.pt_upper_norm
  )
  SELECT * FROM final
")

# QC: how many DIKI
print(DBI::dbGetQuery(con, "SELECT SUM(flag_diki) AS n_diki, COUNT(*) AS n_all FROM base_with_diki"))

# ==== 4) Restrict to PS/SS drug roles and build exposures (CNI/mTOR/Anti/CS/Bio) ====
DBI::dbExecute(con, "CREATE OR REPLACE VIEW psss_rows AS
  SELECT primaryid, drugname_std, role_cod_std, tto_days, flag_diki
  FROM base_with_diki
  WHERE role_cod_std IN ('PS','SS')
")

# Dictionaries (can extend as needed)
CNI  <- c('TACROLIMUS','FK506','PROGRAF','ADVAGRAF','CYCLOSPORINE','NEORAL','SANDIMMUNE')
MTOR <- c('EVEROLIMUS','AFINITOR','CERTICAN','SIROLIMUS','RAPAMUNE','TEMSIROLIMUS','TORISEL')
ANTI <- c('MYCOPHENOLATE MOFETIL','CELLCEPT','MYCOPHENOLIC ACID','MYFORTIC','AZATHIOPRINE','IMURAN')
CS   <- c('PREDNISONE','DELTASONE','METHYLPREDNISOLONE','MEDROL','SOLU-MEDROL','DEXAMETHASONE','DECADRON','PREDNISOLONE')
BIO  <- c('RITUXIMAB','MABTHERA','RITUXAN','BASILIXIMAB','SIMULECT','BELATACEPT','NULOJIX')

like_or <- function(vec, col="drugname_std") paste(sprintf("%s LIKE '%%%s%%'", col, gsub("'", "''", vec)), collapse=" OR ")

DBI::dbExecute(con, sprintf("
  CREATE OR REPLACE VIEW exposure_by_case AS
  SELECT
    primaryid,
    MAX(CASE WHEN %s THEN 1 ELSE 0 END) AS flag_cni,
    MAX(CASE WHEN %s THEN 1 ELSE 0 END) AS flag_mtor,
    MAX(CASE WHEN %s THEN 1 ELSE 0 END) AS flag_anti,
    MAX(CASE WHEN %s THEN 1 ELSE 0 END) AS flag_cs,
    MAX(CASE WHEN %s THEN 1 ELSE 0 END) AS flag_bio,
    ANY_VALUE(tto_days) AS tto_days,
    MAX(flag_diki) AS flag_diki
  FROM psss_rows
  GROUP BY primaryid
",
                            like_or(CNI), like_or(MTOR), like_or(ANTI), like_or(CS), like_or(BIO)
))

DBI::dbGetQuery(con, "
  SELECT 
    SUM(flag_cni)  AS cni,
    SUM(flag_mtor) AS mtor,
    SUM(flag_anti) AS anti_n,
    SUM(flag_cs)   AS cs,
    SUM(flag_bio)  AS bio,
    COUNT(*)       AS n
  FROM exposure_by_case;
")


# ==== 5) TTO bins ====
DBI::dbExecute(con, "
  CREATE OR REPLACE VIEW exposure_tto AS
  SELECT *,
         CASE
           WHEN tto_days IS NULL THEN 'Missing'
           WHEN tto_days <= 7   THEN '0-7d'
           WHEN tto_days <= 30  THEN '8-30d'
           WHEN tto_days <= 90  THEN '31-90d'
           WHEN tto_days <= 365 THEN '91-365d'
           ELSE '>365d'
         END AS tto_bin
  FROM exposure_by_case
")

print(DBI::dbGetQuery(con, "SELECT tto_bin, COUNT(*) n FROM exposure_tto GROUP BY tto_bin ORDER BY n DESC"))

# ==== 6) Class × TTO table  ====
sql_summary <- "
  WITH labeled AS (
    SELECT
      CASE
        WHEN flag_cni = 1  THEN 'Calcineurin inhibitors'
        WHEN flag_mtor = 1 THEN 'mTOR inhibitors'
        WHEN flag_anti = 1 THEN 'Antiproliferative agents'
        WHEN flag_cs = 1   THEN 'Corticosteroids'
        WHEN flag_bio = 1  THEN 'Biologics'
        ELSE NULL
      END AS drug_category,
      tto_bin, flag_diki
    FROM exposure_tto
  )
  SELECT drug_category, tto_bin,
         COUNT(*) AS total_reports,
         SUM(flag_diki) AS diki_cases,
         ROUND(100.0 * SUM(flag_diki)/NULLIF(COUNT(*),0), 2) AS diki_rate_pct
  FROM labeled
  WHERE drug_category IS NOT NULL
  GROUP BY 1,2
  ORDER BY drug_category, tto_bin
"
tto_summary <- DBI::dbGetQuery(con, sql_summary)
fwrite(tto_summary, file.path(OUT_DIR, "TTO_Distribution_by_DrugCategory.csv"))

print(head(tto_summary, 10))
dbDisconnect(con)
cat("\nSaved: ", file.path(OUT_DIR, "TTO_Distribution_by_DrugCategory.csv"), "\n")

# ------------------------------- Wrap up -------------------------------------
try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
message("All subgroup analyses completed. Outputs in: ", OUT_DIR)