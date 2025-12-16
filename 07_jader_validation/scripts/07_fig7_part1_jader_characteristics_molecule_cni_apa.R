#!/usr/bin/env Rscript
## ============================================================
## 07_fig7_part1_jader_characteristics_molecule_cni_apa.R
## Objectives: Generate descriptive tables for CNI + APA 
## (CyA / Tac / MMF / AZA) in the JADER database:
##    1) Sex distribution table
##    2) Age median + IQR
##    3) Age group distribution table
##    4) Top 10 country distribution table
## ============================================================

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(data.table)
})

## ---------------- Global Parameters ----------------
MASTER  <- "D:/JADER/MASTER/JADER_MASTER_PT_English_SUPERMASTER_v4_withDrugCode.parquet"
OUT_DIR <- "D:/JADER/ExternalValidation_molecule_CNI_APA/Descriptives"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

norm_path <- function(x) gsub("\\\\", "/", x)

## Target molecules (standardized to uppercase)
target_molecules <- c(
  "CYCLOSPORINE",
  "TACROLIMUS",
  "MYCOPHENOLATE MOFETIL",
  "AZATHIOPRINE"
)

## ---------------- DuckDB Connection + Data Loading ----------------
con <- dbConnect(duckdb(), read_only = TRUE)
on.exit(try(dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)

sql_load <- sprintf("
  SELECT
      ID AS PRIMARYID,
      UPPER(DRUG_MOLECULE)   AS MOLECULE,
      DRUG_CLASS_STD,
      PT_CODE,
      ROLE_STD,
      SEX,
      AGE_NUM,
      COUNTRY_USE
  FROM read_parquet('%s')
", norm_path(MASTER))

dt <- as.data.table(dbGetQuery(con, sql_load))
cat("Original dt dimensions:", nrow(dt), "rows x", ncol(dt), "cols\n")

## Retain only PS (Primary Suspect) or SS (Secondary Suspect) roles
dt <- dt[ROLE_STD %in% c("PS","SS")]
cat("Dimensions after PS/SS filtering:", nrow(dt), "rows\n")

## Filter for target molecules
dt_mol <- dt[MOLECULE %in% target_molecules]
cat("Dimensions after limiting to 4 CNI+APA molecules:", nrow(dt_mol), "rows\n")
cat("MOLECULE distribution:\n")
print(dt_mol[, .N, by = MOLECULE][order(-N)])

## ============================================================
## Part 1. Sex Distribution Table
## ============================================================

## Standardize sex encoding
dt_mol[, SEX_STD := toupper(trimws(as.character(SEX)))]

dt_mol[SEX_STD %in% c("M","MALE","男性"),   SEX_CAT := "Male"]
dt_mol[SEX_STD %in% c("F","FEMALE","女性"), SEX_CAT := "Female"]
dt_mol[is.na(SEX_CAT),                      SEX_CAT := "Other/Unknown"]

cat("\nSEX_CAT distribution:\n")
print(dt_mol[, .N, by = SEX_CAT])

## molecule × sex counts + proportions
sex_tab <- dt_mol[
  ,
  .(n = .N),
  by = .(MOLECULE, SEX_CAT)
][
  order(MOLECULE, -n)
]

## Calculate proportions within each molecule
sex_tab[, n_total_mol := sum(n), by = MOLECULE]
sex_tab[, pct := n / n_total_mol]

cat("\n=== JADER molecule × sex distribution ===\n")
print(sex_tab)

data.table::fwrite(
  sex_tab,
  file = file.path(OUT_DIR, "JADER_molecule_sex_distribution.csv")
)

## ============================================================
## Part 2. Age Median + IQR
## ============================================================

## AGE_NUM > 0 is considered a valid age
dt_mol[, AGE_NUM := suppressWarnings(as.numeric(AGE_NUM))]
dt_age_valid <- dt_mol[!is.na(AGE_NUM) & AGE_NUM > 0]

cat("\nValid age records count:", nrow(dt_age_valid), "\n")

age_summary <- dt_age_valid[
  ,
  .(
    n_age      = .N,
    median_age = median(AGE_NUM, na.rm = TRUE),
    q1_age       = quantile(AGE_NUM, 0.25, na.rm = TRUE),
    q3_age       = quantile(AGE_NUM, 0.75, na.rm = TRUE)
  ),
  by = MOLECULE
][order(MOLECULE)]

cat("\n=== JADER molecule Age Median + IQR ===\n")
print(age_summary)

data.table::fwrite(
  age_summary,
  file = file.path(OUT_DIR, "JADER_molecule_age_summary.csv")
)

## ============================================================
## Part 3. Age Group Distribution (<18 / 18–39 / 40–64 / 65+)
## ============================================================

dt_age_valid[
  ,
  age_group := cut(
    AGE_NUM,
    breaks = c(-Inf, 18, 40, 65, Inf),
    labels = c("<18","18-39","40-64","65+"),
    right = FALSE
  )
]

age_group_tab <- dt_age_valid[
  !is.na(age_group),
  .(n = .N),
  by = .(MOLECULE, age_group)
][
  order(MOLECULE, age_group)
]

age_group_tab[, n_total_mol := sum(n), by = MOLECULE]
age_group_tab[, pct := n / n_total_mol]

cat("\n=== JADER molecule × age_group distribution ===\n")
print(age_group_tab)

data.table::fwrite(
  age_group_tab,
  file = file.path(OUT_DIR, "JADER_molecule_age_group_distribution.csv")
)

## ============================================================
## Part 4. Country Distribution (Top 10 countries per molecule)
## ============================================================

dt_mol[, COUNTRY_USE := toupper(trimws(as.character(COUNTRY_USE)))]
dt_country_valid <- dt_mol[COUNTRY_USE != "" & !is.na(COUNTRY_USE)]

country_tab <- dt_country_valid[
  ,
  .(n = .N),
  by = .(MOLECULE, COUNTRY_USE)
][
  order(MOLECULE, -n)
]

## Select Top 10 countries for each molecule
country_top10 <- country_tab[, head(.SD, 10), by = MOLECULE]

## Merge with total exposure counts per molecule to calculate proportion
mol_total <- dt_mol[, .(n_total_mol = .N), by = MOLECULE]
country_top10 <- merge(country_top10, mol_total, by = "MOLECULE", all.x = TRUE)
country_top10[, pct := n / n_total_mol]

cat("\n=== JADER molecule × country Top 10 distribution ===\n")
print(country_top10)

data.table::fwrite(
  country_top10,
  file = file.path(OUT_DIR, "JADER_molecule_country_top10.csv")
)

cat("\n✅ JADER molecule-level descriptive tables generated. Output directory:\n", OUT_DIR, "\n")