#!/usr/bin/env Rscript
## ============================================================================
## 06_fig6_generate_risk_portraits.R
##
## Objective:
##    - Construct a population risk portrait table based on the Refit OR results
##      of the "molecule-level + indication + noYear/noRegion/noTTO/noSerious" model.
##
##    - Portrait Dimensions:
##        * Molecular Exposure: None / Tacrolimus / Cyclosporine / Everolimus / Sirolimus / Temsirolimus
##        * Sex: Female / Male / Unknown (optional)
##        * Age: <18 / 18–39 / 40–64 / 65+
##        * Indication: Other / Transplantation / Autoimmune diseases
##
##    - Output:
##        * portrait_risk_table_molecule_withIndi_noyear_lambdaMIN.csv
##
## Assumptions:
##    - Refit OR file originates from the lambda.min model
##    - Factor Baselines:
##        sex baseline        = Female (no explicit term)
##        age_cat baseline    = "<18" (no explicit term)
##        indication baseline = "Other" (no explicit term)
##        molecule baseline   = No molecular flag set to 1
##
## Paths:
##    - DATA_PATH = "D:/FAERS/Inhibitors"
##    - OR_FILE   = file.path(OUT_MAIN, "lasso_refit_logistic_results_molecule_withIndication_noRegion_noTTO_noSerious_noyear_lambdaMIN.csv")
##    - OUT_MAIN  = "D:/FAERS/Inhibitors/Derived_molecule_withIndi_noyear"
## ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

## ---------------------------- Configuration ----------------------------

DATA_PATH <- "D:/FAERS/Inhibitors"
OUT_MAIN  <- file.path(DATA_PATH, "Derived_molecule_withIndi_noyear")

OR_FILE   <- file.path(
  OUT_MAIN,
  "lasso_refit_logistic_results_molecule_withIndication_noRegion_noTTO_noSerious_noyear_lambdaMIN.csv"
)

OUT_FILE  <- file.path(
  OUT_MAIN,
  "portrait_risk_table_molecule_withIndi_noyear_lambdaMIN.csv"
)

## Risk Level Thresholds (Absolute Risk) — Adjust as needed
RISK_BREAKS <- c(0, 0.02, 0.05, 1)  # <2% Low, 2–5% Moderate, >5% High
RISK_LABELS <- c("Low", "Moderate", "High")

## Whether to include sex 'Unknown' in portraits (TRUE = Generate portraits for Unknown sex)
INCLUDE_SEX_UNKNOWN <- TRUE

## ------------------------- Step 1. Input & QC -------------------

stopifnot(file.exists(OR_FILE))

or_dt <- fread(OR_FILE)
cat("Preview of Refit OR results:\n")
print(head(or_dt, 10))

## Required columns check
req_cols <- c("term", "OR", "CI_low", "CI_high", "p_value")
if (!all(req_cols %in% names(or_dt))) {
  stop("Refit OR file is missing required columns: ", paste(setdiff(req_cols, names(or_dt)), collapse = ", "))
}

## Remove useless items like id_case
or_dt <- or_dt[term != "id_case"]

## Extract Intercept
if (!"(Intercept)" %in% or_dt$term) {
  stop("(Intercept) not found in Refit OR file.")
}
or_intercept <- or_dt[term == "(Intercept)", OR][1]
beta0        <- log(or_intercept)
base_risk    <- or_intercept / (1 + or_intercept)

cat(sprintf("\n>>> Intercept OR = %.6f, corresponding Baseline Risk ≈ %.3f%%\n",
            or_intercept, base_risk * 100))

## Remaining Covariates
coef_dt <- or_dt[term != "(Intercept)"]
coef_dt[, beta := log(OR)]

## -------------------- Step 2. Identify Terms -------------------

## Molecule term names (consistent with LASSO)
mol_terms <- c(
  "flag_mol_tacrolimus",
  "flag_mol_cyclosporine",
  "flag_mol_everolimus",
  "flag_mol_sirolimus",
  "flag_mol_temsirolimus"
)

## Check which molecule terms are actually present
mol_terms_present <- intersect(mol_terms, coef_dt$term)
cat("\nMolecule terms present:\n")
print(mol_terms_present)

## Sex terms
sex_terms <- coef_dt[grep("^sex", term), term]
cat("\nSex-related terms:\n")
print(sex_terms)

## Age terms
age_terms <- coef_dt[grep("^age_cat", term), term]
cat("\nAge-related terms:\n")
print(age_terms)

## Indication terms
indi_terms <- coef_dt[grep("^indication_grp", term), term]
cat("\nIndication-related terms:\n")
print(indi_terms)

## Sanity check: verify existence of core factors
needed_terms <- c("flag_mol_cyclosporine",
                  "flag_mol_everolimus",
                  "flag_mol_temsirolimus",
                  "sexMale",
                  "age_cat18.39",
                  "age_cat40.64",
                  "age_cat65.",
                  "indication_grpTransplantation",
                  "indication_grpAutoimmune.diseases")
missing_needed <- setdiff(needed_terms, coef_dt$term)
if (length(missing_needed) > 0) {
  cat("Warning: The following expected terms were not found in Refit OR results and will be treated as 0 (no effect) in portraits:\n")
  print(missing_needed)
}

## Create a named vector beta_vec for convenience
beta_vec <- coef_dt$beta
names(beta_vec) <- coef_dt$term

## Helper function to safely retrieve beta (returns 0 for missing terms)
get_beta <- function(term_name) {
  if (is.na(term_name) || term_name == "") return(0)
  if (!(term_name %in% names(beta_vec))) return(0)
  beta_vec[[term_name]]
}

## -------------------- Step 3. Define Portrait Combination Space -----------------

## 1) Molecule Patterns (can be expanded: currently set as single-agent exposure + None)
##    pattern_id is a label used to describe the portrait
mol_patterns <- list(
  "None"         = character(0),
  "Tacrolimus"   = "flag_mol_tacrolimus",
  "Cyclosporine" = "flag_mol_cyclosporine",
  "Everolimus"   = "flag_mol_everolimus",
  "Sirolimus"    = "flag_mol_sirolimus",
  "Temsirolimus" = "flag_mol_temsirolimus"
)

## If a molecule term does not exist in the model:
## The script still generates the portrait with beta=0 (baseline effect).

## 2) Sex Level mapping to terms
## baseline (no term) = Female
sex_levels <- c("Female", "Male")
if (INCLUDE_SEX_UNKNOWN && "sexUnknown" %in% coef_dt$term) {
  sex_levels <- c(sex_levels, "Unknown")
}

map_sex_to_term <- function(sex_label) {
  if (sex_label == "Male") {
    return(if ("sexMale" %in% coef_dt$term) "sexMale" else NA_character_)
  } else if (sex_label == "Unknown") {
    return(if ("sexUnknown" %in% coef_dt$term) "sexUnknown" else NA_character_)
  } else {
    ## Female = baseline
    return(NA_character_)
  }
}

## 3) Age Stratification mapping to terms
## baseline = "<18"
age_levels <- c("<18", "18-39", "40-64", "65+")

map_age_to_term <- function(age_label) {
  if (age_label == "18-39") {
    return(if ("age_cat18.39" %in% coef_dt$term) "age_cat18.39" else NA_character_)
  } else if (age_label == "40-64") {
    return(if ("age_cat40.64" %in% coef_dt$term) "age_cat40.64" else NA_character_)
  } else if (age_label == "65+") {
    return(if ("age_cat65." %in% coef_dt$term) "age_cat65." else NA_character_)
  } else {
    return(NA_character_)  # "<18" baseline
  }
}

## 4) Indication Stratification mapping to terms
## baseline = "Other"
indi_levels <- c("Other", "Transplantation", "Autoimmune diseases")

map_indi_to_term <- function(indi_label) {
  if (indi_label == "Transplantation") {
    return(if ("indication_grpTransplantation" %in% coef_dt$term)
      "indication_grpTransplantation" else NA_character_)
  } else if (indi_label == "Autoimmune diseases") {
    return(if ("indication_grpAutoimmune.diseases" %in% coef_dt$term)
      "indication_grpAutoimmune.diseases" else NA_character_)
  } else {
    return(NA_character_)  # "Other" baseline
  }
}

## -------------------- Step 4. Build Portrait Grid and Calculate Risk -------------

cat("\n>>> Building portrait combinations and calculating logOR / OR / Absolute Risk...\n")

portrait_list <- list()

for (mol_name in names(mol_patterns)) {
  mol_terms_vec <- mol_patterns[[mol_name]]
  
  for (sex_lbl in sex_levels) {
    sex_term <- map_sex_to_term(sex_lbl)
    
    for (age_lbl in age_levels) {
      age_term <- map_age_to_term(age_lbl)
      
      for (indi_lbl in indi_levels) {
        indi_term <- map_indi_to_term(indi_lbl)
        
        ## Collect all terms involved in this specific portrait
        terms_this <- c(mol_terms_vec, sex_term, age_term, indi_term)
        terms_this <- terms_this[!is.na(terms_this) & nzchar(terms_this)]
        
        ## Sum the betas
        beta_sum <- 0
        if (length(terms_this) > 0) {
          beta_sum <- sum(vapply(terms_this, get_beta, numeric(1)))
        }
        
        logOR_combo  <- beta_sum
        OR_combo     <- exp(logOR_combo)
        
        ## Calculate Absolute Risk: logit(p) = beta0 + beta_sum
        lin_pred     <- beta0 + beta_sum
        risk_abs     <- exp(lin_pred) / (1 + exp(lin_pred))
        
        portrait_list[[length(portrait_list) + 1L]] <- data.table(
          sex        = sex_lbl,
          age        = age_lbl,
          indication = indi_lbl,
          molecule   = mol_name,
          terms_used = paste(terms_this, collapse = ";"),
          logOR_combo = logOR_combo,
          OR_combo    = OR_combo,
          risk_abs    = risk_abs
        )
      }
    }
  }
}

portrait_dt <- rbindlist(portrait_list)

## Assign Risk Categories
portrait_dt[, risk_cat := cut(
  risk_abs,
  breaks = RISK_BREAKS,
  labels = RISK_LABELS,
  include.lowest = TRUE, right = FALSE
)]

## Sort by absolute risk (highest to lowest) for better visualization
setorder(portrait_dt, -risk_abs)

cat("\nPortrait Table Preview (Top 20 by descending absolute risk):\n")
print(head(portrait_dt, 20))

## Save to CSV
fwrite(
  portrait_dt,
  OUT_FILE
)

cat("\n>>> Population risk portraits have been exported to:\n  ", OUT_FILE, "\n")