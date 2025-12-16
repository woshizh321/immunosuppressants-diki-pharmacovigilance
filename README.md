# Immunosuppressant-associated DIKI Pharmacovigilance
# Immunosuppressant-associated Drug-Induced Kidney Injury (DIKI)
### A pharmacovigilance study based on FAERS and JADER databases

## Overview
This repository contains the complete analysis pipeline for a pharmacovigilance study
investigating drug-induced kidney injury (DIKI) associated with major immunosuppressant
classes and molecules.

The study integrates large-scale spontaneous reporting data from:
- **FAERS (FDA Adverse Event Reporting System)**
- **JADER (Japanese Adverse Drug Event Report database)**

and follows a figure-driven, fully reproducible workflow corresponding to the main figures
of the manuscript.

---

## Analysis Pipeline (Figure-oriented)

### **Figure 1. Overall characteristics**
- Temporal trends of reporting volume
- Demographic distributions (age, sex, region)
- Drug-class–level exposure patterns

📁 `01_overall_characteristics/`

---

### **Figure 2. Hierarchical signal detection**
- Disproportionality analyses at SOC, HLGT, HLT, and PT levels
- Metrics: ROR, PRR, χ²

📁 `02_hierarchical_signal_detection/`

---

### **Figure 3. Stratified subgroup analyses**
- Stratification by sex, age, region, indication, seriousness, and time-to-onset (TTO)

📁 `03_subgroup_analysis/`

---

### **Figure 4. Drug-level LASSO modeling**
- L1-penalized logistic regression at drug-class level
- Excluding region and TTO
- Model performance and interpretability analyses

📁 `04_drug_level_lasso/`

---

### **Figure 5. Molecule-level LASSO modeling**
- Fine-grained molecular predictors within CNI and APA classes
- SHAP-like contribution analysis
- Internal cross-validation

📁 `05_molecule_level_lasso/`

---

### **Figure 6. Individualized risk portraits**
- Integration of molecular, demographic, and clinical context
- Absolute and relative DIKI risk estimation

📁 `06_external_validation_JADER/`

---

### **Figure 7. External validation and TTO integration**
- External validation using JADER
- Reintroduction of time-to-onset into risk modeling

📁 `07_TTO_integration/`

---

## Data availability
Due to data-use agreements, raw FAERS and JADER datasets are **not publicly redistributed**.
All analyses assume access to locally stored, preprocessed MASTER files.

---

## Reproducibility
- All scripts are written in **R (≥4.2)** with `duckdb`, `data.table`, and `arrow`
- Each analysis module can be executed independently
- Random seeds are fixed where applicable

---

## Disclaimer
This study is based on spontaneous reporting data and does not establish causality.
Findings should be interpreted as hypothesis-generating.
