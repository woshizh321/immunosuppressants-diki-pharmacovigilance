# Figure 2 – Hierarchical signal detection (ROR/PRR/Chi-square)

## Goal
Compute disproportionality signals for DIKI across hierarchical MedDRA levels:
SOC, HLGT, HLT, and PT.

## Data
FAERS MASTER parquet is required locally (not included in this repository).

## Methods
- Construct 2x2 tables (a,b,c,d) for each level.
- Compute ROR, PRR and Chi-square statistics.
- Apply minimum cell count threshold (configurable).

## Outputs
- tables/: signal tables per MedDRA level (include a/b/c/d and effect estimates).
- figures/: manuscript-ready plots.

## Reproducibility notes
- Outputs are generated locally and not tracked by git.
- Configuration file `config_local.R` is intentionally ignored.
