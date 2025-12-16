# Figure 4 – Drug-level LASSO regression

## Objective
Identify drug-class–level predictors of DIKI using penalized logistic regression.

## Data
FAERS case-level dataset with drug-class exposure indicators.

## Methods
- L1-penalized logistic regression (LASSO) for variable selection
- Refit logistic regression using selected predictors
- Internal validation via cross-validation and calibration

## Notes
Region and TTO variables were intentionally excluded to
preserve compatibility with external validation.
