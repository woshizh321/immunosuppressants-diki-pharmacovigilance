# Figure 5 – Molecule-level LASSO regression

## Objective
Identify individual immunosuppressant molecules associated with DIKI
within CNI and antiproliferative agent classes.

## Rationale
Molecule-level effects are expected to be sparse and correlated.
LASSO was applied to reduce overfitting and enhance interpretability.

## Methods
- L1-penalized logistic regression
- Lambda selection by 1-SE rule
- Refit logistic regression on selected molecules

## Notes
Region and TTO were excluded to ensure comparability with
external validation using JADER.
