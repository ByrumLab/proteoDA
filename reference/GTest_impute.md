# G-test for Detection Imbalance and Conditional Imputation

Performs a G-test (likelihood ratio test) on per-protein detection
counts across two groups based on observed and expected detection
values. Missing values are imputed with the group-wise minimum if either
(1) the p-value is below a defined threshold, or (2) one group has all
missing values.

## Usage

``` r
GTest_impute(
  data_matrix,
  group_labels,
  p_threshold = 0.01,
  return_imputed_flags = TRUE
)
```

## Arguments

- data_matrix:

  A numeric matrix of protein intensities (rows = proteins, columns =
  samples).

- group_labels:

  A character or factor vector specifying the group for each sample
  (e.g., "N" or "P").

- p_threshold:

  A numeric p-value threshold for applying imputation (default = 0.05).

- return_imputed_flags:

  Logical; if TRUE, return a logical vector indicating which proteins
  were imputed.

## Value

A list with:

- p_values:

  Named numeric vector of G-test p-values per protein.

- imputed_data:

  Matrix of the same dimension as `data_matrix` with imputations
  applied.

- imputed_flags:

  (Optional) Logical vector indicating which proteins were imputed.
