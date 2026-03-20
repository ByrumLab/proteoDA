# Perseus-style MNAR imputation (left-censoring) on log2 data

Applies the common Perseus imputation where, per sample, missing values
are drawn from a narrow normal distribution shifted down from the
observed distribution. Works on a plain log2 matrix/data.frame *or* on a
DAList:

- If `DAList$data_per_contrast` exists and is non-empty, imputes each
  contrast's matrix and stores the result in
  `DAList$data_per_contrast[[contrast]]`.

- Otherwise, imputes the global `DAList$data` matrix in-place.

## Usage

``` r
perseus_impute(
  x,
  shift = 1.8,
  width = 0.3,
  robust = TRUE,
  min_obs_per_sample = 5,
  seed = NULL,
  save_before_after = FALSE,
  store_mask = TRUE
)
```

## Arguments

- x:

  A numeric log2 matrix/data.frame (rows=features, cols=samples), or a
  DAList with `$data` and optionally `$data_per_contrast`.

- shift:

  numeric, how far to shift down the mean in SD units (default 1.8)

- width:

  numeric, fraction of SD used as imputation SD (default 0.3)

- robust:

  logical, use median/MAD (TRUE) or mean/SD (FALSE) per sample

- min_obs_per_sample:

  Integer. Minimum number of observed (non-missing) log2 values required
  in a sample before estimating its own mean (mu_j) and standard
  deviation (sd_j) for imputation. If a sample has fewer than this many
  observed values, the function falls back to using the pooled
  distribution across all samples. Default = 5.

- seed:

  integer or NULL, set for reproducible draws

- save_before_after:

  logical, store full before/after matrices and mask into the DAList (or
  return attrs for matrix mode). Default FALSE.

- store_mask:

  logical, store/log an imputed logical mask. Default TRUE.

## Value

If `x` is a matrix/data.frame, returns the imputed matrix (same dims);
when `store_mask=TRUE`, attaches `attr(result, "imputed_mask")`. If `x`
is a DAList, returns the updated DAList with imputed matrices written to
the appropriate slots; when `save_before_after=TRUE`, persists
`before_log2`, `after_log2`, and `imputed_mask`.

## Details

Optionally stores before/after matrices and an imputed mask for
diagnostics.

For each sample j, let μ_j and σ_j be location and scale of its observed
(log2) values. Missing entries in sample j are drawn from Normal( μ_j -
shift*σ_j, (width*σ_j)^2 ). Features missing in *all* samples remain NA
(uninformative).

This function assumes that the input data are on the **log2** scale. As
a safety check, if the maximum absolute value in a matrix exceeds 1e3, a
warning is issued suggesting that the data may not be log2-transformed.

## Examples

``` r
if (FALSE) { # \dontrun{
# Matrix use
X_imp <- perseus_impute(X, seed=123)
mask  <- attr(X_imp, "imputed_mask")

# DAList use (per-contrast if present)
results <- perseus_impute(results, save_before_after = TRUE, seed=1)
} # }
```
