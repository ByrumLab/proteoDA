# Plot histograms highlighting imputed values (Perseus-style imputation)

If you ran `perseus_impute(DAList, save_before_after=TRUE)`, you can
pass:
`logDat_before = DAList$data_per_contrast[[contrast]]$imputation$before_log2`
`logDat_after = DAList$data_per_contrast[[contrast]]$imputation$after_log2`
or use the global `DAList$imputation$before_log2/after_log2` when no
per-contrast data exists.

## Usage

``` r
plot_perseus_imputation(
  logDat_before,
  logDat_after,
  samples = NULL,
  bins = 50,
  facet_ncol = 4,
  overlay = TRUE
)
```

## Arguments

- logDat_before:

  matrix/data.frame of log2 intensities with NAs (pre-imputation)

- logDat_after:

  matrix/data.frame of log2 intensities after
  [`perseus_impute()`](https://byrumlab.github.io/proteoDA/reference/perseus_impute.md)

- samples:

  character/integer vector of columns to plot (default: all)

- bins:

  number of histogram bins

- facet_ncol:

  number of columns in facet wrap

- overlay:

  if TRUE, overlay observed+imputed per sample panel; else stacked

## Value

ggplot object (faceted histograms across selected samples)
