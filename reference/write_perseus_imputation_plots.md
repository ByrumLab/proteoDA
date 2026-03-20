# Write Perseus imputation histograms for each contrast

Uses `DAList$imputation_per_contrast[[ct]]$before_log2` and
`DAList$imputation_per_contrast[[ct]]$after_log2` (created by
`perseus_impute(..., save_before_after = TRUE)`) to generate diagnostic
plots for each contrast.

## Usage

``` r
write_perseus_imputation_plots(
  DAList,
  out_dir = "Imputation_Plots",
  contrasts = NULL,
  samples = NULL,
  bins = 50,
  facet_ncol = 4,
  overlay = TRUE,
  width = 7,
  height = 5,
  dpi = 300,
  device = "png"
)
```

## Arguments

- DAList:

  A DAList with `imputation_per_contrast` diagnostics.

- out_dir:

  Directory to save plots (created if missing). If `NULL`, plots are not
  written to disk.

- contrasts:

  Character vector of contrast names to plot; default is all available
  contrasts in `DAList$imputation_per_contrast`.

- samples:

  Optional subset of sample columns to plot (names or indices); default
  is all samples.

- bins:

  Number of histogram bins.

- facet_ncol:

  Number of columns in the facet layout.

- overlay:

  Logical; if `TRUE`, overlay observed and imputed intensities in the
  same panel; if `FALSE`, use stacked facets.

- width, height:

  Plot size (inches) when saving.

- dpi:

  Resolution (dots per inch) when saving.

- device:

  Graphics device for
  [`ggplot2::ggsave`](https://ggplot2.tidyverse.org/reference/ggsave.html),
  e.g. `"png"`.

## Value

A named list of `ggplot` objects (returned invisibly).

## Details

Each plot is a faceted histogram per sample, comparing observed vs
imputed log2 intensities. The helper
[`plot_perseus_imputation()`](https://byrumlab.github.io/proteoDA/reference/plot_perseus_imputation.md)
is used internally.
