# Compute MD data for all normalization methods

Internal helper: builds a long data.frame with mean intensity (A) and
log2 fold change (M) per feature and normalization method.

## Usage

``` r
compute_MD_long(normList, groups)
```

## Arguments

- normList:

  List of normalized expression matrices (rows = features, columns =
  samples), as returned by
  [`apply_all_normalizations_contrast()`](https://byrumlab.github.io/proteoDA/reference/apply_all_normalizations_contrast.md).

- groups:

  Character or factor vector of group labels, length equal to the number
  of columns in each matrix in `normList`.

## Value

A data.frame with columns `A`, `M`, `method`.

## Details

For per-contrast mode, `groups` is already restricted to the two levels
involved in the contrast. If more than two groups are supplied (global
mode), the first two levels are used (with a warning).
