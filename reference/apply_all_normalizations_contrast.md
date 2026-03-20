# Apply all normalization methods to a matrix (contrast-aware helper)

Takes in a matrix of intensities and applies 8 normalization methods,
returning a named list of matrices. `input_is_log2` controls whether to
log2 internally for methods that expect log2 input. `groups` is only
used for `cycloess`.

## Usage

``` r
apply_all_normalizations_contrast(data, input_is_log2 = FALSE, groups = NULL)
```

## Arguments

- data:

  Matrix/data.frame (rows = features, cols = samples).

- input_is_log2:

  Logical, data are already on log2 scale for log2 methods?

- groups:

  Character vector of group labels (length = ncol(data)); used for
  cycloess.

## Value

A named list of 8 matrices.
