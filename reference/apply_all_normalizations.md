# Apply all normalization methods (legacy wrapper)

Thin wrapper around
[`apply_all_normalizations_contrast()`](https://byrumlab.github.io/proteoDA/reference/apply_all_normalizations_contrast.md)
to preserve the original function name used in tests and older code.

## Usage

``` r
apply_all_normalizations(data, input_is_log2 = FALSE, groups = NULL)
```

## Arguments

- data:

  Matrix/data.frame of intensities (features x samples).

- input_is_log2:

  Logical, whether `data` are already on log2 scale.

- groups:

  Optional character vector of group labels, used for cycloess.

## Value

A named list of matrices as returned by
[`apply_all_normalizations_contrast()`](https://byrumlab.github.io/proteoDA/reference/apply_all_normalizations_contrast.md).
