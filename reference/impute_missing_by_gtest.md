# Run G-test-based imputation on DAList data or contrast-specific data

Run G-test-based imputation on DAList data or contrast-specific data

## Usage

``` r
impute_missing_by_gtest(
  DAList,
  contrast = NULL,
  grouping_column = "group",
  p_threshold = 0.01
)
```

## Arguments

- DAList:

  A DAList object after filtering.

- contrast:

  Optional character string specifying contrast name (for use with
  data_per_contrast). If NULL, loops through all.

- grouping_column:

  Name of column in metadata to define group membership (default =
  "group").

- p_threshold:

  P-value threshold for triggering imputation (default = 0.05).

## Value

The same DAList with imputed data written into \$data or
\$data_per_contrast.
