# Check percentage of DA genes

Internal utility function, used in
[`extract_DA_results`](https://byrumlab.github.io/proteoDA/reference/extract_DA_results.md)
to check if assumptions are met.

## Usage

``` r
check_DA_perc(
  DA_outcomes_table,
  DA_warn_threshold = 0.2,
  pval_thresh,
  lfc_thresh,
  adj_method
)
```

## Arguments

- DA_outcomes_table:

  DA results data frame. Should be the output of
  [`limma::decideTests`](https://rdrr.io/pkg/limma/man/decideTests.html),
  coerced to a data frame.

- DA_warn_threshold:

  Proportion of DA genes at which we warn user.

- pval_thresh:

  P-value threshold used.

- lfc_thresh:

  logFC threshold used.

- adj_method:

  P-value adjustment method used.

## Value

A vector of numeric values giving the % of significant DA proteins
within each contrast.
