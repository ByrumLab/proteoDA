# Evaluate a normalization metric across a list of normalized data

Applies a metric (see
[`norm_metrics`](https://byrumlab.github.io/proteoDA/reference/norm_metrics.md))
across a list of normalized data matrices (most likely produced by
[`apply_all_normalizations`](https://byrumlab.github.io/proteoDA/reference/apply_all_normalizations.md))
and outputs the results in a data frame ready to be used with internal
plotting functions (see
[`pn_plots_generic`](https://byrumlab.github.io/proteoDA/reference/pn_plots_generic.md)).

## Usage

``` r
eval_pn_metric_for_plot(
  normList,
  grouping,
  metric = c("PCV", "PMAD", "PEV", "COR", "log2ratio")
)
```

## Arguments

- normList:

  A named list of normalized data matrices.

- grouping:

  A character or factor vector, listing the group(s) the samples belong
  to.

- metric:

  The normalization metric to calculate. Can be "PCV", "PMAD", "PEV", or
  "COR". See
  [`norm_metrics`](https://byrumlab.github.io/proteoDA/reference/norm_metrics.md).

## Value

A data frame containing the selected metric, ready to be used for
plotting.
